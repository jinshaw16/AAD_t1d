#imputation_readin.R

library(plyr)
library(dplyr)
library(snpStats)
library(annotSnpStats)
library(jimisc)
library(qqman)
library(humarray)
library(ggplot2)
library(stringr)

#want to have it so each index SNP from onengut is included in our analysis, not just their proxies:

t1dsnps<-c("rs2476601","rs12068671",  "rs6691977", "rs3024505", "rs13415583",
           "rs4849135", "rs2111485", "rs35667974", "rs72871627",
"rs3087243", "rs113010081", "rs6819058","rs67797421", "rs2611215", "rs11954020",
"rs72975913","rs72928038", "rs1538171", "rs62447205", "rs10277986", "rs6476839",
"rs61839660", "rs11594656","rs6602437", "rs41295121",
"rs12416116", "rs689", "rs72853903", "rs917911", "rs705705",
"rs653178", "rs9585056", "rs1456988", "rs56994090", "rs72727394",
"rs34593439", "rs151234", "rs12927355",
"rs193778", "rs8056814", "rs12453507",
"rs757411", "rs1052553", "rs1893217", "rs12971201",
"rs1615504", "rs34536443", "rs12720356", "rs402072", "rs516246",
"rs6043409", "rs11203202", "rs6518350", "rs4820830", "rs229533")

t1dloci<-c("PTPN22", "FASLG","CAMSAP2",
"IL10", "AFF3", "ACOXL",
"IFIH1 (1)", "IFIH1 (2)", "IFIH1 (3)",
"CTLA4", "CCR5", "IL2/IL21 (1)", "IL2/IL21 (2)", "CPE",
"IL7R","PTPRK/THEMIS", "BACH2", "CENPW",
"IKZF1", "COBL",
"GLIS3", "IL2RA (1)", "IL2RA (2)", "IL2RA (3)","IL2RA (4)",
"PTEN","INS (1)","INS (2)", "CD69", "IKZF4",
"SH2B3", "GPR183", "LINC01550",
"MEG3","RASGRP1", "CTSH", "IL27", "DEXI (1)",
"DEXI (2)", "BCAR1","IKZF3",
"CCR7", "MAPT","PTPN2 (1)", "PTPN2 (2)", "CD226", "TYK2 (1)", "TYK2 (2)", "PRKD2",
"FUT2", "SIRPG",
"UBASH3A", "ICOSLG", "HORMAD2",
"C1QTNF6")

t1dsnps<-data.frame(snp=t1dsnps, loci=t1dloci)
t1dsnps$id<-rs.to.id(t1dsnps$snp)
t1dsnps$id<-gsub(".*,","",t1dsnps$id)
t1dsnps$id<-ifelse(substr(t1dsnps$id,1,3)=="seq",gsub("_","-", t1dsnps$id),t1dsnps$id)
t1dsnps$id<-ifelse(substr(t1dsnps$id,1,4)=="X1kg", substr(t1dsnps$id,2,100000),t1dsnps$id)
t1dsnps$ord<-c(1:nrow(t1dsnps))
t1dsnps$snp<-as.character(t1dsnps$snp)


d<-"/well/todd/users/jinshaw/t1d_risk/immunochip/"
#read SNP and phenotype data:
load(file=paste0(d,"all_inds_unrel_postqc.RData"))

p<-t1dsnps[!t1dsnps$id %in% colnames(g),]
p$chromosome=sub("(^.*)[_](.*)[_](.*)","\\2",p$id)

samples<-all@samples

#Read the imputed data in:
d<-"/well/todd/users/jinshaw/aad/under_7/imputation/"

getsnps<-function(snp){
#check allignment with the dataset submitted to impute2:
#read in data from IMPUTE2:
imputed<-read.impute(file=paste0(d,snp,"_n_out"), rownames=samples$V2)
DATA<-imputed
DATA<-DATA[rownames(DATA) %in% rownames(samples),]
cs1<-col.summary(DATA)
cs1$rsid<-gsub(":.*","",rownames(cs1))
cs1$rs_id<-rownames(cs1)

#Check imputation quality:
info<-read.table(file=paste0(d,snp,"_n_out_info"), header=TRUE)
info$snp<-c(1:nrow(info))
cs1<-cbind(cs1,info)
cs1$diff<-abs(cs1$exp_freq_a1-cs1$RAF)

cs2<-cs1[cs1$info>=0.8 & cs1$diff<0.05 & cs1$MAF>0.005,]
p1<-p[p$id==snp,]
chrom=p1$chromosome
a<-grep(p1$snp,cs1$rs_id)
if(length(a)==0){
l<-read.table(file=paste0("/well/1000G/mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/1000GP_Phase3_chr",chrom,".legend.gz"),header=T,as.is=T)
a<-l[grepl(p1$snp,l$id),]
p1$snp<-paste0(chrom,":",a$position,":",a$a1)
}
inf<-info[grepl(p1$snp,info$rs_id),]
cs2<-cs1[grepl(p1$snp, cs1$rs_id),]
if (length(a)>0){
message(paste0(snp,"included with info score ",inf$info,", MAF = ",cs2$MAF,", Certain calls = ",cs2$Certain.calls,", HWE z= ",cs2$z.HWE))
}
if (length(a)==0){
message(paste0(snp,"excluded due to either info score ",inf$info,", MAF = ",cs2$MAF,", Certain calls = ",cs2$Certain.calls,", HWE z= ",cs2$z.HWE))
}
DATA<-DATA[,colnames(DATA) %in% cs2$rs_id]


save(cs2, DATA, file=paste0(d,snp,"_n_clean.RData"))
return(DATA)
}

out<-lapply(p$id,getsnps)


