#impute_missing_snps.R
#There are some SNPs that failed QC - need to impute these, doing this here:

library(snpStats)
library(annotSnpStats)
library(humarray)
library(jimisc)
library(snpStatsWriter)

#first seeing which SNPs weren't dropped in QC:
t1dsnps<-c("rs2476601","rs78037977",  "rs6691977", "rs3024505", "rs13415583",
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

t1dloci<-c("PTPN22", "TNFSF4","CAMSAP2",
"IL10", "AFF3", "ACOXL",
"IFIH1 (1)", "IFIH1 (2)", "IFIH1 (3)",
"CTLA4", "CCR5", "IL2/IL21 (1)", "IL2/IL21 (2)", "CPE",
"IL7R","PTPRK/THEMIS", "BACH2", "CENPW",
"IKZF1", "COBL",
"GLIS3", "IL2RA (1)", "IL2RA (2)", "IL2RA (3)","IL2RA (4)",
"PTEN","INS (1)","INS (2)", "CD69", "IKZF4",
"SH2B3", "GPR183", "LINC01550",
"MEG3","RASGRP1", "CTSH", "IL27", "DEXI (1)",
"DEXI (2)", "CTRB1","IKZF3",
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

hits<-t1dsnps[!t1dsnps$id %in% colnames(all),"id"]
h<-all@snps
r<-all
h1<-read.table(file=paste0(d,"ic-sanger-b58c.bim"),header=F, as.is=T)
colnames(h1)<-c("chromosome", "snp.name", "cM", "position", "allele.1", "allele.2")

#Readin the PLINK file of all the individuals and define a 0.5MB region around the lead SNP for imputation:
imputethem<-function(snp){
h2<-h1[h1$snp.name==snp,]
chr<-h2$chromosome
pos<-h2$position
min<-pos-250000
max<-pos+250000
reg<-h[h$chromosome==chr & h$position>min & h$position<max,]
#allign with 1000G for imputation:
tg<-read.table(file=paste0("/well/1000G/WTCHG/1000GP_Phase3/1000GP_Phase3_chr",chr,".legend.gz"), header=T,as.is=T)
reg<-liftthem(reg, chain="hg18ToHg19.over.chain",
updateto="37")
reg$altstrand1<-ifelse(reg$allele.1=="A","T",
ifelse(reg$allele.1=="T","A",
ifelse(reg$allele.1=="G","C",
ifelse(reg$allele.1=="C","G",NA))))
reg$altstrand2<-ifelse(reg$allele.2=="A","T",
ifelse(reg$allele.2=="T","A",
ifelse(reg$allele.2=="G","C",
ifelse(reg$allele.2=="C","G",NA))))
t<-tg[tg$position %in% reg$position37,]
p<-r[,colnames(r) %in% reg$snp.name]
cs<-col.summary(p)
reg<-reg[reg$snp.name %in% rownames(cs),]
rownames(reg)<-reg$snp.name
reg<-reg[rownames(cs),]
reg<-cbind(reg,cs)
#remove SNPs that are duplicates in the 1000 genomes if they have a MAF<0.01 in europeans:
d<-t[duplicated(t$position),]
d<-t[t$position %in% d$position,]
drop<-d[d$ALL<0.001,]
t<-t[!t$id %in% drop$id,]
#and whatever is left, remove from both our data and reference:
d<-t[duplicated(t$position),]
d<-t[t$position %in% d$position,]
reg<-reg[!reg$position37 %in% d$position,]
t<-t[!t$position %in% d$position,]
reg$position<-reg$position37
b<-merge(reg,t,by="position")
p<-p[,b$snp.name]
#correct the strand allignment problems:
switch=which(b$allele.1==b$a1 & b$allele.2==b$a0)
p<-switch.alleles(p,snps=switch)
cs<-col.summary(p)
reg<-p@snps
reg<-liftthem(reg, chain="hg18ToHg19.over.chain",
updateto="37")
reg$position<-reg$position37
reg<-cbind(reg,cs)
b<-merge(reg,t,by="position")
b$diff<-abs(b$RAF-b$EUR)
message(paste0("removing ",nrow(b[b$diff>0.05,]),"/",nrow(b)," SNPs due to strand ambiguity"))
b<-b[b$diff<0.05,]
rownames(b)<-b$snp.name
p<-p[,p@snps$snp.name %in% b$snp.name]
b<-b[colnames(p),]
colnames(p)<-b$id
write.impute(pedfile=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/",snp,"_n"),
as(p,"SnpMatrix"),
a1=b$allele.1,
a2=b$allele.2,
bp=b$position)
min<-min(b$position)
max<-max(b$position)

#now write a script to run this through impute2:
sink(file=paste0("~/programs/aad/under_7/imputation/",snp,"_n"))
cat(paste0("#!/bin/bash
#$ -cwd -V
#$ -N ",snp," -j y
#$ -P todd.prjc -q long.qc

/apps/well/impute2/2.3.0/impute2 -g /well/todd/users/jinshaw/aad/under_7/imputation/",snp,
"_n -m /well/1000G/WTCHG/1000GP_Phase3/genetic_map_chr",chr,"_combined_b37.txt -int ",min-10000," ",max+10000,
" -h /well/1000G/WTCHG/1000GP_Phase3/1000GP_Phase3_chr",chr,
".hap.gz -l /well/1000G/WTCHG/1000GP_Phase3/1000GP_Phase3_chr",chr,".legend.gz -o /well/todd/users/jinshaw/aad/under_7/imputation/",
snp,"_n_out\n"))
sink()

system(paste0("chmod a=rwx ~/programs/aad/under_7/imputation/",snp,"_n"))
system(paste0("qsub ~/programs/aad/under_7/imputation/",snp,"_n"))
return(p)
}

lapply(hits, imputethem)
