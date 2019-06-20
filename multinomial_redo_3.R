#multinomial_redo_3.R
library(snpStats)
library(annotSnpStats)
library(snpStatsWriter)
library(humarray)
library(gridExtra)
library(multinomRob)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

d<-"/well/todd/users/jinshaw/t1d_risk/immunochip/"
#read SNP and phenotype data:
#load(file=paste0(d,"all_inds_unrel_postqc_3.RData"))

#all@samples$onset<-as.numeric(all@samples$onset)
#all@samples$group<-ifelse(all@samples$affected==1,0,
#ifelse(all@samples$onset<7 & !is.na(all@samples$onset),1,
#ifelse(all@samples$onset>=7 & all@samples$onset<13 & !is.na(all@samples$onset),2,
#ifelse(all@samples$onset>=13 & !is.na(all@samples$onset),3,NA))))


#Define the 55 T1D regions we want to examine (from various literature):
#t1dsnps<-c("rs2476601","rs78037977",  "rs6691977", "rs3024505", "rs13415583",
#           "rs4849135", "rs2111485", "rs35667974", "rs72871627",
#"rs3087243", "rs113010081", "rs6819058","rs67797421", "rs2611215", "rs11954020",
#"rs72975913","rs72928038", "rs1538171", "rs62447205", "rs10277986", "rs6476839",
#"rs61839660", "rs11594656","rs6602437", "rs41295121",
#"rs12416116", "rs689", "rs72853903", "rs917911", "rs705705",
#"rs653178", "rs9585056", "rs1456988", "rs56994090", "rs72727394",
#"rs34593439", "rs151234", "rs12927355",
#"rs193778", "rs8056814", "rs12453507",
#"rs757411", "rs1052553", "rs1893217", "rs12971201",
#"rs1615504", "rs34536443", "rs12720356", "rs402072", "rs516246",
#"rs6043409", "rs11203202", "rs6518350", "rs4820830", "rs229533")

#t1dloci<-c("PTPN22", "TNFSF4","CAMSAP2",
#"IL10", "AFF3", "ACOXL",
#"IFIH1 (1)", "IFIH1 (2)", "IFIH1 (3)",
#"CTLA4", "CCR5", "IL2/IL21 (1)", "IL2/IL21 (2)", "CPE",
#"IL7R","THEMIS", "BACH2", "CENPW",
#"IKZF1", "COBL",
#"GLIS3", "IL2RA (1)", "IL2RA (2)", "IL2RA (3)","IL2RA (4)",
#"PTEN","INS (1)","INS (2)", "CD69", "IKZF4",
#"SH2B3", "GPR183", "LINC01550",
#"MEG3","RASGRP1", "CTSH", "IL27", "DEXI (1)",
#"DEXI (2)", "CTRB1","IKZF3",
#"CCR7", "MAPT","PTPN2 (1)", "PTPN2 (2)", "CD226", "TYK2 (1)", "TYK2 (2)", "PRKD2",
#"FUT2", "SIRPG",
#"UBASH3A", "ICOSLG", "HORMAD2",
#"C1QTNF6")


#t1dsnps<-data.frame(snp=t1dsnps, loci=t1dloci)
#t1dsnps$id<-rs.to.id(t1dsnps$snp)
#t1dsnps$id<-gsub(".*,","",t1dsnps$id)
#t1dsnps$id<-ifelse(substr(t1dsnps$id,1,3)=="seq",gsub("_","-", t1dsnps$id),t1dsnps$id)
#t1dsnps$id<-ifelse(substr(t1dsnps$id,1,4)=="X1kg", substr(t1dsnps$id,2,100000),t1dsnps$id)
#t1dsnps$ord<-c(1:nrow(t1dsnps))
#t1dsnps$snp<-as.character(t1dsnps$snp)
#p<-t1dsnps[!t1dsnps$id %in% colnames(all),]


#keep only the ones we want to test:
#g<-all[,colnames(all) %in% t1dsnps$id]
#g<-as(g,"SnpMatrix")
#cs<-col.summary(g)
#w<-which(cs$RAF>0.5)
#g<-switch.alleles(g,snps=w)
#b<-as(g,"numeric")

#and the imputed ones:
getsnp<-function(snp){
load(file=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/",snp,"_3_clean.RData"))
cs<-col.summary(DATA)
w<-which(cs$RAF>0.5)
DATA<-switch.alleles(DATA,snps=w)
d<-as(DATA[rownames(DATA) %in% rownames(g),grepl(snp,colnames(DATA))],"numeric")
colnames(d)<-snp
b<<-cbind2(b,d)
return(d)
}
#l<-lapply(p$id,getsnp)
#colnames(b)<-gsub(":.*","",colnames(b))
#pheno<-all@samples

#table(rownames(b)==rownames(pheno))
#pheno<-cbind(pheno,b)
#pheno<-pheno[!is.na(pheno$sex) & pheno$sex!=0,]
#p1<-pheno
#save(b,pheno,t1dsnps, file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_3.R")

load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_3.R")
t1dsnps$altid<-ifelse(substr(t1dsnps$id,1,1)=="1",paste0("X",t1dsnps$id),t1dsnps$id)
colnames(pheno)<-ifelse(substr(colnames(pheno),1,1)=="1",paste0("X",colnames(pheno)),colnames(pheno))
t1dsnps$altid<-ifelse(!t1dsnps$altid %in% colnames(pheno),t1dsnps$snp,t1dsnps$altid)

#now perform multinomial regression (remove 7-13s for the first 2 models):
pheno$group<-ifelse(pheno$affected==1,0,
ifelse(pheno$onset<7 & !is.na(pheno$onset),1,
ifelse(pheno$onset>=7 &pheno$onset<13 & !is.na(pheno$onset),2,
ifelse(pheno$onset>=13 & !is.na(pheno$onset),3,NA))))
pheno$g0<-ifelse(pheno$group==0,1,ifelse(pheno$group!=0,0,NA))
pheno$g1<-ifelse(pheno$group==1,1,ifelse(pheno$group!=1,0,NA))
pheno$g2<-ifelse(pheno$group==2,1,ifelse(pheno$group!=2,0,NA))
pheno$g3<-ifelse(pheno$group==3,1,ifelse(pheno$group!=3,0,NA))



#carry out the multinomial regressions:
getlikelihoods<-function(snpname){
#model 1: allow the beta for the SNP to vary for both:
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + sex + `",snpname,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + sex + `",snpname,"`"))
form3<-as.formula(paste0("g3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + sex + `",snpname,"`"))

one<-multinomRob(model=list(form0,form1,form2,form3),
data=pheno, print.level=1, MLEonly=T)
llk1<-one$value

equals1<-as.formula(paste0("g1 ~ `",snpname,"` + 0"))
equals2<-as.formula(paste0("g3 ~ `",snpname,"` + 0"))
#model 2: constrain betas to equal each other:
two<-multinomRob(model=list(form0,form1,form2,form3),
equality=list(list(equals1,equals2)), data=pheno, print.level=1, MLEonly=T)
llk2<-two$value

l<-data.frame(unconstrained=llk1, constrained=llk2,
logor1=one$coefficients[nrow(one$coefficients),2], logse1=one$se[nrow(one$coefficients),2],
logor2=one$coefficients[nrow(one$coefficients),3], logse2=one$se[nrow(one$coefficients),3],
logor3=one$coefficients[nrow(one$coefficients),4],logse3=one$se[nrow(one$coefficients),4])
if (snpname=="seq-rs35667974"){
l$logor1<-l$logor1*-1
l$logor2<-l$logor2*-1
l$logor3<-l$logor3*-1
}
l$lb1<-l$logor1-(qnorm(0.975)*l$logse1)
l$ub1<-l$logor1+(qnorm(0.975)*l$logse1)
l$lb2<-l$logor2-(qnorm(0.975)*l$logse2)
l$ub2<-l$logor2+(qnorm(0.975)*l$logse2)
l$lb3<-l$logor3-(qnorm(0.975)*l$logse3)
l$ub3<-l$logor3+(qnorm(0.975)*l$logse3)
l$constrained=l$constrained*-1
l$unconstrained=l$unconstrained*-1
l$loglambda<-l$constrained-l$unconstrained
#compare using likelihood ratio test:
l$chisq<-l$loglambda*-2
l$p<-pchisq(l$chisq,1, lower.tail=F)
save(l,file=paste0("/well/todd/users/jinshaw/aad/under_7/results/",snpname,"_3.RData"))
}
pheno$`seq-rs35667974`<-ifelse(pheno$`seq-rs35667974`==2,0,ifelse(pheno$`seq-rs35667974`==0,2,pheno$`seq-rs35667974`))
getlikelihoods(t1dsnps$altid[as.numeric(args)])


