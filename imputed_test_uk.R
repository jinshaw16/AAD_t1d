#imputed_test_uk.R
#univariable tests to get log-odds ratios for the under 7s at each of the credible SNPs:

library(ggplot2)
library(ggbio)
library(GUESSFM)
library(R2GUESS)
library(reshape2)
library(ggplot2)
library(snpStats)
library(speedglm)

args = commandArgs(trailingOnly=TRUE)


#set paths
mydir <-"/well/todd/users/jinshaw/aad/under_7/guessfm/uk/"
outdir<-"/well/todd/users/jinshaw/output/aad/under_7/guessfm/uk/"

load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_uk_2.R")

#keep only those diagnosed at <7:
pheno<-pheno[(pheno$group==0 | pheno$group==1) & !is.na(pheno$group),]
pheno$t1d<-ifelse(pheno$affected==2,1,ifelse(pheno$affected==1,0,NA))

getodds<-function(snp){
load(file=paste0(mydir,snp, "/data.RData"))

mydir<-paste0(mydir,snp,"/")
colnames(DATA)<-gsub(":",".",colnames(DATA))
colnames(DATA)<-gsub("<",".",colnames(DATA))
colnames(DATA)<-gsub(">",".",colnames(DATA))

colnames(DATA)<-ifelse(substr(colnames(DATA),1,1) %in% c("1","2","3","4","5","6","7","8","9"),
paste0("X",colnames(DATA)),colnames(DATA))

load(file=paste0(mydir,"summx.RData"))

#get the alleles:
system(paste0("awk -F \' \' \' {print $1,$2,$3,$4,$5}\' /well/todd/users/jinshaw/aad/under_7/imputation/",
snp,"_out_uk_n > /well/todd/users/jinshaw/aad/under_7/imputation/",snp,"alleles"))
alls<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/",snp,"alleles"),header=F,as.is=T)
alls$V2<-gsub(":",".",alls$V2)
alls$V2<-gsub("<",".",alls$V2)
alls$V2<-gsub(">",".",alls$V2)
alls$V2<-ifelse(substr(alls$V2,1,1) %in% c("1","2","3","4","5","6","7","8","9"),
paste0("X",alls$V2),alls$V2)

DATA<-DATA[,colnames(DATA) %in% summx$snp]
DATA<-DATA[pheno$uniqueID,]
alls<-alls[alls$V2 %in% colnames(DATA),]
alls<-alls[!duplicated(alls$V2),]
rownames(alls)<-alls$V2
alls<-alls[colnames(DATA),]
s<-snp.rhs.estimates(formula=group ~ PC1 + PC2 + PC3 + PC4 + PC5, data=pheno, link="logit", family="binomial",
snp.data=DATA, uncertain=TRUE)
s<-do.call("rbind",s)
s1<-as.data.frame(s)
for(i in 1:4){
s1[,i]<-unlist(s1[,i])
}
s1<-s1[alls$V2,]
s1$ref=alls$V4
s1$effect<-alls$V5
s1$snp=rownames(s1)
cs<-col.summary(DATA)
if (snp %in% c("imm_17_35306733","imm_9_4280823")){
w<-which(s1$beta>0)
}
if (!snp %in% c("imm_17_35306733","imm_9_4280823")){
w<-which(cs$RAF>0.5)
}
s1[w,"beta"]<-s1[w,"beta"]*-1
s1[w,"ref"]<-alls[w,"V5"]
s1[w,"effect"]<-alls[w,"V4"]
write.table(s1, file=paste0("/well/todd/users/jinshaw/aad/under_7/results/",snp,"_effects_uk.txt"),
quote=F, col.names=T, row.names=F, sep="\t")
}

getodds(args[1])
