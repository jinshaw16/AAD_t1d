#multinomial_split_iter.R
#re doing the multinomial but splitting the cohort into 2 groups to see if results are consistent between groups
#doing this 1000 times (in the cluster) and counting the proportion of times nominally significant for each region.

args = commandArgs(trailingOnly=TRUE)

library(snpStats)
library(annotSnpStats)
library(snpStatsWriter)
library(humarray)
library(gridExtra)
library(multinomRob)
library(ggplot2)
library(dplyr)
library(withr)

#load data:
load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_n.R")
#sample 50% of the collection (cohort and disease status invariant for now):

pheno$group<-ifelse(pheno$group==2,3,pheno$group)
pheno$group<-ifelse(pheno$onset>=7 & pheno$onset<13 & !is.na(pheno$onset),2,pheno$group)
pheno<-pheno[!is.na(pheno$group),]
sampthem<-function(countr,gp){
p1<- pheno[pheno$country==countr & pheno$group==gp,]
with_seed(args, p1<-sample_n(p1, size=ceiling(nrow(p1)/2)))
return(p1)
}
s<-mapply(sampthem,countr=c("AP","USA","UK","EUR","Finland","NI"),
gp=c(0,0,0,0,0,0),SIMPLIFY=F)
s<-do.call("rbind",s)
su<-mapply(sampthem,countr=c("AP","USA","UK","EUR","Finland","NI"),
gp=c(1,1,1,1,1,1),SIMPLIFY=F)
su<-do.call("rbind",su)
sm<-mapply(sampthem,countr=c("AP","USA","UK","EUR","Finland","NI"),
gp=c(2,2,2,2,2,2),SIMPLIFY=F)
sm<-do.call("rbind",sm)
so<-mapply(sampthem,countr=c("AP","USA","UK","EUR","Finland","NI"),
gp=c(3,3,3,3,3,3),SIMPLIFY=F)
so<-do.call("rbind",so)

s<-rbind(s,su)
s<-rbind(s,sm)
s<-rbind(s,so)
 
s1<-pheno[!pheno$uniqueID %in% s$uniqueID,]
save(s,s1, file=paste0("/well/todd/users/jinshaw/aad/under_7/split_in_half_",args,"_n.RData"))
load(file=paste0("/well/todd/users/jinshaw/aad/under_7/split_in_half_",args,"_n.RData"))

#now perform multinomial regression (remove 7-13s):
domult<-function(pheno,half){
p<-pheno

pheno$g0<-ifelse(pheno$group==0,1,ifelse(pheno$group!=0,0,NA))
pheno$g1<-ifelse(pheno$group==1,1,ifelse(pheno$group!=1,0,NA))
pheno$g2<-ifelse(pheno$group==2,1,ifelse(pheno$group!=2,0,NA))
pheno$g3<-ifelse(pheno$group==3,1,ifelse(pheno$group!=3,0,NA))
t1dsnps$altid<-ifelse(substr(t1dsnps$id,1,1)=="1",paste0("X",t1dsnps$id),t1dsnps$id)
colnames(pheno)<-ifelse(substr(colnames(pheno),1,1)=="1",paste0("X",colnames(pheno)),colnames(pheno))
t1dsnps$altid<-ifelse(!t1dsnps$altid %in% colnames(pheno),t1dsnps$snp,t1dsnps$altid)


#carry out the multinomial regressions:
getlikelihoods<-function(snpname){
#model 1: allow the beta for the SNP to vary for both:
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form3<-as.formula(paste0("g3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))


if(snpname %in% c("imm_2_162845188","seq-rs35667974")){
pheno[,snpname]<-ifelse(pheno[,snpname]==2,0,ifelse(pheno[,snpname]==0,2,pheno[,snpname]))
}
one <- tryCatch(multinomRob(model=list(form0,form1,form2, form3),data=pheno, print.level=1, MLEonly=T), error=function(err) NA)
if(is.na(one)){
l<-data.frame(unconstrained=NA, constrained=NA, logor1=NA, logse1=NA,
logor2=NA, logse2=NA)
}
if(!is.na(one)){
llk1<-one$value

equals0<-as.formula(paste0("g1 ~ `",snpname,"` + 0"))
equals1<-as.formula(paste0("g3 ~ `",snpname,"` + 0"))
#model 2: constrain betas to equal each other:
two<-multinomRob(model=list(form0,form1,form2, form3),
equality=list(list(equals0,equals1)), data=pheno, print.level=1, MLEonly=T)
llk2<-two$value
l<-data.frame(unconstrained=llk1, constrained=llk2, logor1=one$coefficients[12,2], logse1=one$se[12,2],
logor2=one$coefficients[12,3], logse2=one$se[12,3])
}
return(l)
}

likelihoods<-lapply(t1dsnps$altid,getlikelihoods)
likelihoods<-do.call("rbind", likelihoods)
likelihoods$snp=t1dsnps$snp
likelihoods$id<-t1dsnps$id
likelihoods$loci<-t1dsnps$loci
likelihoods$altid<-t1dsnps$altid
likelihoods$constrained=likelihoods$constrained*-1
likelihoods$unconstrained=likelihoods$unconstrained*-1
likelihoods$loglambda<-likelihoods$constrained-likelihoods$unconstrained
likelihoods$chisq<-likelihoods$loglambda*-2
likelihoods$p<-pchisq(likelihoods$chisq,1, lower.tail=F)
rownames(likelihoods)<-likelihoods$altid

likelihoods$lb1<-likelihoods$logor1-(qnorm(0.975)*likelihoods$logse1)
likelihoods$ub1<-likelihoods$logor1+(qnorm(0.975)*likelihoods$logse1)
likelihoods$lb2<-likelihoods$logor2-(qnorm(0.975)*likelihoods$logse2)
likelihoods$ub2<-likelihoods$logor2+(qnorm(0.975)*likelihoods$logse2)
likelihoods$logp<-log10(likelihoods$p)*-1

write.table(likelihoods, file=paste0("/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_iter_",args,"_",half,"_n.txt"),
col.names=T, row.names=F, quote=F, sep="\t")
}

mapply(domult, pheno=list(s,s1), half=c(1,2))
