#multinomial_uk_only_sex_adj_3.R
library(snpStats)
library(annotSnpStats)
library(snpStatsWriter)
library(humarray)
library(gridExtra)
library(multinomRob)
library(ggplot2)

args=commandArgs(trailingOnly=TRUE)
d<-"/well/todd/users/jinshaw/t1d_risk/immunochip/"
#read SNP and phenotype data:
load(file=paste0(d,"all_inds_unrel_postqc_3.RData"))

uk<-all@samples
#output snptest file for the use in GUESSFM (which imputes one additional poerson by accident)
uk<-uk[uk$country %in% c("UK","NI"),]
#s1<-uk[,c("uniqueID","affected","PC1","PC2","PC3","PC4","PC5","sex")]
#headit<-data.frame(ID_1=0,ID_2=0,missing=0,t1d="B",PC1="C",PC2="C",PC3="C",PC4="C",PC5="C",sex="D")
#s1$ID_1=s1$uniqueID
#s1$ID_2=s1$uniqueID
#s1$missing<-0
#s1$t1d<-ifelse(s1$affected==2,1,ifelse(s1$affected==1,0,NA))
#s1<-s1[,c("ID_1","ID_2","missing","t1d","PC1","PC2","PC3","PC4","PC5","sex")]
#for (vars in colnames(s1)){
#s1[,vars]<-as.character(s1[,vars])
#}
#s1<-rbind(headit,s1)
#write.table(s1,file="/well/todd/users/jinshaw/aad/under_7/imputation/geno_sexadj.sample",col.names=T,row.names=F,sep=" ",quote=F)

#uk<-uk[uk$country %in% c("UK","NI") & !is.na(uk$country) & !is.na(uk$sex) & uk$sex!=0,]
#all<-all[rownames(all) %in% rownames(uk),]
#all<-all[!is.na(all@samples$onset) | all@samples$affected==1,]
#samples<-all@samples
#snps<-all@snps
#write.plink(file.base="/well/todd/users/jinshaw/aad/geno_unrel_uk_3",
#snps=as(all,"SnpMatrix"),
#pedigree=samples$pedigree,
#id=rownames(samples),
#father=samples$father,
#mother=samples$mother,
#sex=samples$sex,
#phenotype=samples$affected,
#chromosome=snps$chromosome,
#position=snps$position,
#allele.1=snps$allele.1,
#allele.2=snps$allele.2)

#get pcs
#system(paste0("plink --bfile /well/todd/users/jinshaw/aad/geno_unrel_uk_3 --indep-pairwise 1000 50 0.2 --out /well/todd/users/jinshaw/aad/under_7/pruned_uk --allow-no-sex"))
#system(paste0("plink --bfile /well/todd/users/jinshaw/aad/geno_unrel_uk_3 --exclude /well/todd/users/jinshaw/aad/under_7/pruned_uk.prune.out --keep-allele-order",
#" --make-bed --out /well/todd/users/jinshaw/aad/under_7/forpcad_uk_3"))
#system(paste0("~/software/plink2 --bfile /well/todd/users/jinshaw/aad/under_7/forpcad_uk_3 ",
#"--exclude range /well/todd/users/jinshaw/aad/under_7/mhc.txt --make-bed --out /well/todd/users/jinshaw/aad/under_7/forpcad_uk_nomhc_3"))
#system(paste0("~/software/plink2 --bfile /well/todd/users/jinshaw/aad/under_7/forpcad_uk_nomhc_3 --pca approx 10 ",
#"--allow-no-sex --out /well/todd/users/jinshaw/aad/under_7/pcad_uk_3"))


#read this genotype data in and the pcas (only SNPs we're interested in), then keep only the hits we're interested in:

#load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_3.R")
#pheno<-pheno[pheno$country %in% c("UK","NI") & !is.na(pheno$country),]
#rownames(pheno)<-pheno$uniqueID
#pheno<-pheno[rownames(uk),]

#pcs<-read.table(file="/well/todd/users/jinshaw/aad/under_7/pcad_uk_3.eigenvec",as.is=T, header=F)
#pcs<-pcs[,c(1:12)]
#colnames(pcs)<-c("ped","mem",paste0("PC",1:10))
#rownames(pcs)<-pcs$mem
#pcs<-pcs[rownames(pheno),]

#replace the PCs from whole cohort with UK specific ones:
#for (i in 1:10){
#pheno[,paste0("PC",i)]<-pcs[,paste0("PC",i)]
#}

#pheno$group<-ifelse(pheno$affected==1,0,
#ifelse(pheno$onset<7 & !is.na(pheno$onset),1,
#ifelse(pheno$onset>=7 & pheno$onset<13 & !is.na(pheno$onset),2,
#ifelse(pheno$onset>=13 & !is.na(pheno$onset),3,NA))))

#save(pheno,t1dsnps,file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_uk_3.R")


load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_uk_3.R")
t1dsnps$altid<-ifelse(substr(t1dsnps$id,1,1)=="1",paste0("X",t1dsnps$id),t1dsnps$id)
colnames(pheno)<-ifelse(substr(colnames(pheno),1,1)=="1",paste0("X",colnames(pheno)),colnames(pheno))
t1dsnps$altid<-ifelse(!t1dsnps$altid %in% colnames(pheno),t1dsnps$snp,t1dsnps$altid)
pheno<-pheno[!is.na(pheno$sex),]

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
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + sex + `",snpname,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + sex + `",snpname,"`"))
form3<-as.formula(paste0("g3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + sex + `",snpname,"`"))

one<-multinomRob(model=list(form0,form1,form2, form3),data=pheno, print.level=1, MLEonly=T)
llk1<-one$value

equals0<-as.formula(paste0("g1 ~ `",snpname,"` + 0"))
equals1<-as.formula(paste0("g3 ~ `",snpname,"` + 0"))
#model 2: constrain betas to equal each other:
two<-multinomRob(model=list(form0,form1,form2,form3),
equality=list(list(equals0,equals1)), data=pheno, print.level=1, MLEonly=T)
llk2<-two$value
l<-data.frame(unconstrained=llk1, constrained=llk2,logor1=one$coefficients[nrow(one$coefficients),2], logse1=one$se[nrow(one$coefficients),2],
logor2=one$coefficients[nrow(one$coefficients),3], logse2=one$se[nrow(one$coefficients),3],
logor3=one$coefficients[nrow(one$coefficients),4],logse3=one$se[nrow(one$coefficients),4])
if (snpname=="seq-rs35667974"){
l$logor1<-l$logor1*-1
l$logor2<-l$logor2*-1
l$logor3<-l$logor3*-1
}
l$lb1=l$logor1-(qnorm(0.975)*l$logse1)
l$lb2=l$logor2-(qnorm(0.975)*l$logse2)
l$lb3=l$logor3-(qnorm(0.975)*l$logse3)
l$ub1=l$logor1+(qnorm(0.975)*l$logse1)
l$ub2=l$logor2+(qnorm(0.975)*l$logse2)
l$ub3=l$logor3+(qnorm(0.975)*l$logse3)
l$constrained=l$constrained*-1
l$unconstrained=l$unconstrained*-1
l$loglambda<-l$constrained-l$unconstrained
l$chisq<-l$loglambda*-2
l$p<-pchisq(l$chisq,1, lower.tail=F)
save(l,file=paste0("/well/todd/users/jinshaw/aad/under_7/results/uk/",snpname,"_sex_adj.RData"))
}
pheno$`seq-rs35667974`<-ifelse(pheno$`seq-rs35667974`==2,0,ifelse(pheno$`seq-rs35667974`==0,2,pheno$`seq-rs35667974`))

getlikelihoods(t1dsnps$altid[as.numeric(args)])

