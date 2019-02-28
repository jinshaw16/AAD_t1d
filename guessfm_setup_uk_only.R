#guessfm_setup_uk_only.R

library(GUESSFM)
library(R2GUESS)
library(reshape2)
library(vcfR)
library(ggplot2)
library(ggbio)
library(snpStats)
library(jimisc)
library(annotSnpStats)

#want to run GUESSFM on each region.
#define regions (those passing Bonferroni correction in primary analysis):
hits<-c("imm_1_170941164", "imm_4_123335627", "imm_6_128328079", "imm_15_77022012",
"imm_17_35158633", "imm_9_4281928",   "imm_1_205006527","imm_16_73809828",
"imm_10_6170083")


#load imp
load(file="/well/todd/users/jinshaw/aad/under_7/imputation/loadeddata.RData")

#read in sample info:
d<-"/well/todd/users/jinshaw/t1d_risk/immunochip/"
#read SNP and phenotype data:
load(file=paste0(d,"all_inds_unrel_postqc.RData"))
pheno<-pall@samples
pheno$onset<-as.numeric(pheno$onset)
pheno$group<-ifelse(pheno$affected==1,0,
ifelse(pheno$onset<7 & !is.na(pheno$onset),1,
ifelse(pheno$onset>=7 & pheno$onset<13 & !is.na(pheno$onset),2,
ifelse(pheno$onset>=13 & !is.na(pheno$onset),3,NA))))

#keep only those diagnosed at <7 and from the UK:
sample<-pheno[(pheno$group==0 | pheno$group==1) & !is.na(pheno$group) & pheno$country %in% c("NI","UK"),]

sample$t1d<-ifelse(sample$affected==2,1,ifelse(sample$affected==1,0,NA))
tcols<-rownames(ps[[1]])

#as part of QC, want to run the association statsictics through SNPTEST so can remove artififactually associated SNPs prior to running GEUSSFM:
sample$missing<-0
s<-sample[,c("uniqueID")]
write.table(s, file="/well/todd/users/jinshaw/aad/under_7/imputation/uk",col.names=F,row.names=F,quote=F)


dosnptest<-function(snp){
sink(file=paste0("~/programs/aad/under_7/snptest_scripts/",snp,"_uk.sh"))
cat(paste0("qctool -g /well/todd/users/jinshaw/aad/under_7/imputation/",snp,
"_out -s /well/todd/users/jinshaw/aad/under_7/imputation/geno.sample -filetype gen -incl-samples ",
"/well/todd/users/jinshaw/aad/under_7/imputation/uk -os /well/todd/users/jinshaw/aad/under_7/imputation/uk.sample_",snp,
" -og /well/todd/users/jinshaw/aad/under_7/imputation/uk_out_",snp," -ofiletype gen\n"))

cat(paste0("/apps/well/snptest/2.5.4-beta3_CentOS6.6_x86_64_dynamic/snptest_v2.5.4-beta3 ",
"-data /well/todd/users/jinshaw/aad/under_7/imputation/uk_out_",snp," /well/todd/users/jinshaw/aad/under_7/imputation/uk.sample_",snp,
" -pheno t1d -frequentist 1 -method newml -cov_all_continuous -o /well/todd/users/jinshaw/aad/under_7/imputation/snptest/",snp,"_uk"))
sink()
create_job(path=paste0("~/programs/aad/under_7/snptest_scripts/"),subname=paste0(snp,"_uk"),
jobname=paste0(snp,"_uk"),projectletter="c", qletter="c", qlength="short")
system(paste0("chmod a=rwx ~/programs/aad/under_7/snptest_scripts/",snp,"_uk.sh"))
system(paste0("qsub ~/programs/aad/under_7/snptest_scripts/",snp,"_uk"))
}
lapply(hits, dosnptest)

#run GUESSFM
runguessfm<-function(snp){
system(paste0("cut -d \" \" -f 2- /well/todd/users/jinshaw/aad/under_7/imputation/uk_out_",snp,
" > /well/todd/users/jinshaw/aad/under_7/imputation/uk_out_",snp,"_1"))

tcols<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/uk.sample_",snp), skip=2, as.is=T)
tcols<-tcols$V1

DATA<-read.impute(file=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/uk_out_",snp,"_1"), rownames=tcols)
DATA<-DATA[rownames(sample),]
info<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/",snp,"_out_info"),header=T)
Y<-data.frame(outcome=sample$t1d)
covariates<-sample[,c("PC1","PC2","PC3","PC4","PC5")]

rownames(Y)<-rownames(sample)
rownames(covariates)<-rownames(sample)


#filter poorly imputed SNPs:
message("Removing poorly imputed SNPs.")
message("Input matrix has ",ncol(DATA)," SNPs.")
cs <- col.summary(DATA)
wh <- which(cs[,"MAF"]<0.005 | cs[,"Call.rate"]<0.99 | cs[,"Certain.calls"]<0.75 | abs(cs[,"z.HWE"])>20 | is.na(cs[,"z.HWE"]) | info$info<0.8)

if(length(wh)) {
  message("Dropping ",length(wh)," SNPs with |z.HWE|>20, certain calls <0.75, MAF < 0.005, call rate <0.99 or info score<0.8")
  DATA <- DATA[,-wh]
}
DATA<-DATA[,!duplicated(colnames(DATA))]
cs$rs_id<-rownames(cs)
cs<-cs[!duplicated(cs$rs_id),]
info<-info[!duplicated(info$rs_id),]
cs<-merge(cs,info,by="rs_id")
rownames(cs)<-cs$rs_id
cs<-cs[colnames(DATA),]
cs$diff<-abs(cs$RAF-cs$exp_freq_a1)
cs<-cs[cs$diff<0.05,]

cs<-cs[cs$rs_id %in% colnames(DATA),]
rownames(cs)<-cs$rs_id
cs<-cs[colnames(DATA),]


#get info for cases and controls:
res<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/snptest/",snp,"_uk"),header=T,as.is=T)
res<-res[res$rsid %in% rownames(cs),]
res<-res[!duplicated(res$rsid),]
rownames(res)<-res$rsid
res<-res[rownames(cs),]
res$diff<-abs(res$cases_info-res$controls_info)
w<-which(res$diff>0.01)

message(paste0("remove an additional ",length(w), "SNPs due to imputation r2 differences between cases and controls of >1%"))
cs<-cs[-w,]
DATA<-DATA[rownames(sample),]
DATA<-DATA[,-w]

#Run bayesian variable selection via GUESS
mydir <-"/well/todd/users/jinshaw/aad/under_7/guessfm/uk/"
save(Y, DATA, covariates, file=paste0(mydir,snp, "/data.RData"))
load(file=paste0(mydir, snp,"/data.RData"))


com<-run.bvs(X=DATA,Y=Y[,"outcome"],gdir=paste0(mydir,snp),
        guess.command="/users/todd/jinshaw/software/GUESS_v1.1/Main/GUESS",
        nexp=3,                # expected number of causal variants, an overestimate
        nsave=10000,            # number of best models to save
        tag.r2=0.99, 	 	#R2 tag value
        family="binomial",
	covars=covariates, run=FALSE)      
sink(file=paste0("~/programs/aad/under_7/guessscripts/uk_",snp,".sh"))
cat(paste0("#!/bin/bash
#$ -cwd -V
#$ -N ",snp," -j y
#$ -P todd.prjc -q long.qc

#uk_",snp,".sh
",com,"\n"))
sink()
system(paste0("chmod a=rwx ~/programs/aad/under_7/guessscripts/uk_",snp,".sh"))
system(paste0("qsub ~/programs/aad/under_7/guessscripts/uk_",snp,".sh"))
}
commands<-lapply(hits,runguessfm)




