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
likelihoods<-read.table(file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_n.txt",
header=T, sep="\t", as.is=T)
likelihoods<-likelihoods[likelihoods$p<(0.05/nrow(likelihoods)),]
hits<-likelihoods$id


#read in sample info (UK only):
#r<-read.plink(fam="/well/todd/users/jinshaw/aad/all_ages_n.fam",
#bed="/well/todd/users/jinshaw/aad/all_ages_n.bed",
#bim="/well/todd/users/jinshaw/aad/all_ages_n.bim")
#r<-annot.plink(r)

load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_uk_n.R")

#keep only those diagnosed at <7:
sample<-pheno[(pheno$group==0 | pheno$group==1) & !is.na(pheno$group),]

sample$t1d<-ifelse(sample$affected==2,1,ifelse(sample$affected==1,0,NA))
#r<-r[rownames(r) %in% rownames(sample),]

loaddata<-function(snp){
h<-r@snps[r@snps$snp.name==snp,]
chr<-h$chromosome
pos<-h$position
min<-pos-250000
max<-pos+250000
reg<-r@snps[r@snps$chromosome==chr & r@snps$position>min & r@snps$position<max,]
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
write.impute(pedfile=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/uk_",snp),
as(p,"SnpMatrix"),
a1=b$allele.1,
a2=b$allele.2,
bp=b$position)
min<-min(b$position)
max<-max(b$position)
#now write a script to run this through impute2:
sink(file=paste0("~/programs/aad/under_7/imputation/uk_",snp))
cat(paste0("#!/bin/bash
#$ -cwd -V
#$ -N ",snp," -j y
#$ -P todd.prjc -q long.qc

/apps/well/impute2/2.3.0/impute2 -g /well/todd/users/jinshaw/aad/under_7/imputation/uk_",snp,
" -m /well/1000G/WTCHG/1000GP_Phase3/genetic_map_chr",chr,"_combined_b37.txt -int ",min-10000," ",max+10000,
" -h /well/1000G/WTCHG/1000GP_Phase3/1000GP_Phase3_chr",chr,
".hap.gz -l /well/1000G/WTCHG/1000GP_Phase3/1000GP_Phase3_chr",chr,".legend.gz -o /well/todd/users/jinshaw/aad/under_7/imputation/",
snp,"_out_uk_n\n"))
sink()
system(paste0("qsub ~/programs/aad/under_7/imputation/uk_",snp))
return(p)
}

#ps<-lapply(hits,loaddata)
#names(ps)<-hits
#save(ps,file="/well/todd/users/jinshaw/aad/under_7/imputation/loadeddata_uk.RData")




#as part of QC, want to run the association statsictics through SNPTEST so can remove artififactually associated SNPs prior to running GEUSSFM:
sample$missing<-0
sample$t1d<-ifelse(sample$affected==2,1,ifelse(sample$affected==1,0,NA))
s<-sample[,c("uniqueID","member","missing","t1d","PC1","PC2","PC3","PC4","PC5")]
colnames(s)<-c("ID_1","ID_2","missing","t1d","PC1","PC2","PC3","PC4","PC5")
for (i in 1:ncol(s)){
s[,i]<-as.character(s[,i])
}
headit<-data.frame(ID_1=0,ID_2=0,missing=0,t1d="B",PC1="C",PC2="C",PC3="C",PC4="C",PC5="C")
s<-rbind(headit,s)
write.table(s, file="/well/todd/users/jinshaw/aad/under_7/imputation/uk.geno",col.names=T,row.names=F,quote=F,sep=" ")


dosnptest<-function(snp){
sink(file=paste0("~/programs/aad/under_7/snptest_scripts/",snp,"_uk.sh"))
cat(paste0("/apps/well/snptest/2.5.4-beta3_CentOS6.6_x86_64_dynamic/snptest_v2.5.4-beta3 ",
"-data /well/todd/users/jinshaw/aad/under_7/imputation/",snp,"_out_uk_n /well/todd/users/jinshaw/aad/under_7/imputation/uk.geno",
" -pheno t1d -frequentist 1 -method newml -cov_all_continuous -o /well/todd/users/jinshaw/aad/under_7/imputation/snptest/",snp,"_uk"))
sink()
create_job(path=paste0("~/programs/aad/under_7/snptest_scripts/"),subname=paste0(snp,"_uk"),
jobname=paste0(snp,"_uk"),projectletter="c", qletter="c", qlength="short")
system(paste0("chmod a=rwx ~/programs/aad/under_7/snptest_scripts/",snp,"_uk.sh"))
system(paste0("qsub ~/programs/aad/under_7/snptest_scripts/",snp,"_uk"))
}
#lapply(hits, dosnptest)

#run GUESSFM
load(file="/well/todd/users/jinshaw/aad/under_7/imputation/loadeddata_uk.RData")
tcols<-rownames(ps[[1]])
runguessfm<-function(snp){
DATA<-read.impute(file=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/",snp,"_out_uk_n"), rownames=tcols)
DATA<-DATA[rownames(sample),]
info<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/",snp,"_out_uk_n_info"),header=T)
Y<-data.frame(outcome=sample$t1d)
covariates<-sample[,c("PC1","PC2","PC3","PC4","PC5")]

rownames(Y)<-rownames(sample)
rownames(covariates)<-rownames(sample)


#filter poorly imputed SNPs:
message("Removing poorly imputed SNPs.")
message("Input matrix has ",ncol(DATA)," SNPs.")
cs <- col.summary(DATA)
wh <- which(cs[,"MAF"]<0.005 | cs[,"Call.rate"]<0.99 | cs[,"Certain.calls"]<0.75 | abs(cs[,"z.HWE"])>5.45 | is.na(cs[,"z.HWE"]) | info$info<0.8)

if(length(wh)) {
  message("Dropping ",length(wh)," SNPs with |z.HWE|>5.45, certain calls <0.75, MAF < 0.005, call rate <0.99 or info score<0.8")
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
wh1<-which(cs$type==2 & cs$r2_type0<0.8)
cs<-cs[-wh1,]
DATA<-DATA[,colnames(DATA) %in% rownames(cs)]
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

#finally comparing variances manually and dropping those where the expected variance is different (higher or lower) than the imputed:
#h<-as(DATA,"numeric")
#va<-apply(h,2,var)
#cs$var<-va
#cs$expvar<-2*cs$MAF*(1-cs$MAF)
#cs$varrat<-cs$var/cs$expvar
#w2<-which(cs$varrat<0.98 | cs$varrat>1.02)
#cs<-cs[-w2,]
#DATA<-DATA[,-w2]

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




