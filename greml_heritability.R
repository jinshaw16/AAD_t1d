#greml_heritability.R
#estimating heritability of <7s, 7-13s and >13s:
library(snpStats)
library(annotSnpStats)
library(snpStatsWriter)
library(humarray)
library(gridExtra)
library(multinomRob)
library(ggplot2)

d<-"/well/todd/users/jinshaw/t1d_risk/immunochip/"
#read SNP and phenotype data:
load(file=paste0(d,"all_inds_unrel_postqc.RData"))

all@samples$onset<-as.numeric(all@samples$onset)
all@samples$group<-ifelse(all@samples$affected==1,0,
ifelse(all@samples$onset<7 & !is.na(all@samples$onset),1,
ifelse(all@samples$onset>=7 & all@samples$onset<13 & !is.na(all@samples$onset),2,
ifelse(all@samples$onset>=13 & !is.na(all@samples$onset),3,NA))))
rownames(all)<-paste0(all@samples$pedigree,".",all@samples$member)

u<-all[rownames(all) %in% rownames(all@samples[all@samples$group %in% c(0,1),]),]
m<-all[rownames(all) %in% rownames(all@samples[all@samples$group %in% c(0,2),]),]
o<-all[rownames(all) %in% rownames(all@samples[all@samples$group %in% c(0,3),]),]
writeit<-function(df,name){
samples<-df@samples
snps<-df@snps
write.plink(file.base=paste0("/well/todd/users/jinshaw/aad/under_7/",name,"_all"),
snps=as(df,"SnpMatrix"),
pedigree=samples$pedigree,
id=samples$member,
father=samples$father,
mother=samples$mother,
sex=samples$sex,
phenotype=samples$affected,
chromosome=snps$chromosome,
genetic.distance=snps$cM,
position=snps$position,
allele.1=snps$allele.1,
allele.2=snps$allele.2)
samples$t1d<-ifelse(samples$affected==2,1,ifelse(samples$affected==1,0,NA))
write.table(samples[,c("pedigree","member","t1d")],file=paste0("/well/todd/users/jinshaw/aad/under_7/pheno_",name),
col.names=F, row.names=F, sep="\t",quote=F)
write.table(samples[,c("pedigree","member","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")],
file=paste0("/well/todd/users/jinshaw/aad/under_7/covars_",name),col.names=F, row.names=F, sep="\t", quote=F)
}
writeit(u,"under_7")
writeit(m,"mid_range")
writeit(o,"over_13")


dogreml<-function(name){
sink(file=paste0("/users/todd/jinshaw/programs/aad/under_7/greml/",name,".sh"))
cat(paste0("/apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 --bfile /well/todd/users/jinshaw/aad/under_7/",name,
"_all --autosome --maf 0.01 --make-grm --out /well/todd/users/jinshaw/aad/under_7/grm_",name," --thread-num 10\n"))
cat(paste0("/apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 --grm /well/todd/users/jinshaw/aad/under_7/grm_",name,
" --pheno /well/todd/users/jinshaw/aad/under_7/pheno_",name,
" --reml --qcovar /well/todd/users/jinshaw/aad/under_7/covars_",name,
" --prevalence 0.004 --out /well/todd/users/jinshaw/aad/under_7/outreml_",name," --thread-num 10\n"))
sink()
}
dogreml("under_7")
dogreml("mid_range")
dogreml("over_13")

dogremlprev<-function(name,prev){
sink(file=paste0("/users/todd/jinshaw/programs/aad/under_7/greml/",name,"_",prev,".sh"))
cat(paste0("/apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 --grm /well/todd/users/jinshaw/aad/under_7/grm_",name,
" --pheno /well/todd/users/jinshaw/aad/under_7/pheno_",name,
" --reml --qcovar /well/todd/users/jinshaw/aad/under_7/covars_",name,
" --prevalence ",prev," --out /well/todd/users/jinshaw/aad/under_7/outreml_",name,"_",prev," --thread-num 10\n"))
sink()
}

dogremlprev("under_7",0.005)
dogremlprev("mid_range", 0.002)
dogremlprev("mid_range", 0.003)
dogremlprev("over_13",0.002)
dogremlprev("over_13",0.003)

#and now excluding the MHC:
s<-all@snps
s<-s[!(s$chromosome==6 & s$position>25000000 & s$position<35000000),]
all<-all[,colnames(all) %in% rownames(s)]
u<-all[rownames(all) %in% rownames(all@samples[all@samples$group %in% c(0,1),]),]
m<-all[rownames(all) %in% rownames(all@samples[all@samples$group %in% c(0,2),]),]
o<-all[rownames(all) %in% rownames(all@samples[all@samples$group %in% c(0,3),]),]


writeit(u,"under_7nomhc")
writeit(m,"mid_rangenomhc")
writeit(o,"over_13nomhc")
dogreml("under_7nomhc")
dogreml("mid_rangenomhc")
dogreml("over_13nomhc")


dogremlprev("under_7nomhc",0.005)
dogremlprev("mid_rangenomhc", 0.002)
dogremlprev("mid_rangenomhc", 0.003)
dogremlprev("over_13nomhc",0.002)
dogremlprev("over_13nomhc",0.003)

