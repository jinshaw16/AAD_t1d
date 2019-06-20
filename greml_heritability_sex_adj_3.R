#greml_heritability_sex_adj_3.R
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
load(file=paste0(d,"all_inds_unrel_postqc_3.RData"))

all<-all[rownames(all) %in% rownames(all@samples[!is.na(all@samples$sex) & all@samples$sex!=0,]),]
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
write.plink(file.base=paste0("/well/todd/users/jinshaw/aad/under_7/",name,"_all_sexadj"),
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
write.table(samples[,c("pedigree","member","t1d")],file=paste0("/well/todd/users/jinshaw/aad/under_7/pheno_",name,"_sexadj"),
col.names=F, row.names=F, sep="\t",quote=F)
write.table(samples[,c("pedigree","member","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","sex")],
file=paste0("/well/todd/users/jinshaw/aad/under_7/covars_",name,"_sexadj"),col.names=F, row.names=F, sep="\t", quote=F)
system(paste0("plink --bfile /well/todd/users/jinshaw/aad/under_7/",name,
"_all_sexadj --indep-pairwise 1000 50 0.2 --out /well/todd/users/jinshaw/aad/under_7/",name,"_all_sexadj_pruned"))
system(paste0("plink --bfile /well/todd/users/jinshaw/aad/under_7/",name,
"_all_sexadj --exclude /well/todd/users/jinshaw/aad/under_7/",name,"_all_sexadj_pruned.prune.out --make-bed --out /well/todd/users/jinshaw/aad/under_7/",
name,"_all_sexadj_pruned"))
}
writeit(u,"under_7")
writeit(m,"mid_range")
writeit(o,"over_13")


dogreml<-function(name){
sink(file=paste0("/users/todd/jinshaw/programs/aad/under_7/greml/",name,"_sexadj.sh"))
cat(paste0("/apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 --bfile /well/todd/users/jinshaw/aad/under_7/",name,
"_all_sexadj --autosome --maf 0.01 --make-grm --out /well/todd/users/jinshaw/aad/under_7/grm_",name,"_sexadj --thread-num 10\n"))
cat(paste0("/apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 --grm /well/todd/users/jinshaw/aad/under_7/grm_",name,"_sexadj",
" --pheno /well/todd/users/jinshaw/aad/under_7/pheno_",name,"_sexadj",
" --reml --qcovar /well/todd/users/jinshaw/aad/under_7/covars_",name,"_sexadj",
" --prevalence 0.004 --out /well/todd/users/jinshaw/aad/under_7/outreml_",name,"_sexadj --thread-num 10\n"))
sink()
system(paste0("bash /users/todd/jinshaw/programs/aad/under_7/greml/",name,"_sexadj.sh"))
}
dogreml("under_7")
dogreml("mid_range")
dogreml("over_13")

dogremlprev<-function(name,prev){
sink(file=paste0("/users/todd/jinshaw/programs/aad/under_7/greml/",name,"_",prev,"_sexadj.sh"))
cat(paste0("/apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 --grm /well/todd/users/jinshaw/aad/under_7/grm_",name,"_sexadj",
" --pheno /well/todd/users/jinshaw/aad/under_7/pheno_",name,"_sexadj",
" --reml --qcovar /well/todd/users/jinshaw/aad/under_7/covars_",name,"_sexadj",
" --prevalence ",prev," --out /well/todd/users/jinshaw/aad/under_7/outreml_",name,"_",prev,"_sexadj --thread-num 10\n"))
sink()
system(paste0("bash /users/todd/jinshaw/programs/aad/under_7/greml/",name,"_",prev,"_sexadj.sh"))
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




#now generate table for supplementary in results:

sink(file="/well/todd/users/jinshaw/output/aad/under_7/greml/heritibilities_sexadj.txt")

cat(paste0("Disease prevalence (%) (<7,7-13,>13);<7 including HLA;7-13 including HLA;>13 including HLA;<7 excluding HLA;7-13 excluding HLA;>13 excluding HLA\n"))

addline<-function(prevs,und,mid,old, youngnomhc, midnomhc,oldnomhc){
l<-system(paste0("awk -F \'\t\' \'NR==8\' /well/todd/users/jinshaw/aad/under_7/outreml_",und,"_sexadj.hsq"),intern=T)
l<-strsplit(l,split="\t")
h<-data.frame(h=as.numeric(l[[1]][2]), se=as.numeric(l[[1]][3]))
h$lb<-h$h-(qnorm(0.975)*h$se)
h$ub<-h$h+(qnorm(0.975)*h$se)
l1<-system(paste0("awk -F \'\t\' \'NR==8\' /well/todd/users/jinshaw/aad/under_7/outreml_",mid,"_sexadj.hsq"),intern=T)
l1<-strsplit(l1,split="\t")
h1<-data.frame(h=as.numeric(l1[[1]][2]), se=as.numeric(l1[[1]][3]))
h1$lb<-h1$h-(qnorm(0.975)*h1$se)
h1$ub<-h1$h+(qnorm(0.975)*h1$se)
l2<-system(paste0("awk -F \'\t\' \'NR==8\' /well/todd/users/jinshaw/aad/under_7/outreml_",old,"_sexadj.hsq"),intern=T)
l2<-strsplit(l2,split="\t")
h2<-data.frame(h=as.numeric(l2[[1]][2]), se=as.numeric(l2[[1]][3]))
h2$lb<-h2$h-(qnorm(0.975)*h2$se)
h2$ub<-h2$h+(qnorm(0.975)*h2$se)


l3<-system(paste0("awk -F \'\t\' \'NR==8\' /well/todd/users/jinshaw/aad/under_7/outreml_",youngnomhc,"_sexadj.hsq"),intern=T)
l3<-strsplit(l3,split="\t")
h3<-data.frame(h=as.numeric(l3[[1]][2]), se=as.numeric(l3[[1]][3]))
h3$lb<-h3$h-(qnorm(0.975)*h3$se)
h3$ub<-h3$h+(qnorm(0.975)*h3$se)
l4<-system(paste0("awk -F \'\t\' \'NR==8\' /well/todd/users/jinshaw/aad/under_7/outreml_",midnomhc,"_sexadj.hsq"),intern=T)
l4<-strsplit(l4,split="\t")
h4<-data.frame(h=as.numeric(l4[[1]][2]), se=as.numeric(l4[[1]][3]))
h4$lb<-h4$h-(qnorm(0.975)*h4$se)
h4$ub<-h4$h+(qnorm(0.975)*h4$se)
l5<-system(paste0("awk -F \'\t\' \'NR==8\' /well/todd/users/jinshaw/aad/under_7/outreml_",oldnomhc,"_sexadj.hsq"),intern=T)
l5<-strsplit(l5,split="\t")
h5<-data.frame(h=as.numeric(l5[[1]][2]), se=as.numeric(l5[[1]][3]))
h5$lb<-h5$h-(qnorm(0.975)*h5$se)
h5$ub<-h5$h+(qnorm(0.975)*h5$se)

cat(paste0(prevs,";",round(h$h,digits=3), " (",round(h$lb,digits=3),", ",round(h$ub,digits=3),");",
round(h1$h,digits=3), " (",round(h1$lb,digits=3),", ",round(h1$ub,digits=3),");",
round(h2$h,digits=3), " (",round(h2$lb,digits=3),", ",round(h2$ub,digits=3),");",
round(h3$h,digits=3), " (",round(h3$lb,digits=3),", ",round(h3$ub,digits=3),");",
round(h4$h,digits=3), " (",round(h4$lb,digits=3),", ",round(h4$ub,digits=3),");",
round(h5$h,digits=3), " (",round(h5$lb,digits=3),", ",round(h5$ub,digits=3),")\n"))
}
addline("0.4,0.4,0.4","under_7","mid_range","over_13","under_7nomhc","mid_rangenomhc","over_13nomhc")
addline("0.5,0.3,0.3","under_7_0.005","mid_range_0.003","over_13_0.003","under_7nomhc_0.005","mid_rangenomhc_0.003","over_13nomhc_0.003")
addline("0.5,0.2,0.2","under_7_0.005","mid_range_0.002","over_13_0.002","under_7nomhc_0.005","mid_rangenomhc_0.002","over_13nomhc_0.002")
sink()




