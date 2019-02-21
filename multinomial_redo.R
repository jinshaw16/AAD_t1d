#multinomial_redo.R
library(snpStats)
library(annotSnpStats)
library(snpStatsWriter)
library(humarray)
library(gridExtra)
library(multinomRob)
library(ggplot2)

#read in data:
load(file="/well/todd/users/jinshaw/aad/under_7/all_patients_combined.RData")
a<-allgeno[rownames(allgeno) %in% rownames(allgeno@samples[!is.na(allgeno@samples$t1d),]),]
load(file="/well/todd/users/jinshaw/aad/under_7/all_combined_alligned.RData")
ages<-allgeno@samples

ages<-ages[,c("uniqueID","onset")]
a@samples<-merge(a@samples,ages,by="uniqueID",all.x=T)
rownames(a@samples)<-a@samples$uniqueID
a@samples<-a@samples[rownames(a@.Data),]
#Want to analyse these by ancestry group:
uk<-a[rownames(a) %in% rownames(a@samples[a@samples$country %in% c("UK","NI"),]),]
fin<-a[rownames(a) %in% rownames(a@samples[a@samples$country %in% c("Finland"),]),]
eur<-a[rownames(a) %in% rownames(a@samples[a@samples$country %in% c("EUR"),]),]
usa<-a[rownames(a) %in% rownames(a@samples[a@samples$country %in% c("USA"),]),]
ap<-a[rownames(a) %in% rownames(a@samples[a@samples$country %in% c("AP"),]),]


writeitout<-function(cohort, name){
write.plink(snps=as(cohort,"SnpMatrix"),
file.base=paste0("/well/todd/users/jinshaw/aad/",name,"_all"),
pedigree=cohort@samples$pedigree,
id=rownames(cohort@samples),
father=cohort@samples$father,
mother=cohort@samples$mother,
sex=cohort@samples$sex,
phenotype=cohort@samples$affected,
chromosome=cohort@snps$chromosome,
position=cohort@snps$position,
allele.1=cohort@snps$allele.1,
allele.2=cohort@snps$allele.2)
}
invisible(mapply(writeitout, cohort=list(uk,fin,eur,usa,ap),
name=c("uk","fin","eur","usa","ap")))


#remove related individuals (to 3rd degree):
calcgenetic<-function(cohort){
#create script which converts the vcf to PLINK, filters the PLINK file, removes the original plink file to save space.
sink(file=paste0("~/programs/aad/under_7/relatescripts/relatednesscalc_allages_",cohort,"_nomhc.sh"))
cat(paste0("#relatednesscalc_",cohort,"_nomhc.sh\n"))
cat(paste0("plink --bfile /well/todd/users/jinshaw/aad/",cohort,"_all --maf 0.01 --hwe 0.0005 --geno 0.05 --make-bed --out ",
"/well/todd/users/jinshaw/aad/filtered_",cohort," --double-id --keep-allele-order\n"))
cat(paste0("~/software/plink2 --bfile /well/todd/users/jinshaw/aad/filtered_",cohort,
" --exclude range /well/todd/users/jinshaw/aad/under_7/mhc.txt --make-bed --out /well/todd/users/jinshaw/aad/filtered1_",cohort," --double-id\n"))
cat(paste0("plink --bfile /well/todd/users/jinshaw/aad/filtered1_",cohort,
" --indep-pairwise 1000 50 0.2 --out /well/todd/users/jinshaw/aad/",cohort,"\n"))
cat(paste0("plink --bfile /well/todd/users/jinshaw/aad/filtered1_",cohort,
" --exclude /well/todd/users/jinshaw/aad/",cohort,".prune.out --make-bed --out /well/todd/users/jinshaw/aad/geno_",cohort,"\n"))

#remove relateds:

cat(paste0("~/software/king -b /well/todd/users/jinshaw/aad/geno_",cohort,".bed --related --degree 3 --prefix /well/todd/users/jinshaw/aad/qc/",cohort,"\n"))
#cat(paste0("plink --bfile /well/todd/users/jinshaw/aad/",cohort,"_all --keep /well/todd/users/jinshaw/aad/qc/",cohort,
#"unrelated.txt --make-bed --out /well/todd/users/jinshaw/aad/geno_unrel_",cohort,"\n"))
cat(paste0("rm /well/todd/users/jinshaw/aad/filtered1_",cohort,"*\n"))
sink()
system(paste0("chmod a=rwx ~/programs/aad/under_7/relatescripts/relatednesscalc_allages_",cohort,"_nomhc.sh"))
system(paste0("bash ~/programs/aad/under_7/relatescripts/relatednesscalc_allages_",cohort,"_nomhc.sh"))
}
invisible(lapply(c("uk","fin","eur","usa","ap"),calcgenetic))


#read in the unrelated indivuals, combine and save:
readplink<-function(cohort){
r<-read.plink(fam=paste0("/well/todd/users/jinshaw/aad/geno_unrel_",cohort,".fam"),
bim=paste0("/well/todd/users/jinshaw/aad/geno_unrel_",cohort,".bim"),
bed=paste0("/well/todd/users/jinshaw/aad/geno_unrel_",cohort,".bed"))
r<-annot.plink(r)
return(r)
}
alls<-lapply(c("uk","fin","eur","usa","ap"),readplink)

all<-rbind2(alls[[1]],alls[[2]])
all<-rbind2(all,alls[[3]])
all<-rbind2(all,alls[[4]])
all<-rbind2(all,alls[[5]])

samples<-all@samples
snps<-all@snps

write.plink(file.base="/well/todd/users/jinshaw/aad/under_7/all_ages",
snps=as(all,"SnpMatrix"),
pedigree=samples$pedigree,
id=rownames(samples),
father=samples$father,
mother=samples$mother,
sex=samples$sex,
phenotype=samples$affected,
chromosome=snps$chromosome,
position=snps$position,
allele.1=snps$allele.1,
allele.2=snps$allele.2)


#get pcs
system(paste0("plink --bfile /well/todd/users/jinshaw/aad/under_7/all_ages --indep-pairwise 1000 50 0.2 --out /well/todd/users/jinshaw/aad/under_7/pruned"))
system(paste0("plink --bfile /well/todd/users/jinshaw/aad/under_7/all_ages --exclude /well/todd/users/jinshaw/aad/under_7/pruned.prune.out",
" --make-bed --out /well/todd/users/jinshaw/aad/under_7/forpcad"))
system(paste0("~/software/plink2 --bfile /well/todd/users/jinshaw/aad/under_7/forpcad ",
"--exclude range /well/todd/users/jinshaw/aad/under_7/mhc.txt --make-bed --out /well/todd/users/jinshaw/aad/under_7/forpcad_nomhc"))
system(paste0("plink --bfile /well/todd/users/jinshaw/aad/under_7/forpcad_nomhc --pca --allow-no-sex --out /well/todd/users/jinshaw/aad/under_7/pcad"))


#read this genotype data in and the pcas, then keep only the hits we're interested in:
g<-read.plink(fam="/well/todd/users/jinshaw/aad/under_7/all_ages.fam",
bed="/well/todd/users/jinshaw/aad/under_7/all_ages.bed",
bim="/well/todd/users/jinshaw/aad/under_7/all_ages.bim")
g<-annot.plink(g)
g@samples$uniqueID<-g@samples$member

pcs<-read.table(file="/well/todd/users/jinshaw/aad/under_7/pcad.eigenvec",as.is=T, header=F)
pcs<-pcs[,c(1:12)]
colnames(pcs)<-c("ped","mem",paste0("PC",1:10))
rownames(pcs)<-pcs$mem

#and combine with the aad data to generate the groups:
load(file="/well/todd/users/jinshaw/aad/under_7/all_combined_alligned.RData")
ages<-allgeno@samples
ages<-ages[,c("uniqueID","onset")]
g@samples<-merge(g@samples,ages,by="uniqueID",all.x=T)
rownames(g@samples)<-g@samples$uniqueID
g@samples<-g@samples[rownames(g@.Data),]
#and country:
load(file="/well/todd/users/jinshaw/aad/under_7/all_patients_combined.RData")
samp<-allgeno@samples[,c("uniqueID","cohort","country")]
g@samples<-merge(g@samples,samp, by="uniqueID", all.x=T)
rownames(g@samples)<-g@samples$uniqueID
g@samples<-g@samples[rownames(g),]

g@samples$group<-ifelse(g@samples$affected==1,0,
ifelse(g@samples$onset<7 & !is.na(g@samples$onset),1,
ifelse(g@samples$onset>=13 & !is.na(g@samples$onset),2,NA)))
pheno<-g@samples
pcs<-pcs[rownames(pheno),]
pheno<-cbind(pheno,pcs)

#Define the 44 T1D regions we want to examine:
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
p<-t1dsnps[!t1dsnps$id %in% colnames(g),]

#get those that we've had to impute to get the index SNP:
imputed<-c("rs7517998","imm_2_162790641",
"imm_6_126795257","imm_7_50995894","imm_10_6170083",
"imm_11_2141424","imm_12_110346958","imm_19_53906282","imm_21_42698791")


#keep only the ones we want to test:
g<-g[,colnames(g) %in% t1dsnps$id]
g<-as(g,"SnpMatrix")
cs<-col.summary(g)
w<-which(cs$RAF>0.5)
g<-switch.alleles(g,snps=w)
b<-as(g,"numeric")

#and the imputed ones:
getsnp<-function(snp,snp1){
load(file=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/",snp,"_clean.RData"))
if(snp=="imm_12_110346958"){
colnames(DATA)<-snp1
}
cs<-col.summary(DATA)
w<-which(cs$RAF>0.5)
DATA<-switch.alleles(DATA,snps=w)
d<-as(DATA[,grepl(snp1,colnames(DATA))],"numeric")
b<<-cbind2(b,d)
return(d)
}
l<-mapply(getsnp, snp=imputed, snp1=p$snp,SIMPLIFY=F)
colnames(b)<-gsub(":.*","",colnames(b))

table(rownames(b)==rownames(pheno))
pheno<-cbind(pheno,b)
p1<-pheno
save(b,pheno,t1dsnps, file="/well/todd/users/jinshaw/aad/under_7/pheno_mult.R")

load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult.R")

#now perform multinomial regression (remove 7-13s):
pheno<-pheno[pheno$group %in% c(0,1,2) & !is.na(pheno$group),]
pheno$g0<-ifelse(pheno$group==0,1,ifelse(pheno$group!=0,0,NA))
pheno$g1<-ifelse(pheno$group==1,1,ifelse(pheno$group!=1,0,NA))
pheno$g2<-ifelse(pheno$group==2,1,ifelse(pheno$group!=2,0,NA))

t1dsnps$altid<-ifelse(substr(t1dsnps$id,1,1)=="1",paste0("X",t1dsnps$id),t1dsnps$id)
colnames(pheno)<-ifelse(substr(colnames(pheno),1,1)=="1",paste0("X",colnames(pheno)),colnames(pheno))
t1dsnps$altid<-ifelse(!t1dsnps$altid %in% colnames(pheno),t1dsnps$snp,t1dsnps$altid)


#carry out the multinomial regressions:
getlikelihoods<-function(snpname){
#model 1: allow the beta for the SNP to vary for both:
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))

one<-multinomRob(model=list(form0,form1,form2),data=pheno, print.level=1, MLEonly=T)
llk1<-one$value

equals0<-as.formula(paste0("g1 ~ `",snpname,"` + 0"))
equals1<-as.formula(paste0("g2 ~ `",snpname,"` + 0"))
#model 2: constrain betas to equal each other:
two<-multinomRob(model=list(form0,form1,form2),
equality=list(list(equals0,equals1)), data=pheno, print.level=1, MLEonly=T)
llk2<-two$value
l<-data.frame(unconstrained=llk1, constrained=llk2, logor1=one$coefficients[12,2], logse1=one$se[12,2],
logor2=one$coefficients[12,3], logse2=one$se[12,3])
l$constrained=l$constrained*-1
l$unconstrained=l$unconstrained*-1
l$loglambda<-l$constrained-l$unconstrained
#compare using likelihood ratio test:
l$chisq<-l$loglambda*-2
l$p<-pchisq(l$chisq,1, lower.tail=F)
l$lb1<-l$logor1-(qnorm(0.975)*l$logse1)
l$ub1<-l$logor1+(qnorm(0.975)*l$logse1)
l$lb2<-l$logor2-(qnorm(0.975)*l$logse2)
l$ub2<-l$logor2+(qnorm(0.975)*l$logse2)
return(l)
}

likelihoods<-lapply(t1dsnps$altid,getlikelihoods)
likelihoods<-do.call("rbind", likelihoods)
likelihoods$snp=t1dsnps$snp
likelihoods$id<-t1dsnps$id
likelihoods$loci<-t1dsnps$loci
likelihoods$altid<-t1dsnps$altid
write.table(likelihoods, file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo.txt",
col.names=T, row.names=F, quote=F, sep="\t")

likelihoods$logp<-log10(likelihoods$p)*-1
likelihoods<-likelihoods[order(-likelihoods$logp),]
likelihoods$loci<-as.factor(likelihoods$loci)
likelihoods$ord<-c(nrow(likelihoods):1)
likelihoods$loci<-reorder(likelihoods$loci,likelihoods$ord)
r<-likelihoods
one<-ggplot(data=r, aes(x=logor1,y=as.factor(loci))) + geom_point(,colour="red") +
geom_point(data=r, aes(x=logor2,as.numeric(loci)+0.3), colour="blue") +
geom_errorbarh(data=r, aes(xmin=lb1,xmax=ub1, y=as.factor(loci)),colour="red", height=0.01) +
geom_errorbarh(data=r, aes(xmin=lb2,xmax=ub2, y=as.numeric(loci)+0.3),colour="blue", height=0.01) +
geom_vline(xintercept=0, colour="red", linetype="dashed") +
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 7 (red) and over 13 (blue)") +
scale_y_discrete(name="Locus")
two<-ggplot(data=likelihoods, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/56)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity"))))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/inds_het_tests_all.png", res=400,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()

#plotting all 3 odds ratios:
#to do this, fitting a third multinomial logistic regression with the 7-13 group included as an outcome:

load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult.R")
pheno$group<-ifelse(pheno$affected==1,0,
ifelse(pheno$onset<7 & !is.na(pheno$onset),1,
ifelse(pheno$onset>=7 &pheno$onset<13 & !is.na(pheno$onset),2,
ifelse(pheno$onset>=13 & !is.na(pheno$onset),3,NA))))
pheno$g0<-ifelse(pheno$group==0,1,ifelse(pheno$group!=0,0,NA))
pheno$g1<-ifelse(pheno$group==1,1,ifelse(pheno$group!=1,0,NA))          
pheno$g2<-ifelse(pheno$group==2,1,ifelse(pheno$group!=2,0,NA))          
pheno$g3<-ifelse(pheno$group==3,1,ifelse(pheno$group!=3,0,NA))          
colnames(pheno)<-ifelse(substr(colnames(pheno),1,1)=="1",paste0("X",colnames(pheno)),colnames(pheno))

t1dsnps$altid<-ifelse(substr(t1dsnps$id,1,1)=="1",paste0("X",t1dsnps$id),t1dsnps$id)
colnames(pheno)<-ifelse(substr(colnames(pheno),1,1)=="1",paste0("X",colnames(pheno)),colnames(pheno))
t1dsnps$altid<-ifelse(!t1dsnps$altid %in% colnames(pheno),t1dsnps$snp,t1dsnps$altid)


getlikelihoods1<-function(snpname){
snpname<-ifelse(substr(snpname,1,1)=="1",paste0("X",snpname),snpname)
colnames(pheno)[ncol(pheno)]<-ifelse(substr(colnames(pheno)[ncol(pheno)],1,1)=="1",paste0("X",colnames(pheno)[ncol(pheno)]),colnames(pheno)[ncol(pheno)])

if(snpname=="imm_2_162845188"){
pheno[,snpname]<-ifelse(pheno[,snpname]==2,0,ifelse(pheno[,snpname]==0,2,pheno[,snpname]))
)

#model 1: allow the beta for the SNP to vary for both:
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form3<-as.formula(paste0("g3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))

one<-multinomRob(model=list(form0,form1,form2,form3),data=pheno, print.level=1, MLEonly=T)

l<-data.frame(logor1=one$coefficients[12,2], logse1=one$se[12,2],
logor2=one$coefficients[12,3], logse2=one$se[12,3], 
logor3=one$coefficients[12,4],logse3=one$se[12,4])
if(snpname=="imm_2_162845188"){
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
return(l)
}

like1<-lapply(t1dsnps$altid,getlikelihoods1)
like1<-do.call("rbind",like1)
like1$id<-t1dsnps$id
like1$loci<-t1dsnps$loci
like1$snp=t1dsnps$snp
like1$altid<-t1dsnps$altid
rownames(likelihoods)<-likelihoods$altid
likelihoods<-likelihoods[like1$altid,]
#steal likelihoods from initial analysis including just <7s and >13s.
like1$constrained=likelihoods$constrained
like1$unconstrained=likelihoods$unconstrained
like1$loglambda=likelihoods$loglambda
like1$p=likelihoods$p
rownames(likelihoods)<-likelihoods$id
#save:
write.table(likelihoods, file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_3.txt",
col.names=T, row.names=F, quote=F, sep="\t")

#plot results:
likelihoods<-read.table(file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_3.txt",
header=T, as.is=T, sep="\t")
rownames(likelihoods)<-likelihoods$altid
likelihoods$logp<-log10(likelihoods$p)*-1
likelihoods<-likelihoods[order(-likelihoods$logp),]
likelihoods$loci<-as.factor(likelihoods$loci)
likelihoods$ord<-c(nrow(likelihoods):1)
likelihoods$loci<-reorder(likelihoods$loci,likelihoods$ord)
r<-likelihoods
one<-ggplot(data=r, aes(x=logor3,y=as.factor(loci))) + geom_point(,colour="blue") +
geom_point(data=r, aes(x=logor2,as.numeric(loci)+0.2), colour="green") +
geom_point(data=r, aes(x=logor1,as.numeric(loci)+0.4), colour="red") +
geom_errorbarh(data=r, aes(xmin=lb3,xmax=ub3, y=as.factor(loci)),colour="blue", height=0.01) +
geom_errorbarh(data=r, aes(xmin=lb2,xmax=ub2, y=as.numeric(loci)+0.2),colour="green", height=0.01) +
geom_errorbarh(data=r, aes(xmin=lb1,xmax=ub1, y=as.numeric(loci)+0.4),colour="red", height=0.01) +
geom_vline(xintercept=0, colour="red", linetype="dashed") +
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 7 (red), 7-13 (green) and over 13 (blue)") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12),
axis.text.y=element_text(size=12))
two<-ggplot(data=likelihoods, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/56)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
geom_vline(aes(xintercept=1.823909), colour="red") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12), 
axis.text.y=element_text(size=12)) +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <7 and >13"))))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/inds_het_tests_all_inc_midrange.png", res=800,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()

#show all loci with p<0.5 just so it isn't such a busy Figure:
r$pfdr<-p.adjust(r$p, method = "BH")
r1<-r[r$p<0.1,]
r1$ord<-c(nrow(r1):1)
r1$loci<-reorder(r1$loci,r1$ord)
one<-ggplot(data=r1, aes(x=logor3,y=as.factor(loci))) + geom_point(,colour="blue") +
geom_point(data=r1, aes(x=logor2,as.numeric(loci)+0.2), colour="green") +
geom_point(data=r1, aes(x=logor1,as.numeric(loci)+0.4), colour="red") +
geom_errorbarh(data=r1, aes(xmin=lb3,xmax=ub3, y=as.factor(loci)),colour="blue", height=0.01) +
geom_errorbarh(data=r1, aes(xmin=lb2,xmax=ub2, y=as.numeric(loci)+0.2),colour="green", height=0.01) +
geom_errorbarh(data=r1, aes(xmin=lb1,xmax=ub1, y=as.numeric(loci)+0.4),colour="red", height=0.01) +
geom_vline(xintercept=0, colour="red", linetype="dashed") +
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 7 (red), 7-13 (green) and over 13 (blue)") +
scale_y_discrete(name="Locus") + 
theme(axis.title.y=element_text(size=12), 
axis.text.y=element_text(size=12))
two<-ggplot(data=r1, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/55)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
geom_vline(aes(xintercept=1.823909), colour="red") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <7 and >13")))) +
coord_cartesian(xlim=c(0,5)) +
theme(axis.title.y=element_text(size=12), 
axis.text.y=element_text(size=12))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/het_tests_midrange_just0.1.png", res=800,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()


#now only those with a hint of heterogeneity (FDR<0.1):
r$pfdr<-p.adjust(r$p, method = "BH")
r1<-r[r$pfdr<0.1,]
r1$ord<-c(nrow(r1):1)
r1$loci<-reorder(r1$loci,r1$ord)
one<-ggplot(data=r1, aes(x=logor3,y=as.factor(loci))) + geom_point(,colour="blue") +
geom_point(data=r1, aes(x=logor2,as.numeric(loci)+0.2), colour="green") +
geom_point(data=r1, aes(x=logor1,as.numeric(loci)+0.4), colour="red") +
geom_errorbarh(data=r1, aes(xmin=lb3,xmax=ub3, y=as.factor(loci)),colour="blue", height=0.01) +
geom_errorbarh(data=r1, aes(xmin=lb2,xmax=ub2, y=as.numeric(loci)+0.2),colour="green", height=0.01) +
geom_errorbarh(data=r1, aes(xmin=lb1,xmax=ub1, y=as.numeric(loci)+0.4),colour="red", height=0.01) +
geom_vline(xintercept=0, colour="red", linetype="dashed") +
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 7 (red), 7-13 (green) and over 13 (blue)") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=15), 
axis.text.y=element_text(size=12))

two<-ggplot(data=r1, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/55)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=1.823909), colour="red", linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <7 and >13")))) +
coord_cartesian(xlim=c(0,5)) +
theme(axis.title.y=element_text(size=15), 
axis.text.y=element_text(size=12))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/het_tests_all_inc_midrange_no_nonsig.png", res=800,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()




########################
#NOW FOR THE MAFS PLOTS#
########################
#read in all the unrelateds - including the mid range individuals
load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult.R")
t1dsnps$loci<-as.character(t1dsnps$loci)
t1dsnps$loci<-ifelse(t1dsnps$loci=="FASLG","TNFSF4",t1dsnps$loci)
t1dsnps$loci<-ifelse(t1dsnps$loci=="BCAR1","CTRB1",t1dsnps$loci)
#get most likely call for the variants:
g<-as(b,"SnpMatrix")
load(file="/well/todd/users/jinshaw/aad/under_7/all_combined_alligned.RData")
all<-g

a<-as(all,"numeric")
samp<-pheno[,c(1:23)]
samp<-cbind(samp,a)
cases<-samp[samp$affected==2,]
cases$logage<-log(cases$onset)



#plots by allele:
library(plyr)
library(dplyr)
hits<-c("imm_1_170947654", "imm_6_128335625", "imm_15_77022012",
"imm_17_35306733", "imm_9_4280823",   "imm_1_205006527","imm_16_73809828",
"imm_10_6170083")

plotagebysnp<-function(snp){
case<-cases[!is.na(cases[,snp]),]

g<-case %>%
  group_by(case[,snp]) %>%
  summarise(mean = mean(logage,na.rm=T), var=var(logage,na.rm=T), N = n())

g$se<-sqrt(g$var/g$N)
g$lb<-g$mean-(qnorm(0.975)*g$se)
g$ub<-g$mean+(qnorm(0.975)*g$se)
g<-as.data.frame(g)

d<-ggplot(data=case, aes(x = factor(case[,snp]), y = logage))+
    #geom_dotplot(binaxis = 'y', stackdir = 'center',
    #             position = position_dodge(),dotsize=0.01) +
   geom_point(data=g, aes(as.factor(g[,1]), mean), colour="red") +
   geom_segment(data=g, aes(x=as.factor(g[,1]),xend=as.factor(g[,1]),y=lb, yend=ub), colour="red") +
   scale_y_continuous(name="log(AAD)") +
   scale_x_discrete(name=paste0("Copies of minor allele at ",t1dsnps[t1dsnps$id==snp,"snp"]," (",t1dsnps[t1dsnps$id==snp,"loci"],")")) +
   coord_cartesian(ylim=c(1.85,2.15)) +
   annotate("text", x=as.numeric(g[,1])+1.2, y=2.15, label=paste0("N=",g$N))
ggsave(d,file=paste0("/well/todd/users/jinshaw/output/aad/under_7/mafs/redo/",gsub("/"," ",t1dsnps[t1dsnps$id==snp,"loci"]),"_by_allele.png"),
dpi=400, height=20, width=20, units="cm")
}
lapply(hits, plotagebysnp)

samp$onset<-ifelse(samp$affected==1,0,samp$onset)

samp$onset<-ifelse(samp$affected==1,0,
ifelse(samp$onset<1 & samp$affected==2,1,
ifelse(samp$onset>1 & samp$onset<=2,2,
ifelse(samp$onset>2 & samp$onset<=3,3,
ifelse(samp$onset>3 & samp$onset<=4,4,
ifelse(samp$onset>4 & samp$onset<=5,5,
ifelse(samp$onset>5 & samp$onset<=6,6,
ifelse(samp$onset>6 & samp$onset<=7,7,
ifelse(samp$onset>7 & samp$onset<=8,8,
ifelse(samp$onset>8 & samp$onset<=9,9,
ifelse(samp$onset>9 & samp$onset<=10,10,
ifelse(samp$onset>10 & samp$onset<=11,11,
ifelse(samp$onset>11 & samp$onset<=12,12,
ifelse(samp$onset>12 & samp$onset<=13,13,
ifelse(samp$onset>13 & samp$onset<=14,14,
ifelse(samp$onset>14 & samp$onset<=15,15,
ifelse(samp$onset>15 & samp$onset<=16,16,
ifelse(samp$onset>16 & samp$onset<=17,17,
ifelse(samp$onset>17 & samp$onset<=18,18,samp$onset)))))))))))))))))))
samp$onset<-ceiling(samp$onset)
rownames(samp)<-samp$uniqueID
samp<-samp[rownames(all),]

cs<-col.summary(all)
w<-which(cs$RAF>0.5)
all<-switch.alleles(all, snps=w)

samp$onset<-ifelse(samp$onset>=16 & samp$onset<=20,18,
ifelse(samp$onset>20 & samp$onset<=25, 23,
ifelse(samp$onset>25 & samp$onset<=30, 28,
ifelse(samp$onset>30 & samp$onset<=35, 33,
ifelse(samp$onset>35, 37,samp$onset)))))
t<-names(table(samp$onset))

getmafsall<-function(snp){
getcs<-function(onset){
if(onset==0){
a<-all[rownames(all) %in% rownames(samp[samp$affected==1,]),snp]
}
if(onset!=0){
a<-all[rownames(all) %in% rownames(samp[samp$affected==2 & samp$onset==onset,]),snp]
}
cs<-col.summary(a)
rs<-row.summary(a)
out<-data.frame(snp=snp,onset=onset, maf=cs$MAF)
out$se<-sqrt(((out$maf)*(1-out$maf))/(2*nrow(rs)))
out$lb<-out$maf-(qnorm(0.975)*out$se)
out$ub<-out$maf+(qnorm(0.975)*out$se)
return(out)
}
l1<-lapply(t,getcs)
l1<-do.call("rbind", l1)

l1$onset<-as(l1$onset,"character")
l1$onset<-as(l1$onset, "numeric")
g<-ggplot(data=l1, aes(x=onset, y=maf)) + geom_point() +
geom_errorbar(data=l1, aes(x=onset, ymin=lb, ymax=ub)) +
labs(title=paste0(t1dsnps[t1dsnps$id==snp,"loci"])) +
scale_y_continuous(name="Minor Allele Frequency") + theme(axis.text.x=element_text(angle=90)) +
scale_x_continuous(name="Age-at-diagnosis", breaks=c(as.numeric(t)),
labels=c("Controls",t[2]:t[length(t)-5],"16-20","20-25","25-30","30-35",">35")) +
geom_hline(yintercept=l1[l1$onset==0,"maf"],colour="red", linetype="dashed")  +
geom_vline(xintercept=6.5,colour="red", linetype="dashed") +
geom_vline(xintercept=12.5,colour="red", linetype="dashed")

ggsave(g,file=paste0("/well/todd/users/jinshaw/output/aad/under_7/mafs/redo/",gsub("/"," ",t1dsnps[t1dsnps$id==snp,"loci"]),"_all.png"),
dpi=400, height=20, width=20, units="cm")
}

lapply(hits,getmafsall)

