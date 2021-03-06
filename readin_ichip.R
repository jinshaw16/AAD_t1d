#readin_ichip.R

library(snpStats)
library(annotSnpStats)
library(ggplot2)

d<-"/well/todd/users/jinshaw/t1d_risk/immunochip/"
j<-read.table(file=paste0(d,"sanger-controls-preqc.bim"),header=F,as.is=T)
j<-j[j$V5!=0 & j$V6!=0,]

#READ IN EACH GENOTYPE DATASET:
readingeno<-function(base, name){
g<-read.plink(fam=paste0(d, base,".fam"),
                 bed=paste0(d, base,".bed"),
                 bim=paste0(d, base,".bim"))

 
#create annotated SNP matrix
g <- annot.plink(g)
g<-g[,j$V2]
s<-g@snps
w<-which(s$allele.1!=j$V5 & !is.na(s$allele.1))
g<-switch.alleles(g,snps=w)
s<-g@snps
s<-s[!is.na(s$allele.1),]
g<-g[,rownames(s)]
g@samples$member[ is.na(g@samples$member) ] <- 3
rownames(g) <- paste(g@samples$pedigree, g@samples$member, sep=".")
 
save(g, file=paste0(d, name,"_initial.RData"))
return(g)
}
 
m<-mapply(readingeno,base=c("t1d-cases-preqc","finn-preqc", "t1dgc-asp-preqc", "sanger-controls-preqc","uva-controls-preqc"),
name=c("uk","finn", "t1dgc","sanger","uva"), SIMPLIFY=FALSE)

names(m)<-c("uk","finn", "t1dgc","sanger","uva")

s1<-m[[1]]@snps
s2<-m[[2]]@snps
s3<-m[[3]]@snps
s4<-m[[4]]@snps
s5<-m[[5]]@snps

uksnps<-intersect(s1$snp.name,s2$snp.name)
uksnps<-intersect(uksnps,s4$snp.name)
uksnps<-intersect(uksnps,s5$snp.name)

m[["uk"]]<-m[["uk"]][,uksnps]
m[["sanger"]]<-m[["sanger"]][,uksnps]
m[["uva"]]<-m[["uva"]][,uksnps]
m[["finn"]]<-m[["finn"]][,uksnps]

##########################################################################
#combining ancestry groups and obtaining the AAD data etc for each group:#
##########################################################################
#1)UK:
#combine UK + NI cases, controls and family collections together:
uk<-rbind2(m[["uk"]],m[["sanger"]])
uk<-rbind2(uk,m[["uva"]])

#load sample to identify the british collections and merge in:
ni<-read.table(file=paste0(d,'ni-case-to-uva-2012-04-19.tab'),header=T, sep="\t",comment.char="", as.is=T)
ni$uniqueID<-paste0("ni.",ni$X.sampleid)
nic<-read.table(file=paste0(d,'yh-control-to-uva-2012-04-23.tab'),header=T, sep="\t", comment.char="", as.is=T)
nic$uniqueID<-paste0("ni.",nic$sampleid)
nic$onset<-NA
ni<-ni[,c("uniqueID","onset")]
nic<-nic[,c("uniqueID","onset")]
ni<-rbind(ni,nic)
ni$country="NI"
ni$collection="NI GRID"
war<-read.table(file=paste0(d,"warren-to-uva-2011-04-27.tab"),header=T, as.is=T, comment.char="", sep="\t")
war$uniqueID<-paste0(war$familyid,".",war$member)
war<-war[,c("uniqueID","onset")]
war$country="UK"
war$collection="WARREN"
ukall<-rbind(ni,war)

othuk<-m[["finn"]][rownames(m[["finn"]]) %in% ukall$uniqueID,]
othuk<-othuk[,colnames(uk)]

uk<-rbind2(uk,othuk)

#merge in the other info into the sample slot:
uk_ids<-read.table(file=paste0(d, "/ic-t1d-cases-2015-08-06.tab"), header=TRUE, comment.char="")
rownames(uk_ids)<-paste0(uk_ids$X.collection, ".", uk_ids$X.sampleid)
uk_ages<-read.table(file=paste0(d,"/t1d-subject-lookup-2015-08-06.csv"), header=TRUE, sep="\t",comment.char = "")
uk_ages<-merge(uk_ages,uk_ids,by="X.subjectid")
uk_ages$uniqueID<-paste0(uk_ages$X.collection,".",uk_ages$X.sampleid)
uk_ages$collection="GRID"
uk_ages$country="UK"
uk_ages<-uk_ages[,c("uniqueID","onset","country","collection")]
conts<-read.table(file=paste0(d,"control-sample-lookup-2012-05-21.tab"),header=TRUE, sep="\t",comment.char="")
conts$site<-ifelse(conts$site=="sanger","sanger-controls",
ifelse(conts$site=="uva","uva-controls",NA))
conts$uniqueID=paste0(conts$site,".",conts$X.sampleid)
conts$country="UK"
conts$collection="1958BC/UKBS"
conts$onset=NA
conts<-conts[,c("uniqueID","onset","country","collection")]
ukcc<-rbind(uk_ages, conts)
ukall<-rbind(ukcc,ukall)

ukall<-ukall[ukall$uniqueID %in% rownames(uk),]
ukall<-ukall[!duplicated(ukall$uniqueID),]
uk<-uk[rownames(uk) %in% ukall$uniqueID,]
rownames(ukall)<-ukall$uniqueID
ukall<-ukall[rownames(uk),]
uk@samples<-cbind(uk@samples,ukall)
#rename pedigrees so KING doesn't assume all part of one family:
w<-which(uk@samples$pedigree=="t1d-cases")
uk@samples[w,"pedigree"]<-c(10000000001:(nrow(uk@samples[uk@samples$pedigree=="t1d-cases",])+10000000000))
w<-which(uk@samples$pedigree=="ni")
uk@samples[w,"pedigree"]<-c(20000000001:(nrow(uk@samples[uk@samples$pedigree=="ni",])+20000000000))
w<-which(uk@samples$pedigree=="sanger-controls")
uk@samples[w,"pedigree"]<-c(30000000001:(nrow(uk@samples[uk@samples$pedigree=="sanger-controls",])+30000000000))
w<-which(uk@samples$pedigree=="uva-controls")
uk@samples[w,"pedigree"]<-c(40000000001:(nrow(uk@samples[uk@samples$pedigree=="uva-controls",])+40000000000))
save(uk,file=paste0(d,"all_uk_preqc.RData"))

#2) T1DGC:
#treat T1DGC ASPs as one collection for sample QC purposes (makes sense also since genotyped in one batch):
#for consistency with UK, get the country, onset and collection info:
t1<-read.table(file=paste0(d,"T1DGC.2011.03_Resources.csv"),sep=",", header=T, as.is=T)
t1$uniqueID<-paste0(t1$Family.ID,".",t1$Analytic.ID)
t1$onset=t1$Age.of.Onset
t1$onset=ifelse(t1$onset==".",NA,t1$onset)
t1$onset<-as.numeric(t1$onset)
t1$country<-t1$Cohort
t1$country<-ifelse(is.na(t1$country),"USA",t1$country)
t1$collection="T1DGC"
t1<-t1[t1$uniqueID %in% rownames(m[["t1dgc"]]),]
m[["t1dgc"]]<-m[["t1dgc"]][rownames(m[["t1dgc"]]) %in% t1$uniqueID,]
t1<-t1[,c("uniqueID","onset","country","collection")]
rownames(t1)<-t1$uniqueID
t1<-t1[rownames(m[["t1dgc"]]),]
t1dgc<-m[["t1dgc"]]
t1dgc@samples<-cbind(t1dgc@samples,t1)
save(t1dgc,file=paste0(d,"all_t1dgc_preqc.RData"))

#3) FINNS:
#milwaukee collection - older onset T1Ds from finland
milwauki_ages<-read.table(file=paste0(d,"milwaukee-to-uva-2012-04-26.tab"), header=TRUE, comment.char = "", as.is=T)
milwauki_ages$uniqueID<-paste0("finn.",milwauki_ages$X.subjectid)
milwauki_ages$country="Finland"
milwauki_ages<-milwauki_ages[,c("uniqueID","onset","country","collection")]
#IDDMGEN:
finn_ages<-read.table(paste0(d,"fin-2016-01-21.csv"), header=TRUE, comment.char="")
finn_ages$country="Finland"
finn_ages$collection="IDDMGEN"
finn_ages<-finn_ages[,c("uniqueID","onset","country","collection")]
#T1DGEN:
t1dgen_ages<-read.table(file=paste0(d,"/t1dgen-to-uva-2011-12-01.tab"), header=TRUE, comment.char="")
t1dgen_ages$uniqueID<-paste0("finn.", t1dgen_ages$X.subjectid)
t1dgen_ages$country="Finland"
t1dgen_ages$collection="T1DGEN"
t1dgen_ages<-t1dgen_ages[,c("uniqueID","onset","country","collection")]
finns<-rbind(milwauki_ages, finn_ages)
finns<-rbind(finns,t1dgen_ages)
finns<-finns[finns$uniqueID %in% rownames(m[["finn"]]),]

finn<-m[["finn"]]
finn<-finn[rownames(finn) %in% finns$uniqueID,]
rownames(finns)<-finns$uniqueID
finns<-finns[rownames(finn),]
finn@samples<-cbind(finn@samples, finns)
w<-which(finn@samples$pedigree=="finn")
finn@samples[w,"pedigree"]<-c(60000000001:(nrow(finn@samples[finn@samples$pedigree=="finn",])+60000000000))
save(finn, file=paste0(d,"all_finns_preqc.RData"))

#############################################################################################################
#Per sample QC: exclude dodgy looking samples before reading data into plink and performing standard SNP QC:#
#############################################################################################################
qcit<-function(name, xhomlim, hetmin, hetmax){
gen<-get(load(paste0(d,"all_",name,"_preqc.RData")))
gen<-gen[rownames(gen) %in% rownames(gen@samples[!is.na(gen@samples$affected),]),]
all<-gen@samples

#1) Sex discordance
support<-gen@snps
xchrom<-support[support$chromosome==23,]
xchromsnps<-gen[,colnames(gen) %in% xchrom$snp.name]

discordance<-row.summary(xchromsnps)
discordance$uniqueID<-rownames(discordance)
discordance<-merge(discordance, all, by="uniqueID")
library(ggplot2)
discordance$homozygosity<-1-discordance$Heterozygosity
discordance1<-discordance[!is.na(discordance$sex),]
sexdisc<-ggplot(data=discordance1, aes(as.factor(sex), homozygosity)) + geom_boxplot() + coord_cartesian(ylim=c(0.5,1)) +
  geom_segment(x=0, xend=3, y=xhomlim, yend=xhomlim, color="Red", linetype="dashed") + scale_x_discrete(name="Sex", breaks=c(1,2), labels=c("Male", "Female")) +
  scale_y_continuous(name="Homozygosity")
#cutting off at xhomlim - any Females that are >xhomlim homozygous are discounted, and any Males <xhomlim are discounted
discordance1$drop<-ifelse(discordance1$sex==1 & discordance1$homozygosity<xhomlim,1,
                         ifelse(discordance1$sex==2 & discordance1$homozygosity>xhomlim, 1, 0))
discordance2<-discordance1[discordance1$drop==1,]
drop<-discordance2$uniqueID


#create a file which will be added to cumulatively of sample exclusions
sample.exclusions<-all[all$uniqueID %in% drop,]
sample.exclusions<-data.frame(uniqueID=sample.exclusions[,c("uniqueID")])
sample.exclusions$sex_discord<-TRUE



#2) Missingness/heteorzygosity of genotype data:
#want to assess heterozygosity without including sex chromosomes:
nonsex<-support[support$chromosome<23,]
gnonsex<-gen[,colnames(gen) %in% nonsex$snp.name]
hets<-row.summary(gnonsex)

#Missingness
#we will drop any samples with >0.1 missingness:
hets$missing<-1-hets$Call.rate
hets$uniqueID<-rownames(hets)
sample.exclusions.call<-subset(hets,hets$missing>0.1)
sample.exclusions.call$lowcall=TRUE
sample.exclusions.call<-sample.exclusions.call[,c("uniqueID", "lowcall")]
sample.exclusions<-merge(sample.exclusions, sample.exclusions.call, by="uniqueID", all=TRUE)


#Heterozygosity
#showing plot of missingness vs hets to obtain reasonable hets cut-offs:
misshet<-ggplot(data=hets, aes(Heterozygosity, missing)) + geom_point() +
  geom_segment(x=0,xend=1, y=0.1,yend=0.1, color="red", linetype="dashed") +
  geom_segment(x=hetmin,xend=hetmin, y=0,yend=1, color="red", linetype="dashed") +
  geom_segment(x=hetmax,xend=hetmax, y=0,yend=1, color="red", linetype="dashed") +
  scale_x_continuous(name="Heterozygosity Rate", breaks=c(seq(0.1,0.7,0.1))) +
  scale_y_continuous(name="Missing call Rate", breaks=c(seq(0.1,0.7,0.1))) +
  coord_cartesian(xlim=c(0.1, 0.7), ylim=c(0, 0.75))


sample.exclusions.hets<-subset(hets,hets$Heterozygosity>hetmax|hets$Heterozygosity<hetmin)
sample.exclusions.hets$het=TRUE
sample.exclusions.hets<-sample.exclusions.hets[,c("uniqueID", "het")]
sample.exclusions<-merge(sample.exclusions, sample.exclusions.hets, by="uniqueID", all=TRUE)


#do the per-marker QC on only patients surviving the per-individual QC:
gen<-gen[!(rownames(gen) %in% sample.exclusions$uniqueID),]
#finally fill in the sex data if obvious from genotype data given the homozygocity rate on the x chromosome:
rownames(discordance)<-discordance$uniqueID
discordance<-discordance[rownames(gen),]
w<-which(is.na(gen@samples$sex))
gen@samples$sex[w]<-ifelse(discordance[w,"homozygosity"]>0.65 & discordance[w,"homozygosity"]<0.85,2,
ifelse(discordance[w,"homozygosity"]>0.93,1,NA))

write.plink(file.base=paste0(d,name,"_post_sampqc"),
snps=as(gen,"SnpMatrix"),
pedigree=gen@samples$pedigree,
id=gen@samples$member,
father=gen@samples$father,
mother=gen@samples$mother,
sex=gen@samples$sex,
phenotype=gen@samples$affected,
chromosome=gen@snps$chromosome,
genetic.distance=c(rep(NA,nrow(gen@snps))),
position=gen@snps$position,
allele.1=gen@snps$allele.1,
allele.2=gen@snps$allele.2)



##################
#B) Per-marker QC#
##################
#Doing this in PLINK:
system(paste0("plink --bfile ",d,name,"_post_sampqc --keep-allele-order --maf 0.005 --geno 0.95 --hwe 0.00005 --make-bed --out ",d,name,"_post_qc"))
message(paste0("Done ",name))
}

qcit(name="uk", xhomlim=0.93, hetmin=0.2,hetmax=0.3)
qcit(name="t1dgc", xhomlim=0.93, hetmin=0.18,hetmax=0.3)
qcit(name="finns", xhomlim=0.93, hetmin=0.2,hetmax=0.3)

#############################
#REMOVE RELATED INDIVIDUALS:#
#############################

#generating bash scripts to prune, remove MHC, then perform relationship inference on the remaining iChip data. 
#removing all related individuals to 3rd degree (randomly, but I think KING leans towards dropping cases)

aaddir<-"/well/todd/users/jinshaw/aad/"
calcgenetic<-function(cohort){
sink(file=paste0("~/programs/aad/under_7/relatescripts/relatednesscalc_",cohort,"_nomhc.sh"))
cat(paste0("#relatednesscalc_",cohort,"_nomhc.sh\n"))
cat(paste0("~/software/plink2 --bfile ",d,cohort,"_post_qc",
" --exclude range /well/todd/users/jinshaw/aad/under_7/mhc.txt --make-bed --out ",aaddir,cohort,"_post_qc_nomhc --double-id\n"))
cat(paste0("plink --bfile ",aaddir,cohort,"_post_qc_nomhc",
" --indep-pairwise 1000 50 0.2 --out ",aaddir,cohort,"_p\n"))
cat(paste0("plink --bfile ",aaddir,cohort,"_post_qc_nomhc",
" --exclude ",aaddir,cohort,"_p.prune.out --keep-allele-order --make-bed --out ",aaddir,"/geno_",cohort,"_n\n"))

#remove relateds:
cat(paste0("~/software/king -b ",aaddir,"/geno_",cohort,"_n.bed --unrelated --prefix ",aaddir,"qc/",cohort,"_n\n"))
cat(paste0("plink --bfile ",d,cohort,"_post_qc --keep-allele-order --keep ",aaddir,"qc/",cohort,
"_nunrelated.txt --make-bed --out ",aaddir,"geno_unrel_",cohort,"_n\n"))
sink()
system(paste0("chmod a=rwx ~/programs/aad/under_7/relatescripts/relatednesscalc_",cohort,"_nomhc.sh"))
system(paste0("bash ~/programs/aad/under_7/relatescripts/relatednesscalc_",cohort,"_nomhc.sh"))
}
invisible(lapply(c("uk","finns","t1dgc"),calcgenetic))



#Now combine all individuals (intersect of all SNPs):
s1<-read.table(file=paste0(aaddir,"geno_unrel_uk_n.bim"), header=F,as.is=T)
s2<-read.table(file=paste0(aaddir,"geno_unrel_t1dgc_n.bim"), header=F,as.is=T)
s3<-read.table(file=paste0(aaddir,"geno_unrel_finns_n.bim"), header=F,as.is=T)
snps<-intersect(s1$V2, s2$V2)
snps<-intersect(snps, s3$V2)

#filter the R objects for those unrelated invididuals:
filter<-function(name){
g<-get(load(file=paste0(d,"all_",name,"_preqc.RData")))
i<-read.table(file=paste0(aaddir,"geno_unrel_",name,"_n.fam"), header=F,as.is=T)
i$ped<-ifelse(substr(i$V1,1,7)=="1000000","t1d-cases",
ifelse(substr(i$V1,1,7)=="2000000","ni",
ifelse(substr(i$V1,1,7)=="3000000","sanger-controls",
ifelse(substr(i$V1,1,7)=="4000000","uva-controls",
ifelse(substr(i$V1,1,7)=="6000000","finn",i$V1)))))
i$uniqueID<-ifelse(!i$V2 %in% rownames(g),paste0(i$ped,".",i$V2),i$V2)
g<-g[rownames(g) %in% i$uniqueID,colnames(g) %in% snps]
rownames(i)<-i$uniqueID
i<-i[rownames(g),]
g@samples$sex<-i$V5
save(g,file=paste0(d,name,"_unrel_postqc.RData"))
return(g)
}
alls<-lapply(c("uk","t1dgc","finns"),filter)

alls[[2]]<-alls[[2]][,colnames(alls[[1]])]
alls[[3]]<-alls[[3]][,colnames(alls[[1]])]

s1<-alls[[1]]@snps
s2<-alls[[2]]@snps
s3<-alls[[3]]@snps
w<-which(s2$allele.1!=s1$allele.1)
w1<-which(s3$allele.1!=s1$allele.1)
alls[[2]]<-switch.alleles(alls[[2]],snps=w)
alls[[3]]<-switch.alleles(alls[[3]], snps=w)

all<-rbind2(alls[[1]],alls[[2]])
all<-rbind2(all,alls[[3]])


all@samples$affected<-ifelse(is.na(all@samples$affected) & !is.na(all@samples$onset),2,all@samples$affected)
all@samples$onset<-ifelse(all@samples$affected==1,NA,all@samples$onset)
all@samples$onset<-ifelse(all@samples$onset<0 & !is.na(all@samples$onset),NA,all@samples$onset)
all<-all[!(is.na(all@samples$onset) & all@samples$affected==2),]
samples=all@samples
snps=all@snps
write.plink(file.base=paste0(aaddir,"all_ages_3"),
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

sink(file="~/programs/aad/under_7/pccalc.sh")
cat(paste0("#$ -N pcs
#$ -pe shmem 4
#$ -P todd.prjc
#$ -q short.qc\n"))
cat(paste0("/apps/well/plink/1.90b3/plink --bfile ",aaddir,"all_ages_3 --indep-pairwise 1000 50 0.2 --out ",aaddir,"pruned\n"))
cat(paste0("/apps/well/plink/1.90b3/plink --bfile ",aaddir,"all_ages_3 --exclude ",aaddir,"pruned.prune.out",
" --make-bed --keep-allele-order --out ",aaddir,"forpcad_3\n"))
cat(paste0("~/software/plink2 --bfile ",aaddir,"forpcad_3 ",
"--exclude range /well/todd/users/jinshaw/aad/under_7/mhc.txt --make-bed --out ",aaddir,"forpcad_nomhc_3\n"))
#calculate principal components
cat(paste0("~/software/plink2 --bfile ",aaddir,"forpcad_nomhc_3 --allow-no-sex --pca approx 10 --out ",aaddir,"pcs3\n"))
sink()

system("qsub ~/programs/aad/under_7/pccalc.sh")

aaddir<-"/well/todd/users/jinshaw/aad/"
pcs<-read.table(file=paste0(aaddir,"pcs3.eigenvec"),header=T, as.is=T,comment.char="")
rownames(pcs)<-pcs$IID
pcs<-pcs[,c(3:12)]
pcs<-pcs[rownames(all),]
colnames(pcs)<-paste0("PC",c(1:10))

all@samples<-cbind(all@samples, pcs)

save(all,file=paste0(d,"all_inds_unrel_postqc_3.RData"))
