#multinomial_sensitivity_to_cutoff.R
library(snpStats)
library(annotSnpStats)
library(snpStatsWriter)
library(humarray)
library(gridExtra)
library(multinomRob)
library(ggplot2)


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

#load data:
load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult.R")

#now perform multinomial regression of <6 compared to >13:
pheno$group<-ifelse(pheno$onset<6 & !is.na(pheno$onset) & pheno$affected==2,1,
ifelse(pheno$onset>=6 & pheno$onset<13 & !is.na(pheno$onset) & pheno$affected==2,NA,
ifelse(pheno$onset>=13 & !is.na(pheno$onset)	& pheno$affected==2,2,
ifelse(pheno$affected==1,0,NA)))) 

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
return(l)
}

likelihoods<-lapply(t1dsnps$altid,getlikelihoods)
likelihoods<-do.call("rbind", likelihoods)
likelihoods$snp=t1dsnps$snp
likelihoods$id<-t1dsnps$id
likelihoods$loci<-t1dsnps$loci
write.table(likelihoods, file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_lessthan_6.txt",
col.names=T, row.names=F, quote=F, sep="\t")

likelihoods<-read.table(file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_lessthan_6.txt",
header=T, as.is=T, sep="\t")
likelihoods$altid<-t1dsnps$altid
likelihoods$constrained=likelihoods$constrained*-1
likelihoods$unconstrained=likelihoods$unconstrained*-1
likelihoods$loglambda<-likelihoods$constrained-likelihoods$unconstrained
likelihoods$chisq<-likelihoods$loglambda*-2
likelihoods$p<-pchisq(likelihoods$chisq,1, lower.tail=F)
rownames(likelihoods)<-likelihoods$altid
likelihoods$loci<-ifelse(likelihoods$loci=="FASLG","TNFSF4",likelihoods$loci)
likelihoods$loci<-ifelse(likelihoods$loci=="BCAR1","CTRB1",likelihoods$loci)

likelihoods$lb1<-likelihoods$logor1-(qnorm(0.975)*likelihoods$logse1)
likelihoods$ub1<-likelihoods$logor1+(qnorm(0.975)*likelihoods$logse1)
likelihoods$lb2<-likelihoods$logor2-(qnorm(0.975)*likelihoods$logse2)
likelihoods$ub2<-likelihoods$logor2+(qnorm(0.975)*likelihoods$logse2)
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
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 6 (red) and over 13 (blue)") +
scale_y_discrete(name="Locus")
two<-ggplot(data=likelihoods, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/55)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity"))))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/inds_het_tests_all_lessthan_6.png", res=400,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()

#plotting all 3 odds ratios:

load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult.R")
pheno$group<-ifelse(pheno$affected==1,0,
ifelse(pheno$onset<6 & !is.na(pheno$onset),1,
ifelse(pheno$onset>=6 &pheno$onset<13 & !is.na(pheno$onset),2,
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

if(snpname %in% c("imm_2_162845188","seq-rs35667974")){
pheno[,snpname]<-ifelse(pheno[,snpname]==2,0,ifelse(pheno[,snpname]==0,2,pheno[,snpname]))
}

#model 1: allow the beta for the SNP to vary for both:
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form3<-as.formula(paste0("g3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))

one<-multinomRob(model=list(form0,form1,form2,form3),data=pheno, print.level=1, MLEonly=T)

l<-data.frame(logor1=one$coefficients[12,2], logse1=one$se[12,2],
logor2=one$coefficients[12,3], logse2=one$se[12,3], 
logor3=one$coefficients[12,4],logse3=one$se[12,4])
if(snpname %in% c("imm_2_162845188", "seq-rs35667974")){
l$logor1<-l$logor1*-1
l$logor2<-l$logor2*-1
l$logor3<-l$logor3*-1
}
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
#steal likelihoods from initial analysis including just <6s and >13s.
like1$constrained=likelihoods$constrained
like1$unconstrained=likelihoods$unconstrained
likelihoods<-like1
likelihoods$loglambda<-likelihoods$constrained-likelihoods$unconstrained
likelihoods$chisq<-likelihoods$loglambda*-2
likelihoods$p<-pchisq(likelihoods$chisq,1, lower.tail=F)
rownames(likelihoods)<-likelihoods$id
write.table(likelihoods, file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_3_lessthan_6.txt",
col.names=T, row.names=F, quote=F, sep="\t")

likelihoods<-read.table(file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_3_lessthan_6.txt",
header=T, as.is=T, sep="\t")
rownames(likelihoods)<-likelihoods$altid
likelihoods$loci<-ifelse(likelihoods$loci=="FASLG","TNFSF4",likelihoods$loci)
likelihoods$loci<-ifelse(likelihoods$loci=="BCAR1","CTRB1",likelihoods$loci)

likelihoods$lb1<-likelihoods$logor1-(qnorm(0.975)*likelihoods$logse1)
likelihoods$ub1<-likelihoods$logor1+(qnorm(0.975)*likelihoods$logse1)
likelihoods$lb2<-likelihoods$logor2-(qnorm(0.975)*likelihoods$logse2)
likelihoods$ub2<-likelihoods$logor2+(qnorm(0.975)*likelihoods$logse2)
likelihoods$lb3<-likelihoods$logor3-(qnorm(0.975)*likelihoods$logse3)
likelihoods$ub3<-likelihoods$logor3+(qnorm(0.975)*likelihoods$logse3)
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
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 6 (red),6-13 (green) and over 13 (blue)") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12),
axis.text.y=element_text(size=12))
two<-ggplot(data=likelihoods, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/56)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
geom_vline(aes(xintercept=1.9), colour="red") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12), 
axis.text.y=element_text(size=12)) +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <6 and >13"))))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/inds_het_tests_all_inc_midrange_lessthan_6.png", res=800,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()

#show all loci with p<0.5 just so it isn't such a busy Figure:
r$pfdr<-p.adjust(r$p, method = "BH")
r1<-r[r$p<0.1,]
r1$ord<-c(nrow(r1):1)
r1$loci<-reorder(r1$loci,r1$ord)
o<-r1[r1$pfdr<0.1,]
o<-o[o$pfdr==max(o$pfdr),]
one<-ggplot(data=r1, aes(x=logor3,y=as.factor(loci))) + geom_point(,colour="blue") +
geom_point(data=r1, aes(x=logor2,as.numeric(loci)+0.2), colour="green") +
geom_point(data=r1, aes(x=logor1,as.numeric(loci)+0.4), colour="red") +
geom_errorbarh(data=r1, aes(xmin=lb3,xmax=ub3, y=as.factor(loci)),colour="blue", height=0.01) +
geom_errorbarh(data=r1, aes(xmin=lb2,xmax=ub2, y=as.numeric(loci)+0.2),colour="green", height=0.01) +
geom_errorbarh(data=r1, aes(xmin=lb1,xmax=ub1, y=as.numeric(loci)+0.4),colour="red", height=0.01) +
geom_vline(xintercept=0, colour="red", linetype="dashed") +
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 6 (red),6-13 (green) and over 13 (blue)") +
scale_y_discrete(name="Locus") + 
theme(axis.title.y=element_text(size=12), 
axis.text.y=element_text(size=12))
two<-ggplot(data=r1, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/55)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
geom_vline(aes(xintercept=1.9), colour="red") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <6 and >13")))) +
coord_cartesian(xlim=c(0,5)) +
theme(axis.title.y=element_text(size=12), 
axis.text.y=element_text(size=12))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/het_tests_midrange_just0.1_lessthan_6.png", res=800,
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
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 6 (red),6-13 (green) and over 13 (blue)") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=15), 
axis.text.y=element_text(size=12))

two<-ggplot(data=r1, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/55)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=1.9), colour="red", linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <6 and >13")))) +
coord_cartesian(xlim=c(0,5)) +
theme(axis.title.y=element_text(size=15), 
axis.text.y=element_text(size=12))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/het_tests_all_inc_midrange_no_nonsig_lessthan_6.png", res=800,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()



########################
#SAME BUT FOR <5 VS >13#
########################

#load data:
load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult.R")

#now perform multinomial regression of <5 compared to >13:
pheno$group<-ifelse(pheno$onset<5 & !is.na(pheno$onset) & pheno$affected==2,1,
ifelse(pheno$onset>=5 & pheno$onset<13 & !is.na(pheno$onset) & pheno$affected==2,NA,
ifelse(pheno$onset>=13 & !is.na(pheno$onset)    & pheno$affected==2,2,
ifelse(pheno$affected==1,0,NA))))

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
return(l)
}

likelihoods<-lapply(t1dsnps$altid,getlikelihoods)
likelihoods<-do.call("rbind", likelihoods)
likelihoods$snp=t1dsnps$snp
likelihoods$id<-t1dsnps$id
likelihoods$loci<-t1dsnps$loci
write.table(likelihoods, file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_lessthan_5.txt",
col.names=T, row.names=F, quote=F, sep="\t")

likelihoods<-read.table(file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_lessthan_5.txt",
header=T, as.is=T, sep="\t")
likelihoods$altid<-t1dsnps$altid
likelihoods$constrained=likelihoods$constrained*-1
likelihoods$unconstrained=likelihoods$unconstrained*-1
likelihoods$loglambda<-likelihoods$constrained-likelihoods$unconstrained
likelihoods$chisq<-likelihoods$loglambda*-2
likelihoods$p<-pchisq(likelihoods$chisq,1, lower.tail=F)
rownames(likelihoods)<-likelihoods$altid
likelihoods$loci<-ifelse(likelihoods$loci=="FASLG","TNFSF4",likelihoods$loci)
likelihoods$loci<-ifelse(likelihoods$loci=="BCAR1","CTRB1",likelihoods$loci)

likelihoods$lb1<-likelihoods$logor1-(qnorm(0.975)*likelihoods$logse1)
likelihoods$ub1<-likelihoods$logor1+(qnorm(0.975)*likelihoods$logse1)
likelihoods$lb2<-likelihoods$logor2-(qnorm(0.975)*likelihoods$logse2)
likelihoods$ub2<-likelihoods$logor2+(qnorm(0.975)*likelihoods$logse2)
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
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 5 (red) and over 13 (blue)") +
scale_y_discrete(name="Locus")
two<-ggplot(data=likelihoods, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/55)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity"))))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/inds_het_tests_all_lessthan_5.png", res=400,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()

#plotting all 3 odds ratios:
load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult.R")
pheno$group<-ifelse(pheno$affected==1,0,
ifelse(pheno$onset<5 & !is.na(pheno$onset),1,
ifelse(pheno$onset>=5 &pheno$onset<13 & !is.na(pheno$onset),2,
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

if(snpname %in% c("imm_2_162845188","seq-rs35667974")){
pheno[,snpname]<-ifelse(pheno[,snpname]==2,0,ifelse(pheno[,snpname]==0,2,pheno[,snpname]))
}

#model 1: allow the beta for the SNP to vary for both:
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form3<-as.formula(paste0("g3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))

one<-multinomRob(model=list(form0,form1,form2,form3),data=pheno, print.level=1, MLEonly=T)

l<-data.frame(logor1=one$coefficients[12,2], logse1=one$se[12,2],
logor2=one$coefficients[12,3], logse2=one$se[12,3], 
logor3=one$coefficients[12,4],logse3=one$se[12,4])
if(snpname %in% c("imm_2_162845188","seq-rs35667974")){
l$logor1<-l$logor1*-1
l$logor2<-l$logor2*-1
l$logor3<-l$logor3*-1
}
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
#steal likelihoods from initial analysis including just <5s and >13s.
like1$constrained=likelihoods$constrained
like1$unconstrained=likelihoods$unconstrained
likelihoods<-like1
likelihoods$loglambda<-likelihoods$constrained-likelihoods$unconstrained
likelihoods$chisq<-likelihoods$loglambda*-2
likelihoods$p<-pchisq(likelihoods$chisq,1, lower.tail=F)
rownames(likelihoods)<-likelihoods$id
write.table(likelihoods, file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_3_lessthan_5.txt",
col.names=T, row.names=F, quote=F, sep="\t")

likelihoods<-read.table(file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_3_lessthan_5.txt",
header=T, as.is=T, sep="\t")
rownames(likelihoods)<-likelihoods$altid
likelihoods$loci<-ifelse(likelihoods$loci=="FASLG","TNFSF4",likelihoods$loci)
likelihoods$loci<-ifelse(likelihoods$loci=="BCAR1","CTRB1",likelihoods$loci)

likelihoods$lb1<-likelihoods$logor1-(qnorm(0.975)*likelihoods$logse1)
likelihoods$ub1<-likelihoods$logor1+(qnorm(0.975)*likelihoods$logse1)
likelihoods$lb2<-likelihoods$logor2-(qnorm(0.975)*likelihoods$logse2)
likelihoods$ub2<-likelihoods$logor2+(qnorm(0.975)*likelihoods$logse2)
likelihoods$lb3<-likelihoods$logor3-(qnorm(0.975)*likelihoods$logse3)
likelihoods$ub3<-likelihoods$logor3+(qnorm(0.975)*likelihoods$logse3)
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
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 5 (red),5-13 (green) and over 13 (blue)") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12),
axis.text.y=element_text(size=12))
two<-ggplot(data=likelihoods, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/56)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
geom_vline(aes(xintercept=2.2), colour="red") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12),
axis.text.y=element_text(size=12)) +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <5 and >13"))))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/inds_het_tests_all_inc_midrange_lessthan_5.png", res=800,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()

#show all loci with p<0.5 just so it isn't such a busy Figure:
r$pfdr<-p.adjust(r$p, method = "BH")
r1<-r[r$p<0.1,]
r1$ord<-c(nrow(r1):1)
r1$loci<-reorder(r1$loci,r1$ord)
o<-r1[r1$pfdr<0.1,]
o<-o[o$pfdr==max(o$pfdr),]

one<-ggplot(data=r1, aes(x=logor3,y=as.factor(loci))) + geom_point(,colour="blue") +
geom_point(data=r1, aes(x=logor2,as.numeric(loci)+0.2), colour="green") +
geom_point(data=r1, aes(x=logor1,as.numeric(loci)+0.4), colour="red") +
geom_errorbarh(data=r1, aes(xmin=lb3,xmax=ub3, y=as.factor(loci)),colour="blue", height=0.01) +
geom_errorbarh(data=r1, aes(xmin=lb2,xmax=ub2, y=as.numeric(loci)+0.2),colour="green", height=0.01) +
geom_errorbarh(data=r1, aes(xmin=lb1,xmax=ub1, y=as.numeric(loci)+0.4),colour="red", height=0.01) +
geom_vline(xintercept=0, colour="red", linetype="dashed") +
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 5 (red),5-13 (green) and over 13 (blue)") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12),
axis.text.y=element_text(size=12))
two<-ggplot(data=r1, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/55)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
geom_vline(aes(xintercept=2.2), colour="red") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <5 and >13")))) +
coord_cartesian(xlim=c(0,5)) +
theme(axis.title.y=element_text(size=12),
axis.text.y=element_text(size=12))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/het_tests_midrange_just0.1_lessthan_5.png", res=800,
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
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 5 (red),5-13 (green) and over 13 (blue)") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=15),
axis.text.y=element_text(size=12))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/het_tests_all_inc_midrange_no_nonsig_lessthan_5.png", res=800,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()

