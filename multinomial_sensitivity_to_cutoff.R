#multinomial_sensitivity_to_cutoff.R
library(snpStats)
library(annotSnpStats)
library(snpStatsWriter)
library(humarray)
library(gridExtra)
library(multinomRob)
library(ggplot2)

#load data:
domultiall<-function(lcutoff){
load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_2.R")

t1dsnps$altid<-ifelse(substr(t1dsnps$id,1,1)=="1",paste0("X",t1dsnps$id),t1dsnps$id)
colnames(pheno)<-ifelse(substr(colnames(pheno),1,1)=="1",paste0("X",colnames(pheno)),colnames(pheno))
t1dsnps$altid<-ifelse(!t1dsnps$altid %in% colnames(pheno),t1dsnps$snp,t1dsnps$altid)

#now perform multinomial regression of <lcutoff compared to >13:
pheno$group<-ifelse(pheno$onset<lcutoff & !is.na(pheno$onset) & pheno$affected==2,1,
ifelse(pheno$onset>=lcutoff & pheno$onset<13 & !is.na(pheno$onset) & pheno$affected==2,2,
ifelse(pheno$onset>=13 & !is.na(pheno$onset)    & pheno$affected==2,3,
ifelse(pheno$affected==1,0,NA))))
pheno$g0<-ifelse(pheno$group==0,1,ifelse(pheno$group!=0,0,NA))
pheno$g1<-ifelse(pheno$group==1,1,ifelse(pheno$group!=1,0,NA))
pheno$g2<-ifelse(pheno$group==2,1,ifelse(pheno$group!=2,0,NA))
pheno$g3<-ifelse(pheno$group==3,1,ifelse(pheno$group!=3,0,NA))
pheno$`seq-rs35667974`<-ifelse(pheno$`seq-rs35667974`==2,0,
ifelse(pheno$`seq-rs35667974`==0,2,pheno$`seq-rs35667974`))

#carry out the multinomial regressions:
getlikelihoods<-function(snpname){
#model 1: allow the beta for the SNP to vary for both:
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form3<-as.formula(paste0("g3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))

one<-multinomRob(model=list(form0,form1,form2,form3),data=pheno, print.level=1, MLEonly=T)
llk1<-one$value
equals0<-as.formula(paste0("g1 ~ `",snpname,"` + 0"))
equals1<-as.formula(paste0("g3 ~ `",snpname,"` + 0"))
#model 2: constrain betas to equal each other:
two<-multinomRob(model=list(form0,form1,form2,form3),
equality=list(list(equals0,equals1)), data=pheno, print.level=1, MLEonly=T)
llk2<-two$value
l<-data.frame(unconstrained=llk1, constrained=llk2, logor1=one$coefficients[nrow(one$coefficients),2], logse1=one$se[nrow(one$coefficients),2],
logor2=one$coefficients[nrow(one$coefficients),3], logse2=one$se[nrow(one$coefficients),3],
logor3=one$coefficients[nrow(one$coefficients),4],logse3=one$se[nrow(one$coefficients),4])
if (snpname %in% c("seq-rs35667974")){
l$logor1<-l$logor1*-1
l$logor2<-l$logor2*-1
l$logor3<-l$logor3*-1
}
l$lb1<-l$logor1-(qnorm(0.975)*l$logse1)
l$lb2<-l$logor2-(qnorm(0.975)*l$logse2)
l$lb3<-l$logor3-(qnorm(0.975)*l$logse3)
l$ub1<-l$logor1+(qnorm(0.975)*l$logse1)
l$ub2<-l$logor2+(qnorm(0.975)*l$logse2)
l$ub3<-l$logor3+(qnorm(0.975)*l$logse3)
l$constrained=l$constrained*-1
l$unconstrained=l$unconstrained*-1
l$loglambda<-l$constrained-l$unconstrained
l$chisq<-l$loglambda*-2
l$p<-pchisq(l$chisq,1, lower.tail=F)
return(l)
}

likelihoods<-lapply(t1dsnps$altid,getlikelihoods)
likelihoods<-do.call("rbind", likelihoods)
likelihoods$snp=t1dsnps$snp
likelihoods$id<-t1dsnps$id
likelihoods$loci<-t1dsnps$loci
write.table(likelihoods, file=paste0("/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_lessthan_",lcutoff,".txt"),
col.names=T, row.names=F, quote=F, sep="\t")

likelihoods<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_lessthan_",lcutoff,".txt"),
header=T, as.is=T, sep="\t")
likelihoods$logp<-log10(likelihoods$p)*-1
likelihoods<-likelihoods[order(-likelihoods$logp),]
likelihoods$loci<-as.factor(likelihoods$loci)
likelihoods$ord<-c(nrow(likelihoods):1)
likelihoods$loci<-reorder(likelihoods$loci,likelihoods$ord)
r<-likelihoods
r$pfdr<-p.adjust(r$p, method = "BH")
r1<-r[r$pfdr<0.1,]
r2<-r[r$pfdr>=0.1,]
fdrline<-(max(r2$logp)+(min(r1$logp)-max(r2$logp))/2)
one<-ggplot(data=r, aes(x=logor3,y=as.factor(loci))) + geom_point(,colour="blue") +
geom_point(data=r, aes(x=logor2,as.numeric(loci)+0.2), colour="green") +
geom_point(data=r, aes(x=logor1,as.numeric(loci)+0.4), colour="red") +
geom_errorbarh(data=r, aes(xmin=lb3,xmax=ub3, y=as.factor(loci)),colour="blue", height=0.01) +
geom_errorbarh(data=r, aes(xmin=lb2,xmax=ub2, y=as.numeric(loci)+0.2),colour="green", height=0.01) +
geom_errorbarh(data=r, aes(xmin=lb1,xmax=ub1, y=as.numeric(loci)+0.4),colour="red", height=0.01) +
geom_vline(xintercept=0, colour="red", linetype="dashed") +
scale_x_continuous(name=paste0("T1D log-odds ratio for those diagnosed under ",lcutoff," (red),",lcutoff,"-13 (green) and over 13 (blue)")) +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12),
axis.text.y=element_text(size=12))
two<-ggplot(data=likelihoods, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/nrow(r))*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
geom_vline(aes(xintercept=fdrline), colour="red") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12), 
axis.text.y=element_text(size=12)) +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <",lcutoff," and >13"))))

png(file=paste0("/well/todd/users/jinshaw/output/aad/under_7/multinom/redo_1/inds_het_tests_all_inc_midrange_lessthan_",lcutoff,".png"), res=800,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()

one<-ggplot(data=r, aes(x=logor3,y=as.factor(loci))) + geom_point(,colour="blue") +
geom_point(data=r, aes(x=logor2,as.numeric(loci)+0.2), colour="green") +
geom_point(data=r, aes(x=logor1,as.numeric(loci)+0.4), colour="red") +
geom_errorbarh(data=r, aes(xmin=lb3,xmax=ub3, y=as.factor(loci)),colour="blue", height=0.01) +
geom_errorbarh(data=r, aes(xmin=lb2,xmax=ub2, y=as.numeric(loci)+0.2),colour="green", height=0.01) +
geom_errorbarh(data=r, aes(xmin=lb1,xmax=ub1, y=as.numeric(loci)+0.4),colour="red", height=0.01) +
geom_vline(xintercept=0, colour="red", linetype="dashed") +
scale_x_continuous(name=paste0("T1D log-odds ratio for those diagnosed under ",lcutoff," (red), ",lcutoff,"-13 (green) and over 13 (blue)")) +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12),
axis.text.y=element_text(size=12),
axis.text.x=element_text(size=9),
axis.title.x=element_text(size=9, hjust=0.94))
two<-ggplot(data=likelihoods, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/nrow(r))*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
geom_vline(aes(xintercept=fdrline), colour="red") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12),
axis.text.y=element_text(size=12),
axis.text.x=element_text(size=9),
axis.title.x=element_text(size=9, hjust=1.5)) +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <",lcutoff," and >13"))))


png(file=paste0("/well/todd/users/jinshaw/output/aad/under_7/multinom/redo_1/inds_het_tests_all_inc_midrange_lessthan_",lcutoff,"_sm.png"), res=800,
width=25, height=20, units="cm")
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
scale_x_continuous(name=paste0("T1D log-odds ratio for those diagnosed under ",lcutoff," (red),",lcutoff,"-13 (green) and over 13 (blue)")) +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=15), 
axis.text.y=element_text(size=12))

two<-ggplot(data=r1, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/nrow(r))*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=fdrline), colour="red", linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <",lcutoff," and >13")))) +
coord_cartesian(xlim=c(0,5)) +
theme(axis.title.y=element_text(size=15), 
axis.text.y=element_text(size=12))

png(file=paste0("/well/todd/users/jinshaw/output/aad/under_7/multinom/redo_1/het_tests_all_inc_midrange_no_nonsig_lessthan_",lcutoff,".png"), res=800,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()

one<-ggplot(data=r1, aes(x=logor3,y=as.factor(loci))) + geom_point(,colour="blue") +
geom_point(data=r1, aes(x=logor2,as.numeric(loci)+0.2), colour="green") +
geom_point(data=r1, aes(x=logor1,as.numeric(loci)+0.4), colour="red") +
geom_errorbarh(data=r1, aes(xmin=lb3,xmax=ub3, y=as.factor(loci)),colour="blue", height=0.01) +
geom_errorbarh(data=r1, aes(xmin=lb2,xmax=ub2, y=as.numeric(loci)+0.2),colour="green", height=0.01) +
geom_errorbarh(data=r1, aes(xmin=lb1,xmax=ub1, y=as.numeric(loci)+0.4),colour="red", height=0.01) +
geom_vline(xintercept=0, colour="red", linetype="dashed") +
scale_x_continuous(name=paste0("T1D log-odds ratio for those diagnosed under ",lcutoff," (red),",lcutoff,"-13 (green) and over 13 (blue)")) +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=15),
axis.text.y=element_text(size=12),
axis.text.x=element_text(size=9),
axis.title.x=element_text(size=9, hjust=0.94))

two<-ggplot(data=r1, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/nrow(r))*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=fdrline), colour="red", linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <",lcutoff," and >13")))) +
coord_cartesian(xlim=c(0,5)) +
theme(axis.title.y=element_text(size=15),
axis.text.y=element_text(size=12),
axis.text.x=element_text(size=9),
axis.title.x=element_text(size=9, hjust=1.5))


png(file=paste0("/well/todd/users/jinshaw/output/aad/under_7/multinom/redo_1/het_tests_all_inc_midrange_no_nonsig_lessthan_",lcutoff,"_sm.png"), res=800,
width=25, height=20, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()
}

domultiall(6)
domultiall(5)
