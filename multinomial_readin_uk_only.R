#multinomial_readin_uk_only.R


load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_uk_2.R")
t1dsnps$altid<-ifelse(substr(t1dsnps$id,1,1)=="1",paste0("X",t1dsnps$id),t1dsnps$id)
colnames(pheno)<-ifelse(substr(colnames(pheno),1,1)=="1",paste0("X",colnames(pheno)),colnames(pheno))
t1dsnps$altid<-ifelse(!t1dsnps$altid %in% colnames(pheno),t1dsnps$snp,t1dsnps$altid)

getlik<-function(snp){
load(file=paste0("/well/todd/users/jinshaw/aad/under_7/results/uk/",snp,".RData"))
return(l)
}
likelihoods<-lapply(t1dsnps$altid,getlik)
likelihoods<-do.call("rbind", likelihoods)
likelihoods$snp=t1dsnps$snp
likelihoods$id<-t1dsnps$id
likelihoods$loci<-t1dsnps$loci
write.table(likelihoods, file="/well/todd/users/jinshaw/aad/under_7/results/uk/inds_likelihoods_redo_n.txt",
col.names=T, row.names=F, quote=F, sep="\t")

likelihoods<-read.table(file="/well/todd/users/jinshaw/aad/under_7/results/uk/inds_likelihoods_redo_n.txt",
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
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 7 (red), 7-13 (green) and over 13 (blue)") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12),
axis.text.y=element_text(size=12))
two<-ggplot(data=likelihoods, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/56)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
geom_vline(aes(xintercept=fdrline), colour="red") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12),
axis.text.y=element_text(size=12)) +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <7 and >13"))))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo_1/uk/inds_het_tests_all_inc_midrange_n.png", res=800,
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
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 7 (red), 7-13 (green) and over 13 (blue)") +
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
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <7 and >13"))))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo_1/uk/inds_het_tests_all_inc_midrange_n_sm.png", res=800,
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
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 7 (red), 7-13 (green) and over 13 (blue)") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=15),
axis.text.y=element_text(size=12))

two<-ggplot(data=r1, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/55)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=fdrline), colour="red", linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <7 and >13")))) +
coord_cartesian(xlim=c(0,5)) +
theme(axis.title.y=element_text(size=15),
axis.text.y=element_text(size=12))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo_1/uk/het_tests_all_inc_midrange_no_nonsig_n.png", res=800,
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
scale_x_continuous(name="T1D log-odds ratio for those diagnosed under 7 (red), 7-13 (green) and over 13 (blue)") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=15),
axis.text.y=element_text(size=12),
axis.text.x=element_text(size=9),
axis.title.x=element_text(size=9, hjust=0.94))

two<-ggplot(data=r1, aes(logp, as.factor(loci))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/55)*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=fdrline), colour="red", linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <7 and >13")))) +
coord_cartesian(xlim=c(0,5)) +
theme(axis.title.y=element_text(size=15),
axis.text.y=element_text(size=12),
axis.text.x=element_text(size=9),
axis.title.x=element_text(size=9, hjust=1.5))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo_1/uk/het_tests_all_inc_midrange_no_nonsig_n_sm.png", res=800,
width=25, height=20, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()
