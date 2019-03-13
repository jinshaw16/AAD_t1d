#multinomial_readin.R
library(ggplot2)
library(gridExtra)
library(snpStats)
library(annotSnpStats)

load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_n.R")
t1dsnps$altid<-ifelse(substr(t1dsnps$id,1,1)=="1",paste0("X",t1dsnps$id),t1dsnps$id)
colnames(pheno)<-ifelse(substr(colnames(pheno),1,1)=="1",paste0("X",colnames(pheno)),colnames(pheno))
t1dsnps$altid<-ifelse(!t1dsnps$altid %in% colnames(pheno),t1dsnps$snp,t1dsnps$altid)

getlik<-function(snpname){
load(file=paste0("/well/todd/users/jinshaw/aad/under_7/results/",snpname,"_n.RData"))
return(l)
}
likelihoods<-lapply(t1dsnps$altid,getlik)
likelihoods<-do.call("rbind", likelihoods)
likelihoods$snp=t1dsnps$snp
likelihoods$id<-t1dsnps$id
likelihoods$loci<-t1dsnps$loci
likelihoods$altid<-t1dsnps$altid
write.table(likelihoods, file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_n.txt",
col.names=T, row.names=F, quote=F, sep="\t")

likelihoods<-read.table(file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_n.txt",
header=T, sep="\t", as.is=T)
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
geom_vline(aes(xintercept=log10(0.05/nrow(r))*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
geom_vline(aes(xintercept=fdrline), colour="red") +
scale_y_discrete(name="Locus") +
theme(axis.title.y=element_text(size=12),
axis.text.y=element_text(size=12)) +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <7 and >13"))))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/inds_het_tests_all_inc_midrange_n.png", res=800,
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
geom_vline(aes(xintercept=fdrline), colour="red", linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <7 and >13")))) +
coord_cartesian(xlim=c(0,5)) +
theme(axis.title.y=element_text(size=15),
axis.text.y=element_text(size=12))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/het_tests_all_inc_midrange_no_nonsig_n.png", res=800,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()

#of those not FDR associateed, how many have strongest signal in the <7s, middle in the 7-13s and weakest in >13s?
not<-r[r$pfdr>0.1,]
not$hit<-ifelse((not$logor1<not$logor2 & not$logor2 < not$logor3 & abs(not$logor1)>abs(not$logor3)) |
(not$logor1 > not$logor2 & not$logor2 > not$logor3 & abs(not$logor1)>abs(not$logor3)),1,0)
t<-table(not$hit)
pbinom(q=t[2],size=(t[1]+t[2]), p=1/6,lower.tail=F)

########################
#NOW FOR THE MAFS PLOTS#
########################
#read in all the unrelateds - including the mid range individuals
load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_n.R")
r<-read.table(file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_n.txt",
header=T, as.is=T, sep="\t")
r$pfdr<-p.adjust(r$p, method = "BH")
r1<-r[r$pfdr<0.1,]

#get most likely call for the variants:
g<-as(b,"SnpMatrix")
a<-as(g,"numeric")
samp<-pheno[,!colnames(pheno) %in% t1dsnps$id]
samp<-cbind(samp,a)
cases<-samp[samp$affected==2,]
cases$logage<-log(cases$onset)



#plots by allele:
library(plyr)
library(dplyr)
hits<-r1$altid

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
ggsave(d,file=paste0("/well/todd/users/jinshaw/output/aad/under_7/mafs/redo/",gsub("/"," ",t1dsnps[t1dsnps$id==snp,"loci"]),"_by_allele_n.png"),
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
samp<-samp[rownames(g),]

cs<-col.summary(g)
w<-which(cs$RAF>0.5)
g<-switch.alleles(g, snps=w)

samp$onset<-ifelse(samp$onset>=16 & samp$onset<=20,18,
ifelse(samp$onset>20 & samp$onset<=25, 23,
ifelse(samp$onset>25 & samp$onset<=30, 28,
ifelse(samp$onset>30 & samp$onset<=35, 33,
ifelse(samp$onset>35, 37,samp$onset)))))
t<-names(table(samp$onset))

getmafsall<-function(snp){
getcs<-function(onset){
if(onset==0){
a<-g[rownames(g) %in% rownames(samp[samp$affected==1,]),snp]
}
if(onset!=0){
a<-g[rownames(g) %in% rownames(samp[samp$affected==2 & samp$onset==onset,]),snp]
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

ggsave(g,file=paste0("/well/todd/users/jinshaw/output/aad/under_7/mafs/redo/",gsub("/"," ",t1dsnps[t1dsnps$id==snp,"loci"]),"_all_n.png"),
dpi=400, height=20, width=20, units="cm")
}

lapply(hits,getmafsall)


