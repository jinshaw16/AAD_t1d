#multinomial_redo.R
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
ifelse(all@samples$onset>=13 & !is.na(all@samples$onset),2,NA)))


#Define the 55 T1D regions we want to examine (from various literature):
t1dsnps<-c("rs2476601","rs78037977",  "rs6691977", "rs3024505", "rs13415583",
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
p<-t1dsnps[!t1dsnps$id %in% colnames(all),]


#keep only the ones we want to test:
g<-all[,colnames(all) %in% t1dsnps$id]
g<-as(g,"SnpMatrix")
cs<-col.summary(g)
w<-which(cs$RAF>0.5)
g<-switch.alleles(g,snps=w)
b<-as(g,"numeric")

#and the imputed ones:
getsnp<-function(snp){
load(file=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/",snp,"_n_clean.RData"))
cs<-col.summary(DATA)
w<-which(cs$RAF>0.5)
DATA<-switch.alleles(DATA,snps=w)
d<-as(DATA[rownames(DATA) %in% rownames(g),grepl(snp,colnames(DATA))],"numeric")
colnames(d)<-snp
b<<-cbind2(b,d)
return(d)
}
l<-lapply(p$id,getsnp)
colnames(b)<-gsub(":.*","",colnames(b))
pheno<-all@samples

table(rownames(b)==rownames(pheno))
pheno<-cbind(pheno,b)
p1<-pheno
save(b,pheno,t1dsnps, file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_n.R")

load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_n.R")
t1dsnps$altid<-ifelse(substr(t1dsnps$id,1,1)=="1",paste0("X",t1dsnps$id),t1dsnps$id)
colnames(pheno)<-ifelse(substr(colnames(pheno),1,1)=="1",paste0("X",colnames(pheno)),colnames(pheno))
t1dsnps$altid<-ifelse(!t1dsnps$altid %in% colnames(pheno),t1dsnps$snp,t1dsnps$altid)

#now perform multinomial regression (remove 7-13s for the first 2 models):
ph<-pheno[pheno$group %in% c(0,1,2) & !is.na(pheno$group),]
ph$g0<-ifelse(ph$group==0,1,ifelse(ph$group!=0,0,NA))
ph$g1<-ifelse(ph$group==1,1,ifelse(ph$group!=1,0,NA))
ph$g2<-ifelse(ph$group==2,1,ifelse(ph$group!=2,0,NA))
pheno$group<-ifelse(pheno$affected==1,0,
ifelse(pheno$onset<7 & !is.na(pheno$onset),1,
ifelse(pheno$onset>=7 &pheno$onset<13 & !is.na(pheno$onset),2,
ifelse(pheno$onset>=13 & !is.na(pheno$onset),3,NA))))
pheno$g0<-ifelse(pheno$group==0,1,ifelse(pheno$group!=0,0,NA))
pheno$g1<-ifelse(pheno$group==1,1,ifelse(pheno$group!=1,0,NA))
pheno$g2<-ifelse(pheno$group==2,1,ifelse(pheno$group!=2,0,NA))
pheno$g3<-ifelse(pheno$group==3,1,ifelse(pheno$group!=3,0,NA))



#carry out the multinomial regressions:
getlikelihoods<-function(snpname){
#model 1: allow the beta for the SNP to vary for both:
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))

equals1<-as.formula(paste0("g1 ~ PC1 + 0"))
equals2<-as.formula(paste0("g2 ~ PC1 + 0"))
equals3<-as.formula(paste0("g1 ~ PC2 + 0"))
equals4<-as.formula(paste0("g2 ~ PC2 + 0"))
equals5<-as.formula(paste0("g1 ~ PC3 + 0"))
equals6<-as.formula(paste0("g2 ~ PC3 + 0"))
equals7<-as.formula(paste0("g1 ~ PC4 + 0"))
equals8<-as.formula(paste0("g2 ~ PC4 + 0"))
equals9<-as.formula(paste0("g1 ~ PC5 + 0"))
equals10<-as.formula(paste0("g2 ~ PC5 + 0"))
equals11<-as.formula(paste0("g1 ~ PC6 + 0"))
equals12<-as.formula(paste0("g2 ~ PC6 + 0"))
equals13<-as.formula(paste0("g1 ~ PC7 + 0"))
equals14<-as.formula(paste0("g2 ~ PC7 + 0"))
equals15<-as.formula(paste0("g1 ~ PC8 + 0"))
equals16<-as.formula(paste0("g2 ~ PC8 + 0"))
equals17<-as.formula(paste0("g1 ~ PC9 + 0"))
equals18<-as.formula(paste0("g2 ~ PC9 + 0"))
equals19<-as.formula(paste0("g1 ~ PC10 + 0"))
equals20<-as.formula(paste0("g2 ~ PC10 + 0"))
one<-multinomRob(model=list(form0,form1,form2),
equality=list(list(equals1,equals2),
list(equals3,equals4),
list(equals5,equals6),
list(equals7,equals8),
list(equals9,equals10),
list(equals11,equals12),
list(equals13,equals14),
list(equals15,equals16),
list(equals17,equals18),
list(equals19,equals20)),data=ph, print.level=1, MLEonly=T)
llk1<-one$value

equals21<-as.formula(paste0("g1 ~ `",snpname,"` + 0"))
equals22<-as.formula(paste0("g2 ~ `",snpname,"` + 0"))
#model 2: constrain betas to equal each other:
two<-multinomRob(model=list(form0,form1,form2),
equality=list(list(equals1,equals2),
list(equals3,equals4),
list(equals5,equals6),
list(equals7,equals8),
list(equals9,equals10),
list(equals11,equals12),
list(equals13,equals14),
list(equals15,equals16),
list(equals17,equals18),
list(equals19,equals20),
list(equals21,equals22),list(equals21,equals22)), data=ph, print.level=1, MLEonly=T)
llk2<-two$value
l<-data.frame(unconstrained=llk1, constrained=llk2)
l$constrained=l$constrained*-1
l$unconstrained=l$unconstrained*-1
l$loglambda<-l$constrained-l$unconstrained
#compare using likelihood ratio test:
l$chisq<-l$loglambda*-2
l$p<-pchisq(l$chisq,1, lower.tail=F)
#model 3: get estimates for all 3 age group strata:
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
form3<-as.formula(paste0("g3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + `",snpname,"`"))
equals1<-as.formula(paste0("g1 ~ PC1 + 0"))
equals2<-as.formula(paste0("g2 ~ PC1 + 0"))
equals3<-as.formula(paste0("g3 ~ PC1 + 0"))
equals4<-as.formula(paste0("g1 ~ PC2 + 0"))
equals5<-as.formula(paste0("g2 ~ PC2 + 0"))
equals6<-as.formula(paste0("g3 ~ PC3 + 0"))
equals7<-as.formula(paste0("g1 ~ PC3 + 0"))
equals8<-as.formula(paste0("g2 ~ PC3 + 0"))
equals9<-as.formula(paste0("g3 ~ PC3 + 0"))
equals10<-as.formula(paste0("g1 ~ PC4 + 0"))
equals11<-as.formula(paste0("g2 ~ PC4 + 0"))
equals12<-as.formula(paste0("g3 ~ PC4 + 0"))
equals13<-as.formula(paste0("g1 ~ PC5 + 0"))
equals14<-as.formula(paste0("g2 ~ PC5 + 0"))
equals15<-as.formula(paste0("g3 ~ PC5 + 0"))
equals16<-as.formula(paste0("g1 ~ PC6 + 0"))
equals17<-as.formula(paste0("g2 ~ PC6 + 0"))
equals18<-as.formula(paste0("g3 ~ PC6 + 0"))
equals19<-as.formula(paste0("g1 ~ PC7 + 0"))
equals20<-as.formula(paste0("g2 ~ PC7 + 0"))
equals21<-as.formula(paste0("g3 ~ PC7 + 0"))
equals22<-as.formula(paste0("g1 ~ PC8 + 0"))
equals23<-as.formula(paste0("g2 ~ PC8 + 0"))
equals24<-as.formula(paste0("g3 ~ PC8 + 0"))
equals25<-as.formula(paste0("g1 ~ PC9 + 0"))
equals26<-as.formula(paste0("g2 ~ PC9 + 0"))
equals27<-as.formula(paste0("g3 ~ PC9 + 0"))
equals28<-as.formula(paste0("g1 ~ PC10 + 0"))
equals29<-as.formula(paste0("g2 ~ PC10 + 0"))
equals30<-as.formula(paste0("g3 ~ PC10 + 0"))

three<-multinomRob(model=list(form0,form1,form2,form3),
equality=list(list(equals1,equals2,equals3),
list(equals4,equals5,equals6),
list(equals7,equals8,equals9),
list(equals10,equals11,equals12),
list(equals13,equals14,equals15),
list(equals16,equals17,equals18),
list(equals19,equals20,equals21),
list(equals22,equals23,equals24),
list(equals25,equals26,equals27),
list(equals28,equals29,equals30)),data=pheno, print.level=1, MLEonly=T)

l1<-data.frame(logor1=three$coefficients[nrow(three$coefficients),2], logse1=one$se[nrow(three$coefficients),2],
logor2=three$coefficients[nrow(three$coefficients),3], logse2=three$se[nrow(three$coefficients),3],
logor3=three$coefficients[nrow(three$coefficients),4],logse3=three$se[nrow(three$coefficients),4])
l1$lb1<-l1$logor1-(qnorm(0.975)*l1$logse1)
l1$ub1<-l1$logor1+(qnorm(0.975)*l1$logse1)
l1$lb2<-l1$logor2-(qnorm(0.975)*l1$logse2)
l1$ub2<-l1$logor2+(qnorm(0.975)*l1$logse2)
l1$lb3<-l1$logor3-(qnorm(0.975)*l1$logse3)
l1$ub3<-l1$logor3+(qnorm(0.975)*l1$logse3)
l<-cbind(l,l1)
return(l)
}

likelihoods<-lapply(t1dsnps$altid,getlikelihoods)
likelihoods<-do.call("rbind", likelihoods)
likelihoods$snp=t1dsnps$snp
likelihoods$id<-t1dsnps$id
likelihoods$loci<-t1dsnps$loci
likelihoods$altid<-t1dsnps$altid
write.table(likelihoods, file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_n.txt",
col.names=T, row.names=F, quote=F, sep="\t")

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
geom_vline(aes(xintercept=1.823909), colour="red", linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) Likelihood ratio test for heterogeneity between <7 and >13")))) +
coord_cartesian(xlim=c(0,5)) +
theme(axis.title.y=element_text(size=15), 
axis.text.y=element_text(size=12))

png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/het_tests_all_inc_midrange_no_nonsig_n.png", res=800,
width=40, height=30, units="cm")
grid.arrange(one,two,ncol=2)
dev.off()




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

