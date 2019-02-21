#hla_imputation_readin_prot_redo_3.R
#redo_3 changes the HLA class I associations to be adjusted for the other class I alleles.
#as well as the class II alleles adjusted for class I
#and finally codes the HLA alleles I think appropriately with regards to DR3 and DR4 independent effects (excluding the DR34s altogether).
#also adjusting here for 10 PCs

library(ggplot2)
library(epicalc)
library(gridExtra)
library(multinomRob)

ord<-read.table(file="/well/todd/users/jinshaw/aad/under_7/all_ages.fam", header=F,as.is=F)

#read imputation results into R:
load(file="/well/todd/users/jinshaw/aad/under_7/hla_all.RData")
h<-hla

hla<-hla[hla$group %in% c(0,1,3),]

library(multinomRob)
hla$g0<-ifelse(hla$group==0,1,0)
hla$g1<-ifelse(hla$group==1,1,0)
hla$g2<-ifelse(hla$group==3,1,0)

hlano34<-hla[hla$dr34_1!=1,]
getlikelihoods<-function(hap, frame, adjusted){
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ",adjusted,"`",hap,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ",adjusted,"`",hap,"`"))
if(hap=="dr11_1"){
l<-data.frame(unconstrained=NA,constrained=NA,logor1=NA, logse1=NA, logor2=NA, logse2=NA,loglambda=NA,
chisq=NA,p=NA,lb1=NA,ub1=NA,lb2=NA,ub2=NA,logp=NA)
}
if(hap!="dr11_1"){
one<-multinomRob(model=list(form0,form1,form2),data=frame, print.level=1, MLEonly=T)
llk1<-one$value

equals0<-as.formula(paste0("g1 ~ `",hap,"` + 0"))
equals1<-as.formula(paste0("g2 ~ `",hap,"` + 0"))
#model 2: constrain betas to equal each other:
two<-multinomRob(model=list(form0,form1,form2),
equality=list(list(equals0,equals1)), data=frame, print.level=1, MLEonly=T)
llk2<-two$value
l<-data.frame(unconstrained=llk1, constrained=llk2, logor1=one$coefficients[nrow(one$coefficients),2], logse1=one$se[nrow(one$coefficients),2],
logor2=one$coefficients[nrow(one$coefficients),3], logse2=one$se[nrow(one$coefficients),3])
l$constrained<-l$constrained*-1
l$unconstrained<-l$unconstrained*-1
l$loglambda<-l$constrained-l$unconstrained
l$chisq<-l$loglambda*-2
l$p<-pchisq(l$chisq,1,lower.tail=F)
l$lb1<-l$logor1-(qnorm(0.975)*l$logse1)
l$ub1<-l$logor1+(qnorm(0.975)*l$logse1)
l$lb2<-l$logor2-(qnorm(0.975)*l$logse2)
l$ub2<-l$logor2+(qnorm(0.975)*l$logse2)
l$logp<-log10(l$p)*-1
}
return(l)
}

l_1<-mapply(getlikelihoods,hap=c("dr3_1","dr4_1","dr34_1","dr11_1","dr13_1","dr7_1","dr14_1","dr15_1"), 
frame=list(hlano34,hlano34,hla,hlano34,hlano34,hlano34,hlano34,hlano34),
adjusted=c("a0201_1 + a2402_1 + b1801_1 + b3906_1 + b4403_1 +",
"a0201_1 + a2402_1 + b1801_1 + b3906_1 + b4403_1 +",
"a0201_1 + a2402_1 + b1801_1 + b3906_1 + b4403_1 +",
"","","","",""),SIMPLIFY=FALSE)
l_1<-do.call("rbind",l_1)
colnames(l_1)<-paste0(colnames(l_1),"_1")


#get the same thing but plotting the midrange betas too:
load(file="/well/todd/users/jinshaw/aad/under_7/hla_all.RData")
hla<-h
hla$g0<-ifelse(hla$group==0,1,0)
hla$g1<-ifelse(hla$group==1,1,0)
hla$g2<-ifelse(hla$group==2,1,0)
hla$g3<-ifelse(hla$group==3,1,0)
hlano34<-hla[hla$dr34_1!=1,]
getlikelihoodall<-function(hap, frame, adjusted){
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ",adjusted," `",hap,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ",adjusted,"`",hap,"`"))
form3<-as.formula(paste0("g3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ",adjusted,"`",hap,"`"))

if(hap %in% c("dr11","dr11_1","dr14_1")){
l<-data.frame(logor1=NA, logse1=NA,
logor2=NA, logse2=NA, logor3=NA,
logse3=NA, lb1=NA,ub1=NA,lb2=NA,
ub2=NA,lb3=NA,ub3=NA)
}
if(!hap %in% c("dr11","dr11_1","dr14_1")){
one<-multinomRob(model=list(form0,form1,form2,form3),data=frame, print.level=1, MLEonly=T)
llk1<-one$value

l<-data.frame(logor1=one$coefficients[nrow(one$coefficients),2], logse1=one$se[nrow(one$coefficients),2],
logor2=one$coefficients[nrow(one$coefficients),3], logse2=one$se[nrow(one$coefficients),3], logor3=one$coefficients[nrow(one$coefficients),4],
logse3=one$se[nrow(one$coefficients),4])
l$lb1<-l$logor1-(qnorm(0.975)*l$logse1)
l$ub1<-l$logor1+(qnorm(0.975)*l$logse1)
l$lb2<-l$logor2-(qnorm(0.975)*l$logse2)
l$ub2<-l$logor2+(qnorm(0.975)*l$logse2)
l$lb3<-l$logor3-(qnorm(0.975)*l$logse3)
l$ub3<-l$logor3+(qnorm(0.975)*l$logse3)
}
return(l)
}

l1_1<-mapply(getlikelihoodall,hap=c("dr3_1","dr4_1","dr34_1","dr11_1","dr13_1","dr7_1","dr14_1","dr15_1"),
frame=list(hlano34, hlano34, hla, hlano34, hlano34, hlano34, hlano34, hlano34),
adjusted=c("a0201_1 + a2402_1 + b1801_1 + b3906_1 + b4403_1 +",
"a0201_1 + a2402_1 + b1801_1 + b3906_1 + b4403_1 +",
"a0201_1 + a2402_1 + b1801_1 + b3906_1 + b4403_1 +",
"","","","",""), SIMPLIFY=FALSE)
l1_1<-do.call("rbind",l1_1)
colnames(l1_1)<-paste0(colnames(l1_1),"_1")


######################
#CLASS I association:#
######################
h<-hla
hla<-hla[hla$group %in% c(0,1,3),]
hla$group<-ifelse(hla$group==0,0,
ifelse(hla$group==1,1,
ifelse(hla$group==3,2,NA)))

hla$g0<-ifelse(hla$group==0,1,0)
hla$g1<-ifelse(hla$group==1,1,0)
hla$g2<-ifelse(hla$group==2,1,0)

getlikelihoodsadj<-function(hap,adj){
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + dr3_1 + dr4_1 + dr34_1 + ",adj," +`",hap,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + dr3_1 + dr4_1 + dr34_1 + ",adj," +`",hap,"`"))

if(hap %in% c("a0224_1")){
l<-data.frame(unconstrained=NA, constrained=NA, logor1=NA, logse1=NA,
logor2=NA, logse2=NA,loglambda=NA,chisq=NA,p=NA,logp=NA)
}
if(!hap %in% c("a0224_1")){
one<-multinomRob(model=list(form0,form1,form2),data=hla, print.level=1, MLEonly=T)
llk1<-one$value

equals0<-as.formula(paste0("g1 ~ `",hap,"` + 0"))
equals1<-as.formula(paste0("g2 ~ `",hap,"` + 0"))
#model 2: constrain betas to equal each other:
two<-multinomRob(model=list(form0,form1,form2),
equality=list(list(equals0,equals1)), data=hla, print.level=1, MLEonly=T)
llk2<-two$value
l<-data.frame(unconstrained=llk1, constrained=llk2, logor1=one$coefficients[nrow(one$coefficients),2], logse1=one$se[nrow(one$coefficients),2],
logor2=one$coefficients[nrow(one$se),3], logse2=one$se[nrow(one$se),3])
l$constrained=l$constrained*-1
l$unconstrained=l$unconstrained*-1
l$loglambda=l$constrained-l$unconstrained
l$chisq=l$loglambda*-2
l$p<-pchisq(l$chisq,1, lower.tail=F)
l$logp<-log10(l$p)*-1
}
return(l)
}
hla$dpb10301_1<-ifelse(hla$dpb10301_1==2,0,ifelse(hla$dpb10301_1==0,2,hla$dpb10301_1))
l2_1<-mapply(getlikelihoodsadj,hap=c("dpb10301_1","dpb10402_1","a0201_1","a2402_1","a3201_1","a1101_1","b1801_1", "b3906_1","b4403_1"),
adj=c("a0201_1 + a2402_1 + b1801_1 + b3906_1 + b4403_1",
"a0201_1 + a2402_1 +b1801_1 + b3906_1 + b4403_1",
"b1801_1 + b3906_1 + b4403_1",
"b1801_1 + b3906_1 + b4403_1",
"b1801_1 + b3906_1 + b4403_1",
"b1801_1 + b3906_1 + b4403_1",
"a0201_1 + a2402_1 + a3201_1",
"a0201_1 + a2402_1 + a3201_1",
"a0201_1 + a2402_1 + a3201_1"),SIMPLIFY=F)
l2_1<-do.call("rbind",l2_1)


#now do the 3 odds ratio comparison:
hla<-h
hla$g0<-ifelse(hla$group==0,1,0)
hla$g1<-ifelse(hla$group==1,1,0)
hla$g2<-ifelse(hla$group==2,1,0)
hla$g3<-ifelse(hla$group==3,1,0)

getlikelihoodsadj1<-function(hap,adj){
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + dr3_1 + dr4_1 + dr34_1 + ",adj," +`",hap,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + dr3_1 + dr4_1 + dr34_1 + ",adj," +`",hap,"`"))
form3<-as.formula(paste0("g3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + dr3_1 + dr4_1 + dr34_1 + ",adj," +`",hap,"`"))
if(hap %in% c("a0224_1")){
l<-data.frame(logor1=NA, logse1=NA,
logor2=NA, logse2=NA,
logor3=NA, logse3=NA,
lb1=NA,ub1=NA,lb2=NA,ub2=NA,lb3=NA,ub3=NA)
}
if(!hap %in% c("a0224_1")){
one<-multinomRob(model=list(form0,form1,form2, form3),data=hla, print.level=1, MLEonly=T)
l<-data.frame(logor1=one$coefficients[nrow(one$coefficients),2], logse1=one$se[nrow(one$coefficients),2],
logor2=one$coefficients[nrow(one$coefficients),3], logse2=one$se[nrow(one$coefficients),3],
logor3=one$coefficients[nrow(one$coefficients),4], logse3=one$se[nrow(one$coefficients),4])
if(hap %in% c("dpb10301_1")){
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
}
return(l)
}
hla$dpb10301_1<-ifelse(hla$dpb10301_1==2,0,ifelse(hla$dpb10301_1==0,2,hla$dpb10301_1))

lalls_1<-mapply(getlikelihoodsadj1,
hap=c("dpb10301_1","dpb10402_1","a0201_1","a2402_1","a3201_1","a1101_1","b1801_1", "b3906_1","b4403_1"),
adj=c("a0201_1 + a2402_1 + b1801_1 + b3906_1 + b4403_1",
"a0201_1 + a2402_1 +b1801_1 + b3906_1 + b4403_1",
"b1801_1 + b3906_1 + b4403_1",
"b1801_1 + b3906_1 + b4403_1",
"b1801_1 + b3906_1 + b4403_1",
"b1801_1 + b3906_1 + b4403_1",
"a0201_1 + a2402_1 + a3201_1",
"a0201_1 + a2402_1 + a3201_1",
"a0201_1 + a2402_1 + a3201_1"),SIMPLIFY=F)
lalls_1<-do.call("rbind",lalls_1)


#now we want to combine the class II results with the class I and plot all on one graph - likely will use for the publication.
l1_1$id<-c("DR3-DQ2","DR4-DQ8","DR3-DQ2/DR4-DQ8","DRB1*11:04-DQB1*03:01",
"DRB1*13:03-DQB1*03:01","DRB1*07:01-DQB1*03:03","DRB1*14:01-DQB1*05:03",
"DRB1*15:01-DQB1*06:02")
l1_1$logp<-l_1$logp
l1_1<-l1_1[order(-l1_1$logp),]
l1_1$ord<-c(nrow(l1_1):1)

lalls_1$id<-c("DPB1*03:01","DPB1*04:02",
"A*02:01","A*24:02","A*32:01", "A*11:01",
"B*18:01", "B*39:06","B44:03")
lalls_1$logp<-l2_1$logp

lout<-lalls_1[order(-lalls_1$logp),]
lout$ord<-c(nrow(lout):1)
lout$id<-as.factor(lout$id)
lout$id<-reorder(lout$id,lout$ord)
colnames(lout)[1:12]<-paste0(colnames(lout)[1:12],"_1")
#add in the class II:
lout<-rbind(lout,l1_1)
lout<-lout[!is.na(lout$logp) & !is.na(lout$logor1_1),]
lout<-lout[order(-lout$logp),]
lout$ord<-c(nrow(lout):1)
lout$id<-as.factor(lout$id)
lout$id<-reorder(lout$id,lout$ord)
save(lout,file="/well/todd/users/jinshaw/aad/under_7/results/hla_multinomial_3.RData")

#show all those with p<0.5:
load(file="/well/todd/users/jinshaw/aad/under_7/results/hla_multinomial_3.RData")
lout$p<-10^(-1*lout$logp)
lout1<-lout[lout$logp>0.301,]
lout1$ord<-c(nrow(lout1):1)
lout1$id<-as.factor(lout1$id)
lout1$id<-reorder(lout1$id,lout1$ord)


one_out<-ggplot(data=lout1, aes(x=logor3_1,y=as.factor(id))) + geom_point(,colour="blue") +
geom_point(data=lout1, aes(x=logor2_1,as.numeric(id)+0.2), colour="green") +
geom_point(data=lout1, aes(x=logor1_1,as.numeric(id)+0.4), colour="red") +
geom_errorbarh(data=lout1, aes(xmin=lb3_1,xmax=ub3_1, y=as.factor(id)),colour="blue", height=0.01) +
geom_errorbarh(data=lout1, aes(xmin=lb2_1,xmax=ub2_1, y=as.numeric(id)+0.2),colour="green", height=0.01) +
geom_errorbarh(data=lout1, aes(xmin=lb1_1,xmax=ub1_1, y=as.numeric(id)+0.4),colour="red", height=0.01) +
geom_vline(xintercept=0, colour="red", linetype="dashed") +
scale_x_continuous(name="T1D log-OR diagnosed < 7 (red), 7-13 (green) and > 13 (blue)") +
scale_y_discrete(name="Locus") +
coord_cartesian(xlim=c(-4,4)) +
theme(axis.title.x=element_text(size=8))
two_out<-ggplot(data=lout1, aes(logp, as.factor(id))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/nrow(lout))*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) heterogeneity between <7 and >13"))))+
theme(axis.title.x=element_text(size=8))

png(file="/well/todd/users/jinshaw/output/aad/under_7/hla/redo/multinom_p_0.5_hla_adjusted_dr34_3.png", res=400,
width=25, height=20, units="cm")
grid.arrange(one_out,two_out,ncol=2)
dev.off()



#now all tested:
lout<-lout[!is.na(lout$logor1_1),]
lout$ord<-c(nrow(lout):1)
lout$id<-as.factor(lout$id)
lout$id<-reorder(lout$id,lout$ord)

one_out<-ggplot(data=lout, aes(x=logor3_1,y=as.factor(id))) + geom_point(,colour="blue") +
geom_point(data=lout, aes(x=logor2_1,as.numeric(id)+0.2), colour="green") +
geom_point(data=lout, aes(x=logor1_1,as.numeric(id)+0.4), colour="red") +
geom_errorbarh(data=lout, aes(xmin=lb3_1,xmax=ub3_1, y=as.factor(id)),colour="blue", height=0.01) +
geom_errorbarh(data=lout, aes(xmin=lb2_1,xmax=ub2_1, y=as.numeric(id)+0.2),colour="green", height=0.01) +
geom_errorbarh(data=lout, aes(xmin=lb1_1,xmax=ub1_1, y=as.numeric(id)+0.4),colour="red", height=0.01) +
geom_vline(xintercept=0, colour="red", linetype="dashed") +
scale_x_continuous(name="T1D log-OR diagnosed < 7 (red), 7-13 (green) and > 13 (blue)") +
scale_y_discrete(name="Locus") +
coord_cartesian(xlim=c(-4,4)) +
theme(axis.title.x=element_text(size=8))
two_out<-ggplot(data=lout, aes(logp, as.factor(id))) + geom_point() +
geom_vline(aes(xintercept=log10(0.05/nrow(lout))*-1),colour="red",linetype="dashed") +
geom_vline(aes(xintercept=log10(0.05)*-1),colour="red",linetype="dotted") +
scale_y_discrete(name="Locus") +
scale_x_continuous(name=bquote("-log"[10]~.(paste0("(p) heterogeneity between <7 and >13"))))+
theme(axis.title.x=element_text(size=8))

png(file="/well/todd/users/jinshaw/output/aad/under_7/hla/redo/multinom_all_hla_adjusted_dr34_3.png", res=800,
width=25, height=20, units="cm")
grid.arrange(one_out,two_out,ncol=2)
dev.off()


