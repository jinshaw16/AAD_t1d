#hla_multinomial_n.R
#carries out multinomial logistic regression for the HLA classical alleles/haplotypes (class II)
#adjusts for top 10 PCs, codes haplotypes/classical alleles additively
#removes DR3/4s when examining other class II alleles as much confounding and want independent effects of e.g. DR3 and DR4 on their own.
#adjusts for as many other HLA haplotypes/alleles as possible for convergence in each analysis
#remove classical alleles/haplotypes from analysis if any group has less than 4 individuals (some of the rare class II haplotypes).

library(ggplot2)
library(epicalc)
library(gridExtra)
library(multinomRob)


#read imputation results into R:
load(file="/well/todd/users/jinshaw/aad/under_7/hla_all_n.RData")

library(multinomRob)
hla$g0<-ifelse(hla$group==0,1,0)
hla$g1<-ifelse(hla$group==1,1,0)
hla$g2<-ifelse(hla$group==2,1,0)
hla$g3<-ifelse(hla$group==3,1,0)

hlano34<-hla[hla$dr34_1!=1,]
getlikelihoods<-function(hap, frame, adjusted){
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ",adjusted,"`",hap,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ",adjusted,"`",hap,"`"))
form3<-as.formula(paste0("g3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + ",adjusted,"`",hap,"`"))

#having a rule that we need at least 5 from each group to even attempt fitting the model otherwise very very unstable (5 is unstable too, of course!)
t<-table(frame$group, frame[,hap])
s2<-sum(t[2,colnames(t) %in% c("1","2")])
s3<-sum(t[3,colnames(t) %in% c("1","2")])
s4<-sum(t[4,colnames(t) %in% c("1","2")])
if(s2<5 | s3<5 | s4<5){
l<-data.frame(unconstrained=NA,constrained=NA,logor1=NA, logse1=NA, logor2=NA, logse2=NA,logor3=NA, logse3=NA,loglambda=NA,
chisq=NA,p=NA,lb1=NA,ub1=NA,lb2=NA,ub2=NA,lb3=NA,ub3=NA,logp=NA)
}
if(s2>=5 & s3>=5 & s4>=5){
one<-multinomRob(model=list(form0,form1,form2,form3),data=frame, print.level=1, MLEonly=T)
llk1<-one$value

equals0<-as.formula(paste0("g1 ~ `",hap,"` + 0"))
equals1<-as.formula(paste0("g3 ~ `",hap,"` + 0"))
#model 2: constrain betas to equal each other:
two<-multinomRob(model=list(form0,form1,form2, form3),
equality=list(list(equals0,equals1)), data=frame, print.level=1, MLEonly=T)
llk2<-two$value
l<-data.frame(unconstrained=llk1, constrained=llk2, logor1=one$coefficients[nrow(one$coefficients),2], logse1=one$se[nrow(one$coefficients),2],
logor2=one$coefficients[nrow(one$coefficients),3], logse2=one$se[nrow(one$coefficients),3],
logor3=one$coefficients[nrow(one$coefficients),4], logse3=one$se[nrow(one$coefficients),4])
l$constrained<-l$constrained*-1
l$unconstrained<-l$unconstrained*-1
l$loglambda<-l$constrained-l$unconstrained
l$chisq<-l$loglambda*-2
l$p<-pchisq(l$chisq,1,lower.tail=F)
l$lb1<-l$logor1-(qnorm(0.975)*l$logse1)
l$ub1<-l$logor1+(qnorm(0.975)*l$logse1)
l$lb2<-l$logor2-(qnorm(0.975)*l$logse2)
l$ub2<-l$logor2+(qnorm(0.975)*l$logse2)
l$lb3<-l$logor3-(qnorm(0.975)*l$logse3)
l$ub3<-l$logor3+(qnorm(0.975)*l$logse3)
l$logp<-log10(l$p)*-1
}
return(l)
}

l_1<-mapply(getlikelihoods,hap=c("dr3_1","dr4_1","dr34_1","dr11_1","dr13_1","dr7_1","dr14_1","dr15_1"), 
frame=list(hlano34,hlano34,hla,hlano34,hlano34,hlano34,hlano34,hlano34),
adjusted=c("a0201_1 + a2402_1 + b1801_1 + b3906_1 + ",
"a0201_1 + a2402_1 + b1801_1 + b3906_1 + ",
"a0201_1 + a2402_1 + b1801_1 + b3906_1 + ",
"","","","",""),SIMPLIFY=FALSE)
l_1<-do.call("rbind",l_1)
l_1$id<-c("DR3-DQ2","DR4-DQ8","DR3-DQ2/DR4-DQ8","DRB1*11:04-DQB1*03:01",
"DRB1*13:03-DQB1*03:01","DRB1*07:01-DQB1*03:03","DRB1*14:01-DQB1*05:03",
"DRB1*15:01-DQB1*06:02")


######################
#CLASS I association:#
######################
getlikelihoodsadj<-function(hap,adj){
form0<-as.formula(paste0("g0 ~ 0"))
form1<-as.formula(paste0("g1 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + dr3_1 + dr4_1 + dr34_1 + ",adj," +`",hap,"`"))
form2<-as.formula(paste0("g2 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + dr3_1 + dr4_1 + dr34_1 + ",adj," +`",hap,"`"))
form3<-as.formula(paste0("g3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + dr3_1 + dr4_1 + dr34_1 + ",adj," +`",hap,"`"))

#having a rule that we need at least 5 from each group to even attempt fitting the model otherwise very very unstable (5 is unstable too, of course!)
t<-table(hla$group, hla[,hap])
s2<-sum(t[2,colnames(t) %in% c("1","2")])
s3<-sum(t[3,colnames(t) %in% c("1","2")])
s4<-sum(t[4,colnames(t) %in% c("1","2")])
if(s2<5 | s3<5 | s4<5){
l<-data.frame(unconstrained=NA, constrained=NA, logor1=NA, logse1=NA,
logor2=NA, logse2=NA,logor3=NA, logse3=NA,
lb1,ub1,lb2,ub2,lb3,ub3,loglambda=NA,chisq=NA,p=NA,logp=NA)
}
if(s2>=5 & s3>=5 & s4>=5){
one<-multinomRob(model=list(form0,form1,form2,form3),data=hla, print.level=1, MLEonly=T)
llk1<-one$value

equals0<-as.formula(paste0("g1 ~ `",hap,"` + 0"))
equals1<-as.formula(paste0("g3 ~ `",hap,"` + 0"))
#model 2: constrain betas to equal each other:
two<-multinomRob(model=list(form0,form1,form2, form3),
equality=list(list(equals0,equals1)), data=hla, print.level=1, MLEonly=T)
llk2<-two$value
l<-data.frame(unconstrained=llk1, constrained=llk2, logor1=one$coefficients[nrow(one$coefficients),2], logse1=one$se[nrow(one$coefficients),2],
logor2=one$coefficients[nrow(one$se),3], logse2=one$se[nrow(one$se),3],
logor3=one$coefficients[nrow(one$se),4], logse3=one$se[nrow(one$se),4])
l$constrained=l$constrained*-1
l$unconstrained=l$unconstrained*-1
l$lb1<-l$logor1*(qnorm(0.975)*l$logse1)
l$ub1<-l$logor1*(qnorm(0.975)*l$logse1)
l$lb2<-l$logor2*(qnorm(0.975)*l$logse2)
l$ub2<-l$logor2*(qnorm(0.975)*l$logse2)
l$lb3<-l$logor3*(qnorm(0.975)*l$logse3)
l$ub3<-l$logor3*(qnorm(0.975)*l$logse3)
l$loglambda=l$constrained-l$unconstrained
l$chisq=l$loglambda*-2
l$p<-pchisq(l$chisq,1, lower.tail=F)
l$logp<-log10(l$p)*-1
}
return(l)
}
l2_1<-mapply(getlikelihoodsadj,hap=c("dpb10301_1","dpb10402_1","a0201_1","a2402_1","a3201_1","a1101_1","b1801_1", "b3906_1","b4403_1"),
adj=c("a0201_1 + a2402_1 + b1801_1 + b3906_1",
"a0201_1 + a2402_1 +b1801_1 + b3906_1",
"b1801_1 + b3906_1",
"b1801_1 + b3906_1",
"b1801_1 + b3906_1",
"b1801_1 + b3906_1",
"a0201_1 + a2402_1 + a3201_1",
"a0201_1 + a2402_1 + a3201_1",
"a0201_1 + a2402_1 + a3201_1"),SIMPLIFY=F)
l2_1<-do.call("rbind",l2_1)
l2_1$id<-c("DPB1*03:01","DPB1*04:02",
"A*02:01","A*24:02","A*32:01", "A*11:01",
"B*18:01", "B*39:06", "B*44:03")


#now we want to combine the class II results with the class I and plot all on one graph - likely will use for the publication.

#add in the class II:
lout<-rbind(l_1,l2_1)
lout<-lout[!is.na(lout$logp) & !is.na(lout$logor1),]
lout<-lout[order(-lout$logp),]
lout$ord<-c(nrow(lout):1)
lout$id<-as.factor(lout$id)
lout$id<-reorder(lout$id,lout$ord)
save(lout,file="/well/todd/users/jinshaw/aad/under_7/results/hla_multinomial_n.RData")


load(file="/well/todd/users/jinshaw/aad/under_7/results/hla_multinomial_n.RData")

one_out<-ggplot(data=lout, aes(x=logor3,y=as.factor(id))) + geom_point(,colour="blue") +
geom_point(data=lout, aes(x=logor2,as.numeric(id)+0.2), colour="green") +
geom_point(data=lout, aes(x=logor1,as.numeric(id)+0.4), colour="red") +
geom_errorbarh(data=lout, aes(xmin=lb3,xmax=ub3, y=as.factor(id)),colour="blue", height=0.01) +
geom_errorbarh(data=lout, aes(xmin=lb2,xmax=ub2, y=as.numeric(id)+0.2),colour="green", height=0.01) +
geom_errorbarh(data=lout, aes(xmin=lb1,xmax=ub1, y=as.numeric(id)+0.4),colour="red", height=0.01) +
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

png(file="/well/todd/users/jinshaw/output/aad/under_7/hla/redo/multinom_all_hla_adjusted_dr34_3_n.png", res=800,
width=25, height=20, units="cm")
grid.arrange(one_out,two_out,ncol=2)
dev.off()


