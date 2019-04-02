#hla_check.R
library(ggplot2)
#Comparing the hla from genotype to imputed in the case-contorl cohort.


hlag<-read.table(file="/well/todd/users/jinshaw/t1d_risk/immunochip/cc-hla-2015-11-13.tab",header=T)
#get the lookup files:
lookup1<-read.table(file="/well/todd/users/jinshaw/t1d_risk/immunochip/support/ic-cc-clean-lookup-2013-02-27.tab",header=T,as.is=T,comment.char="")
colnames(lookup1)<-c("collection","member","uniqueID","dilID")
hlag<-merge(hlag, lookup1, by="uniqueID")
hlag$sample.id=hla$uniqueID


for (i in c("DQA1_4digit_1", "DQA1_4digit_2",
"HLAB_4digit_1","HLAB_4digit_2",
"DRB1_4digit_1","DRB1_4digit_2",
"DQB1_4digit_1","DQB1_4digit_2",
"HLAA_4digit_1","HLAA_4digit_2")
){
hlag[,i]<-ifelse(nchar(hlag[,i])==3,paste0("0",hlag[,i]),
ifelse(nchar(hlag[,i])<3,NA,hlag[,i]))
hlag[,i]<-ifelse(!is.na(hlag[,i]),paste0(substr(hlag[,i],1,2),":",substr(hlag[,i],3,4)),NA)
}


load(file="/well/todd/users/jinshaw/aad/under_7/hla_all_2.RData")
lookup1$sample.id<-paste0(lookup1$collection,".",lookup1$member)
hlaout<-merge(hla,lookup1,by="sample.id")
hlaout$uniqueID=hlaout$V3
hlag$uniqueID.y=hlag$uniqueID
hlaout<-merge(hlaout,hlag,by="uniqueID.y")

doinitialprop<-function(name1,name2){
dummy<-ifelse((hlaout[,paste0(name1,"1")]==hlaout[,paste0(name2,"1")]  & hlaout[,paste0(name1,"2")]==hlaout[,paste0(name2,"2")]) |
(hlaout[,paste0(name1,"1")]==hlaout[,paste0(name2,"2")]  & hlaout[,paste0(name1,"2")]==hlaout[,paste0(name2,"1")]),1,
ifelse(is.na(hlaout[,paste0(name1,"1")]) | is.na(hlaout[,paste0(name1,"2")]) | is.na(hlaout[,paste0(name2,"1")]) | is.na(hlaout[,paste0(name2,"2")]),NA,0))
t1<-data.frame(locus=name1,match=dummy)
return(t1)
}
ts<-mapply(doinitialprop, name1=c("A","B","DRB1", "DQB1","DQA1"), 
name2=c("HLAA_4digit_","HLAB_4digit_","DRB1_4digit_", "DQB1_4digit_","DQA1_4digit_"),SIMPLIFY=F)
ts<-do.call("rbind",ts)

p<-ts[!is.na(ts$match),]
t1<-table(p$locus)
png(file="/well/todd/users/jinshaw/output/aad/under_7/hla/redo_1/concordance_4digits_n.png",
height=20,width=20, units="cm", res=400)
ggplot(data=ts[!is.na(ts$match),], aes(as.factor(locus), fill=as.factor(match))) + geom_bar(position="fill",na.rm =TRUE) +
scale_fill_hue(name="HLA 4 digit\n concordance", breaks=c(0,1),labels=c("FALSE","TRUE")) +
scale_x_discrete(name="Locus") +
scale_y_continuous(name="Proportion concordance on both chromosomes") +
annotate("text",x=c(1,2,3,4,5), y=c(0.1,0.1,0.1,0.1,0.1), label=c(paste0("N=",t1)))
dev.off()


#and without DQA1:
ts<-mapply(doinitialprop, name1=c("A","B","DRB1", "DQB1"),
name2=c("HLAA_4digit_","HLAB_4digit_","DRB1_4digit_", "DQB1_4digit_"),SIMPLIFY=F)
ts<-do.call("rbind",ts)

p<-ts[!is.na(ts$match),]
t1<-table(p$locus)
png(file="/well/todd/users/jinshaw/output/aad/under_7/hla/redo_1/concordance_4digits_nodqa1_n.png",
height=20,width=20, units="cm", res=400)
ggplot(data=ts[!is.na(ts$match),], aes(as.factor(locus), fill=as.factor(match))) + geom_bar(position="fill",na.rm =TRUE) +
scale_fill_hue(name="HLA 4 digit\n concordance", breaks=c(0,1),labels=c("FALSE","TRUE")) +
scale_x_discrete(name="Locus") +
scale_y_continuous(name="Proportion concordance on both chromosomes") +
annotate("text",x=c(1,2,3,4), y=c(0.1,0.1,0.1,0.1), label=c(paste0("N=",t1)))
dev.off()


png(file="/well/todd/users/jinshaw/output/aad/under_7/hla/redo_1/concordance_4digits_nodqa1_n_sm.png",
height=20,width=25, units="cm", res=800)
ggplot(data=ts[!is.na(ts$match),], aes(as.factor(locus), fill=as.factor(match))) + geom_bar(position="fill",na.rm =TRUE) +
scale_fill_hue(name="HLA 4 digit\n concordance", breaks=c(0,1),labels=c("FALSE","TRUE")) +
scale_x_discrete(name="Locus") +
scale_y_continuous(name="Proportion concordance on both chromosomes") +
annotate("text",x=c(1,2,3,4), y=c(0.1,0.1,0.1,0.1), label=c(paste0("N=",t1)))
dev.off()

