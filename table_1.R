#table_1.R
#produce table which summarises patient characteristics (discovery and replication):

library(plyr)
library(dplyr)

#load discovery data info:
d<-"/well/todd/users/jinshaw/t1d_risk/immunochip/"
#read SNP and phenotype data:
load(file=paste0(d,"all_inds_unrel_postqc_3.RData"))
pheno<-all@samples
pheno$onset<-as.numeric(pheno$onset)
pheno$group<-ifelse(pheno$affected==1,0,
ifelse(pheno$onset<7 & !is.na(pheno$onset),1,
ifelse(pheno$onset>=7 & pheno$onset<13 & !is.na(pheno$onset),2,
ifelse(pheno$onset>=13 & !is.na(pheno$onset),3,NA))))
pheno<-pheno[!is.na(pheno$group)| pheno$affected==1,]
pheno<-pheno[!is.na(pheno$sex) & pheno$sex!=0,]


g<-pheno %>%
group_by(group)	%>%
summarise(mon=mean(onset),N=n(),
med=median(onset),lb=quantile(onset,na.rm=T)[2],
ub=quantile(onset,na.rm=T)[4],
min=min(onset,na.rm=T),max=max(onset,na.rm=T))

sink(file="/well/todd/users/jinshaw/output/aad/under_7/summary/table_1_n_1.txt")
cat(paste0("Discovery;Controls;<7;7-13;>13\n"))

cat(paste0("N;",g[1,3],";",g[2,3],";",g[3,3],";",g[4,3],"\n"))
cat(paste0("Median (IQR) age-at-diagnosis;",g[1,4]," (",g[1,5],", ",g[1,6],") [",g[1,7],", ",g[1,8],"];",
g[2,4]," (",g[2,5],", ",g[2,6],") [",g[2,7],", ",g[2,8],"];",
g[3,4]," (",g[3,5],", ",g[3,6],") [",g[3,7],", ",g[3,8],"];",
g[4,4]," (",g[4,5],", ",g[4,6],") [",g[4,7],", ",g[4,8],"]\n"))

addcat<-function(df,cat,nums){
t<-table(df[,cat],df$group)

for (i in nums){
cat(paste0(rownames(t)[i],";",t[i,1]," (",round(100*t[i,1]/(sum(t[,1])),digits=1),"%);",
t[i,2]," (",round(100*t[i,2]/(sum(t[,2])),digits=1),"%);",
t[i,3]," (",round(100*t[i,3]/(sum(t[,3])),digits=1),"%);",
t[i,4]," (",round(100*t[i,4]/(sum(t[,4])),digits=1),"%)\n"))
}
}
addcat("sex",nums=2,df=pheno)
addcat("country",nums=c(1:6),df=pheno)
sink()

