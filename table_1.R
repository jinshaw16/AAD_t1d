#table_1.R
#produce table which summarises patient characteristics (discovery and replication):

library(plyr)
library(dplyr)

#load discovery data info:
d<-"/well/todd/users/jinshaw/t1d_risk/immunochip/"
#read SNP and phenotype data:
load(file=paste0(d,"all_inds_unrel_postqc.RData"))
pheno<-all@samples
pheno$onset<-as.numeric(pheno$onset)
pheno$group<-ifelse(pheno$affected==1,0,
ifelse(pheno$onset<7 & !is.na(pheno$onset),1,
ifelse(pheno$onset>=7 & pheno$onset<13 & !is.na(pheno$onset),2,
ifelse(pheno$onset>=13 & !is.na(pheno$onset),3,NA))))
pheno<-pheno[!is.na(pheno$group)| pheno$affected==1,]



g<-pheno %>%
group_by(group)	%>%
summarise(mon=mean(onset),N=n())

sink(file="/well/todd/users/jinshaw/output/aad/under_7/summary/table_1_n.txt")
cat(paste0("Discovery;Controls;<7;7-13;>13\n"))

cat(paste0("N;",g[1,3],";",g[2,3],";",g[3,3],";",g[4,3],"\n"))
cat(paste0("Mean age-at-diagnosis;",round(g[1,2],digits=2),
";",round(g[2,2],digits=2),";",round(g[3,2],digits=2),
";",round(g[4,2],digits=2),"\n"))

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

