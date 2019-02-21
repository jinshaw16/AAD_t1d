#table_1.R
#produce table which summarises patient characteristics (discovery and replication):

library(plyr)
library(dplyr)

#load discovery data info:
load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult.R")
pheno$group<-ifelse(pheno$onset>=7 & pheno$onset<13 & !is.na(pheno$onset),2,ifelse(pheno$onset>=13 & !is.na(pheno$onset),3,pheno$group))
pheno<-pheno[!is.na(pheno$group),]



g<-pheno %>%
group_by(group)	%>%
summarise(mon=mean(onset),N=n())

sink(file="/well/todd/users/jinshaw/output/aad/under_7/summary/table_1.txt")
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


#and load replication data info:
load(file="/well/todd/users/jinshaw/aad/under_7/additionals_snps.RData")
load(file="/well/todd/users/jinshaw/aad/under_7/all_combined_alligned.RData")
k<-allgeno@samples
k<-k[k$uniqueID %in% pheno$uniqueID,]
k<-k[,c("uniqueID","country")]
k<-k[rownames(pheno),]
pheno$country<-k$country
pheno$group<-ifelse(pheno$onset<7,1,ifelse(pheno$onset>=7 & pheno$onset<13,2,ifelse(pheno$onset>=13,3,NA)))

pheno<-pheno[,!duplicated(colnames(pheno))]
g1<-pheno %>%
group_by(group) %>%
summarise(mon=mean(onset),N=n())

cat(paste0("N;0;",g1[1,3],";",g1[2,3],";",g1[3,3],"\n"))
cat(paste0("Mean age-at-diagnosis;NA;",round(g1[1,2],digits=2),
";",round(g1[2,2],digits=2),";",round(g1[3,2],digits=2),"\n"))


addcat1<-function(df,cat,nums){
t<-table(df[,cat],df$group)

for (i in nums){
cat(paste0(rownames(t)[i],";0 (0%);",t[i,1]," (",round(100*t[i,1]/(sum(t[,1])),digits=1),"%);",
t[i,2]," (",round(100*t[i,2]/(sum(t[,2])),digits=1),"%);",
t[i,3]," (",round(100*t[i,3]/(sum(t[,3])),digits=1),"%)\n"))
}
}

addcat1("sex",nums=2,df=pheno)
addcat1("country",nums=c(1:6),df=pheno)
sink()

