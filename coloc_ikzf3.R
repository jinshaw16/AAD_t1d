#coloc_ikzf3.R
#looking at colocalisation of disease association and whole blood eQTL for IKZF3, ORMDL3 and GASDMB
library(stringr)
library(snpStats)
library(coloc)
library(ggplot2)
#load the disease associations with the <7s:

mydir <-"/well/todd/users/jinshaw/aad/under_7/guessfm/uk/"
outdir<-"/well/todd/users/jinshaw/output/aad/under_7/guessfm/uk/"

load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_uk_2.R")
pheno<-pheno[(pheno$group==0|pheno$group==1) & !is.na(pheno$group),]


#and the imputed data for this region:
getdis<-function(snp){
load(file=paste0(mydir,snp, "/data.RData"))
#Generate summary stats for each SNP in the 0.5 Mb region:
DATA<-DATA[rownames(pheno),]
colnames(DATA)<-gsub(":",".",colnames(DATA))
colnames(DATA)<-gsub("<",".",colnames(DATA))
colnames(DATA)<-gsub(">",".",colnames(DATA))
colnames(DATA)<-ifelse(substr(colnames(DATA),1,1) %in% c("1","2","3","4","5","6","7","8","9"),
paste0("X",colnames(DATA)),colnames(DATA))
#get the alleles:
system(paste0("awk -F \' \' \' {print $1,$2,$3,$4,$5}\' /well/todd/users/jinshaw/aad/under_7/imputation/",
snp,"_n_out > /well/todd/users/jinshaw/aad/under_7/imputation/",snp,"alleles"))
alls<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/",snp,"alleles"),header=F,as.is=T)
alls$V2<-gsub(":",".",alls$V2)
alls$V2<-gsub("<",".",alls$V2)
alls$V2<-gsub(">",".",alls$V2)
alls$V2<-ifelse(substr(alls$V2,1,1) %in% c("1","2","3","4","5","6","7","8","9"),
paste0("X",alls$V2),alls$V2)
#load the whole blood eqtls for the 3 genes in the same region:
alls<-alls[alls$V2 %in% colnames(DATA),]
alls<-alls[!duplicated(alls$V2),]
rownames(alls)<-alls$V2
alls<-alls[colnames(DATA),]

s<-snp.rhs.estimates(formula=group ~ PC1 + PC2 + PC3 + PC4 + PC5, data=pheno, link="logit", family="binomial",
snp.data=DATA, uncertain=TRUE)
s<-do.call("rbind",s)
s1<-as.data.frame(s)
for(i in 1:4){
s1[,i]<-unlist(s1[,i])
}
s1<-s1[alls$V2,]
s1$ref=alls$V4
s1$effect<-alls$V5
s1$SNPPos<-alls$V3

cs<-col.summary(DATA)
s1$MAF=cs$MAF
s1$type="cc"
s1$snp<-gsub("\\..*","",rownames(s1))
s1$z<-s1$beta/sqrt(s1$Var.beta)
s1$p<-2*pnorm(-abs(s1$z))
s1<-s1[order(s1$SNP,s1$p),]
s1<-s1[!duplicated(s1$SNP),]
return(s1)
}


#load the whole blood eqtls for the 3 genes in the same region:
testcoloc<-function(snp, gene, plot=TRUE){
message(paste0("Testing colocalisation of disease and ",gene,"\n"))
s1<-getdis(snp)
k<-system(paste0("zgrep ",gene," /well/todd/users/jinshaw/eqtl/eqtlgen/cis-eQTLs_full_20180905.txt.gz"),intern=T)

k1<-str_split_fixed(k, "\t",14)
k1<-as.data.frame(k1)
k1$V1<-as.character(k1$V1)
k1$V1<-gsub("E","e",k1$V1)
k1$V1<-as.numeric(k1$V1)
colnames(k1)<-c("Pvalue","SNP","SNPChr","SNPPos","Zscore","AssessedAllele",
"OtherAllele","Gene","GeneSymbol","GeneChr","GenePos","NrCohorts","NrSamples","FDR")

for(vars in c("SNPChr","Zscore","SNPPos","GeneChr","GenePos","NrCohorts","NrSamples","FDR")){
k1[,vars]<-as.character(k1[,vars])
k1[,vars]<-as.numeric(k1[,vars])
}

for(vars in c("SNP","AssessedAllele","OtherAllele","Gene","GeneSymbol")){
k1[,vars]<-as.character(k1[,vars])
}

b<-merge(s1,k1,by="SNPPos")
b<-b[order(b$SNPPos),]
s1<-s1[s1$SNPPos %in% b$SNPPos,]
s1<-s1[order(s1$SNPPos),]
k1<-k1[k1$SNPPos %in% b$SNPPos,]
k1<-k1[order(k1$SNPPos),]
#get alligned to correct allele:
wh<-which(k1$OtherAllele==s1$ref & k1$AssessedAllele==s1$effect)
wh1<-which(k1$OtherAllele==s1$effect & k1$AssessedAllele==s1$ref)
s1[wh1,"beta"]<-s1[wh1,"beta"]*-1
s1[wh1,"z"]<-s1[wh1,"z"]*-1
wall<-c(wh,wh1)
s1<-s1[wall,]
k1<-k1[wall,]
k1$N<-k1$NrSamples
k1$pvalues<-k1$Pvalue
k1$MAF<-s1$MAF
k1$type="quant"
k1$snp=k1$SNP
s1$pvalues=s1$p
s1$s=nrow(pheno[pheno$affected==2,])/nrow(pheno)
df1=s1[,c("snp","pvalues","type","N","MAF","s")]
df2=k1[,c("snp","pvalues","type","N","MAF")]
rownames(df2)<-rownames(df1)
df2$SNP<-df1$SNP
coloc<-coloc.abf(df1, df2, p1 = 1e-04, p2 = 1e-04,
      p12 = 1e-05)
probs<-coloc[[2]]
l<-s1[,c("snp","SNPPos")]
probs<-merge(probs,l,by="snp")
if (plot==TRUE){
ggplot(data=probs, aes(SNPPos, SNP.PP.H4)) + geom_point()
}
return(coloc)
}
ikz<-testcoloc(snp="imm_17_35306733","IKZF3")
orm<-testcoloc(snp="imm_17_35306733","ORMDL3")
gsdm<-testcoloc(snp="imm_17_35306733","GSDMB")

ctsh<-testcoloc(snp="imm_15_77022012", "CTSH")

glis3<-testcoloc(snp="imm_9_4280823", "GLIS3")
