#coloc_ikzf3_sexadj_3.R
#looking at colocalisation of disease association and whole blood eQTL for IKZF3, ORMDL3 and GASDMB
library(stringr)
library(snpStats)
#library(coloc)
library(ggplot2)
library(ggbio)
#load the disease associations with the <7s:

mydir <-"/well/todd/users/jinshaw/aad/under_7/guessfm/uk/"
outdir<-"/well/todd/users/jinshaw/output/aad/under_7/guessfm/uk/"

load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_uk_3.R")
pheno<-pheno[(pheno$group==0|pheno$group==1) & !is.na(pheno$group),]
pheno<-pheno[!is.na(pheno$sex),]

#and the imputed data for this region:
getdis<-function(snp){
load(file=paste0(mydir,snp, "_sexadj/data.RData"))
#Generate summary stats for each SNP in the 0.5 Mb region:
DATA<-DATA[rownames(pheno),]
colnames(DATA)<-gsub(":",".",colnames(DATA))
colnames(DATA)<-gsub("<",".",colnames(DATA))
colnames(DATA)<-gsub(">",".",colnames(DATA))
colnames(DATA)<-ifelse(substr(colnames(DATA),1,1) %in% c("1","2","3","4","5","6","7","8","9"),
paste0("X",colnames(DATA)),colnames(DATA))
#get the alleles:
system(paste0("zcat /well/todd/users/jinshaw/aad/under_7/imputation/",
snp,"_3_out.gz | awk -F \' \' \' {print $1,$2,$3,$4,$5}\' > /well/todd/users/jinshaw/aad/under_7/imputation/",snp,"alleles"))
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

s<-snp.rhs.estimates(formula=group ~ PC1 + PC2 + PC3 + PC4 + PC5 + sex, data=pheno, link="logit", family="binomial",
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

#editting the coloc package slightly as want to just use the z scores (this is more accurate with the minimum floating point precision number)
approx.bf.z <- function(z,f,type, N, s, suffix=NULL) {
  if(type=="quant") {
    sd.prior <- 0.15
    V <- Var.data(f, N)
  } else {
    sd.prior <- 0.2
    V <- Var.data.cc(f, N, s)
  }
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- sd.prior^2 / (sd.prior^2 + V)
  ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
  lABF = 0.5 * (log(1-r) + (r * z^2))
  ret <- data.frame(V,z,r,lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  return(ret)  
}

logdiff <- function(x,y) {
  my.max <- max(x,y)                              ##take out the maximum value in log form
  my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  return(my.res)
}

logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}

Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
}

Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}

sdY.est <- function(vbeta, maf, n) {
    warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
    oneover <- 1/vbeta
    nvx <- 2 * n * maf * (1-maf)
    m <- lm(nvx ~ oneover - 1)
    cf <- coef(m)[['oneover']]
    if(cf < 0)
        stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
    return(sqrt(cf))
}

combine.abf <- function(l1, l2, p1, p2, p12) {
  lsum <- l1 + l2
  lH0.abf <- 0
  lH1.abf <- log(p1) + logsum(l1)
  lH2.abf <- log(p2) + logsum(l2)
  lH3.abf <- log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(lsum))
  lH4.abf <- log(p12) + logsum(lsum)

  all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
  my.denom.log.abf <- logsum(all.abf)
  pp.abf <- exp(all.abf - my.denom.log.abf)
  names(pp.abf) <- paste("PP.H", (1:length(pp.abf)) - 1, ".abf", sep = "")
  print(signif(pp.abf,3))
  print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
  return(pp.abf)
}

approx.bf.estimates <- function (z, V, type, suffix=NULL, sdY=1) {
  sd.prior <- if (type == "quant") { 0.15*sdY } else { 0.2 }
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}

process.dataset <- function(d, suffix) {
  #message('Processing dataset')

  nd <- names(d)
  if (! 'type' %in% nd)
    stop("dataset ",suffix,": ",'The variable type must be set, otherwise the Bayes factors cannot be computed')

  if(!(d$type %in% c("quant","cc")))
      stop("dataset ",suffix,": ","type must be quant or cc")
  
  if(d$type=="cc") {
      if(! "s" %in% nd)
          stop("dataset ",suffix,": ","please give s, proportion of samples who are cases")
      if("pvalues" %in% nd && !( "MAF" %in% nd))
          stop("dataset ",suffix,": ","please give MAF if using p values")
      if(d$s<=0 || d$s>=1)
          stop("dataset ",suffix,": ","s must be between 0 and 1")
  }
  
  if(d$type=="quant") {
      if(!("sdY" %in% nd || ("MAF" %in% nd && "N" %in% nd )))
          stop("dataset ",suffix,": ","must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
  }
  
  if("beta" %in% nd && "varbeta" %in% nd) {  ## use beta/varbeta.  sdY should be estimated by now for quant
    if(length(d$beta) != length(d$varbeta))
      stop("dataset ",suffix,": ","Length of the beta vectors and variance vectors must match")
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$beta))
    if(length(d$snp) != length(d$beta))
      stop("dataset ",suffix,": ","Length of snp names and beta vectors must match")
 
    if(d$type=="quant" && !('sdY' %in% nd)) 
          d$sdY <- sdY.est(d$varbeta, d$MAF, d$N)
   df <- approx.bf.estimates(z=d$beta/sqrt(d$varbeta),
                              V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
    df$snp <- as.character(d$snp)
    return(df)
  }

  if("z" %in% nd & "MAF" %in% nd & "N" %in% nd) { ## no beta/varbeta: use p value / MAF approximation
    if (length(d$z) != length(d$MAF))
      stop('Length of the z-value vectors and MAF vector must match')
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$z))
    df <- data.frame(z = d$z,
                     MAF = d$MAF,
                     snp=as.character(d$snp))    
    colnames(df)[-3] <- paste(colnames(df)[-3], suffix, sep=".")
#   df <- subset(df, df$MAF>0 & df$z>0) # all p values and MAF > 0
    abf <- approx.bf.z(z=df$z, f=df$MAF, type=d$type, N=d$N, s=d$s, suffix=suffix)
    df <- cbind(df, abf)
    return(df)  
  }
  else{
  stop("Must give, as a minimum, one of:\n(beta, varbeta, type, sdY)\n(beta, varbeta, type, MAF)\n(z, MAF, N, type)")
}
}

coloc.abf <- function(dataset1, dataset2, MAF=NULL, 
                      p1=1e-4, p2=1e-4, p12=1e-5) {

  if(!is.list(dataset1) || !is.list(dataset2))
    stop("dataset1 and dataset2 must be lists.")
  if(!("MAF" %in% names(dataset1)) & !is.null(MAF))
    dataset1$MAF <- MAF
  if(!("MAF" %in% names(dataset2)) & !is.null(MAF))
    dataset2$MAF <- MAF
  
  df1 <- process.dataset(d=dataset1, suffix="df1")
  df2 <- process.dataset(d=dataset2, suffix="df2")
  merged.df <- merge(df1,df2)

   if(!nrow(merged.df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")

  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
  ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  my.denom.log.abf <- logsum(merged.df$internal.sum.lABF)
  merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)
  
 
############################## 

  pp.abf <- combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)  
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.abf)
  
  output<-list(summary=results, results=merged.df)
  return(output)
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
k1$z<-abs(k1$Zscore)
s1$pvalues=s1$p
s1$z<-abs(s1$z)
s1$s=nrow(pheno[pheno$affected==2,])/nrow(pheno)
df1=s1[,c("snp","z","type","N","MAF","s")]
df2=k1[,c("snp","z","type","N","MAF")]
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

plotit<-function(snp, gene){
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
k1$absz<-abs(k1$Zscore)
s1$pvalues=s1$p
s1$absz<-abs(s1$z)
df1=s1[,c("snp","pvalues","type","N","MAF","absz")]
df1$trait="T1D"
df2=k1[,c("snp","pvalues","type","N","MAF","absz")]
df2$trait=paste0(gene," eQTL")
rownames(df2)<-rownames(df1)
df2$SNP<-rownames(df2)
df1$SNP<-rownames(df1)
dfb<-rbind(df1,df2)
l<-s1[,c("snp","SNPPos")]
dfb<-merge(dfb,l,by="snp")
dfb$SNP<-gsub("\\.",":",dfb$SNP)
load(file=paste0(mydir,snp, "_sexadj/data.RData"))
DATA<-DATA[,colnames(DATA) %in% dfb$SNP]
top<-dfb[dfb$trait==paste0(gene," eQTL"),]
top<-top[top$absz==max(top$absz),]
tsnp<-top$SNP
lds<-ld(DATA[,!colnames(DATA) %in% tsnp],DATA[,tsnp],stats="R.squared")
next1<-data.frame(tsnp=1)
colnames(next1)=tsnp
lds<-rbind(lds,next1)
rownames(lds)[nrow(lds)]<-tsnp
lds$SNP=rownames(lds)
o<-ggplot(data=dfb, aes(SNPPos, absz, colour=as.factor(trait))) + geom_point()
dfb1<-merge(dfb,lds,by="SNP")
colnames(dfb1)[ncol(dfb1)]<-"ld"
o1<-ggplot(data=dfb1[dfb1$trait==paste0(gene," eQTL"),], aes(SNPPos, absz, colour=ld)) + geom_point() +
scale_y_continuous(name="Absolute eQTL z-score") +
scale_colour_continuous(name=paste0("LD with top \n",gene," eQTL"))
o2<-ggplot(data=dfb1[dfb1$trait==paste0("T1D"),], aes(SNPPos, absz, colour=ld)) + geom_point() +
scale_y_continuous(name="Absolute T1D z-score") +
scale_colour_continuous(name=paste0("LD with top \n",gene," eQTL"))
t<-tracks(o1,o2, xlab=paste0("Position along chromosome ",k1$GeneChr[1]))
ggsave(t, file=paste0("/well/todd/users/jinshaw/output/aad/under_7/coloc/",snp,"_",gene,"_sexadj_3.png"),
width=25,height=20, dpi=800, units="cm")
ggsave(t, file=paste0("/well/todd/users/jinshaw/output/aad/under_7/coloc/",snp,"_",gene,"_sexadj_3.pdf"),
width=10, height=7)
}

il10<-testcoloc(snp="imm_1_205006527", "IL10")
plotit(snp="imm_1_205006527","IL10")
il24<-testcoloc(snp="imm_1_205006527", "IL24")
plotit(snp="imm_1_205006527","IL24")
faim3<-testcoloc(snp="imm_1_205006527", "FAIM3")
plotit(snp="imm_1_205006527","FAIM3")
ikz<-testcoloc(snp="imm_17_35306733","IKZF3")
plotit(snp="imm_17_35306733", "IKZF3")
orm<-testcoloc(snp="imm_17_35306733","ORMDL3")
plotit(snp="imm_17_35306733", "ORMDL3")
gsdm<-testcoloc(snp="imm_17_35306733","GSDMB")
plotit(snp="imm_17_35306733", "GSDMB")

ctsh<-testcoloc(snp="imm_15_77022012", "CTSH")
plotit(snp="imm_15_77022012", "CTSH")

glis3<-testcoloc(snp="imm_9_4280823", "GLIS3")
plotit(snp="imm_9_4280823", "GLIS3")

themis<-testcoloc(snp="imm_6_128335625","THEMIS")
plotit(snp="imm_6_128335625", "THEMIS")
ptprk<-testcoloc(snp="imm_6_128335625","PTPRK")
plotit(snp="imm_6_128335625", "PTPRK")

