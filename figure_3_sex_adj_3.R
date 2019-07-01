#figure_3_sex_adj_3.R

#faster (but still quite slow!) way of generating Figure 3

library(ggplot2)
library(ggbio)
library(GUESSFM)
library(R2GUESS)
library(reshape2)
library(ggplot2)
library(snpStats)
library(speedglm)
library(jimisc)
library(GenomicRanges)
library(stringr)
library(coloc)

#define regions:
likelihoods<-read.table(file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_3.txt",
header=T, sep="\t", as.is=T)
likelihoods$pfdr<-p.adjust(likelihoods$p,method="BH")
likelihoods<-likelihoods[likelihoods$pfdr<0.1,]
hits<-likelihoods$id

#set paths
mydir <-"/well/todd/users/jinshaw/aad/under_7/guessfm/uk/"
outdir<-"/well/todd/users/jinshaw/output/aad/under_7/guessfm/uk/"


e<-read.table(file="/well/todd/users/jinshaw/eqtl/eqtlgen/cis-eQTLs_full_20180905.txt.gz",header=T,as.is=T)
genes<-read.table(file="/well/todd/users/jinshaw/aad/under_7/refFlat.txt.gz",as.is=T)
genes<-genes[!duplicated(genes$V1),]
genes<-GRanges(seqnames=genes$V3,
ranges=IRanges(genes$V5,end=genes$V6),
gene=genes$V1)


getfigs<-function(snp, gene){

message(paste0("READING IN FINE MAPPING RESULTS"))
load(file=paste0(mydir, snp,"_sexadj/data.RData"))
mydir<-paste0(mydir,snp,"_sexadj/")
colnames(DATA)<-gsub(":",".",colnames(DATA))
colnames(DATA)<-gsub("<",".",colnames(DATA))
colnames(DATA)<-gsub(">",".",colnames(DATA))

colnames(DATA)<-ifelse(substr(colnames(DATA),1,1) %in% c("1","2","3","4","5","6","7","8","9"),
paste0("X",colnames(DATA)),colnames(DATA))

chromosome<-substr(snp,5,6)
chromosome<-ifelse(substr(chromosome,2,2)=="_",substr(chromosome,1,1),chromosome)

load(file=paste0(mydir,"summx.RData"))
summx$tag<-as.character(summx$tag)
summx$tag<-gsub("\\..*","",summx$tag)
summx$tag<-as.factor(summx$tag)
le<-length(unique(summx$tag))
summx$tag<-reorder(summx$tag,-summx$ppsum)
hues = c("red","blue","green","yellow","turquoise1","orange","purple","slategrey1","tan","pink","lightblue")
hues =hues[1:le]
#seq(from=10, to=(10+((le-1)*22)), by=22)
names(hues)<-levels(summx$tag)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Homo.sapiens)

g<-GRanges(seqnames=c(paste0("chr",chromosome)),
ranges=IRanges(min(summx$position),end=max(summx$position)))
o1<-subsetByOverlaps(genes,g)
o1<-as.data.frame(o1)
if (snp=="imm_17_35306733"){
o1<-o1[o1$gene %in% c("IKZF3", "ORMDL3", "GSDMB"),]
}
if(nrow(o1)>0){
m<-(0.2*nrow(o1))
o1$y=c(seq(0.2,m,0.2))
}
if(nrow(o1)==0){
g<-GRanges(seqnames=c(paste0("chr",chromosome)),
ranges=IRanges(min(summx$position)-200000,end=max(summx$position)+100000))
o1<-subsetByOverlaps(genes,g)
o1<-as.data.frame(o1)
m<-(0.2*nrow(o1))
o1$y=c(seq(0.2,m,0.2))
}
huesg = c("green","purple","turquoise1","orange","slategrey1","tan","pink","lightblue")
huesg =huesg[1:nrow(o1)]
o1$gene<-as.factor(o1$gene)
names(huesg)<-levels(o1$gene)
t1<-ggplot(data=o1, aes(x=start, xend=end, y=y, yend=y,colour=as.factor(gene))) + geom_segment(arrow=arrow(length = unit(0.03, "npc"))) +
theme(axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.title.y=element_blank(),
axis.text.x=element_blank(),
axis.title.x=element_blank(),
axis.ticks.x=element_blank()) +
scale_colour_manual(name="Gene name",values=c(huesg)) + scale_y_continuous(limits=c(0, (m+0.2))) +
theme(legend.text=element_text(size=6),  
legend.title=element_text(size=6))


e1<-e[e$SNPChr==chromosome & e$SNPPos>min(summx$position)-50000 & e$SNPPos< max(summx$position)+50000,]

e1$logzscore<-ifelse(e1$Zscore>0,log10(e1$Zscore),
ifelse(e1$Zscore<=0, log10(-e1$Zscore)*-1,NA))
if(snp %in% c("imm_17_35306733")){
e1<-e1[e1$Pvalue<5*10^-150,]
}
if(snp %in% c("imm_15_77022012","imm_1_205006527","imm_20_1564206")){
e1<-e1[e1$Pvalue<5*10^-50,]
}
if(snp %in% c("imm_10_6170083","imm_6_128335625")){
e1<-e1[e1$Pvalue<5*10^-25,]
}
if(snp %in% c("imm_9_4280823","imm_4_123335627","imm_14_97557760")){
e1<-e1[e1$Pvalue<5*10^-8,]
}
summx$SNPPos<-summx$position

direction<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/results/",snp,"_effects_uk_sex_adj_3.txt"),header=T,
sep="\t")
direction$se<-sqrt(direction$Var.beta)
direction$lb<-direction$beta-(qnorm(0.975)*direction$se)
direction$ub<-direction$beta+(qnorm(0.975)*direction$se)
direction<-merge(summx,direction,by="snp",all.x=T)

e1<-merge(e1,direction, by="SNPPos",all.x=T)

e1$logzscore<-ifelse(e1$OtherAllele==e1$effect & !is.na(e1$effect), e1$logzscore*-1,e1$logzscore)

#make the credible SNPs the same colours as they are in the guessfm analysis (might be tricky) but think will make much clearer):
e1$tag1<-factor(e1$tag,levels=c(levels(e1$tag),"nogroup"))
e1[is.na(e1$tag),"tag1"]<-"nogroup"
e1$credible<-ifelse(!is.na(e1$a1),1,0)

hues1 = c(hues,"black")
names(hues1)<-levels(e1$tag1)


direct1<-ggplot(data=direction, aes(position, beta,colour=as.factor(tag))) + geom_point() +
geom_errorbar(data=direction, aes(ymin=lb, ymax=ub, x=position,colour=as.factor(tag)),size=0.1) +
scale_y_continuous(name="T1D (<7) log-odds ratio") +
#scale_x_continuous(name=paste0("Position along chromosome ",chromosome)) +
geom_hline(yintercept=0, colour="red", linetype="dashed") +
scale_colour_manual(name="Tag group", values = c(hues)) +
theme(legend.text=element_text(size=6),  
legend.title=element_text(size=6),
axis.title.x=element_text(size=7),
axis.title.y=element_text(size=9))



#Now colocalisation:
message("Performing the T1D association analyses...")
#DISEASE:
load(file="/well/todd/users/jinshaw/aad/under_7/pheno_mult_uk_3.R")
pheno<-pheno[(pheno$group==0|pheno$group==1) & !is.na(pheno$group),]
pheno<-pheno[!is.na(pheno$sex),]

load(file=paste0(mydir, "data.RData"))
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


#EQTL:
message("Gathering eQTL results...")
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
load(file=paste0(mydir, "/data.RData"))
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
dfb1<-merge(dfb,lds,by="SNP")
colnames(dfb1)[ncol(dfb1)]<-"ld"
ou1<-ggplot(data=dfb1[dfb1$trait==paste0(gene," eQTL"),], aes(SNPPos, absz, colour=ld)) + geom_point() +
scale_y_continuous(name="Absolute eQTL z-score") +
scale_colour_continuous(name=paste0("LD with top \n",gene," eQTL")) +
theme(legend.text=element_text(size=6),  
legend.title=element_text(size=6), 
axis.title.x=element_text(size=7), 
axis.title.y=element_text(size=9))
ou2<-ggplot(data=dfb1[dfb1$trait==paste0("T1D"),], aes(SNPPos, absz, colour=ld)) + geom_point() +
scale_y_continuous(name="Absolute T1D z-score") +
scale_colour_continuous(name=paste0("LD with top \n",gene," eQTL")) +
scale_x_continuous(breaks=c(min(o1$start,o1$end,min(e1$SNPPos)), max(e1$SNPPos))) +
theme(legend.text=element_text(size=6),
legend.title=element_text(size=6), 
axis.title.x=element_text(size=7), 
axis.title.y=element_text(size=9))


trying2<-ggplot(data=e1, aes(SNPPos, logzscore, colour=tag1,shape=as.factor(GeneSymbol))) + geom_point() +
scale_y_continuous(name="eQTL log z-score") +
scale_shape_discrete(name="Gene") +
scale_x_continuous(breaks=c(min(o1$start,o1$end,min(e1$SNPPos)), max(e1$SNPPos))) +
scale_colour_manual(name="Tag group",values = c(hues1),guide=FALSE) +
geom_hline(yintercept=0, colour="red",linetype="dashed") +
theme(legend.text=element_text(size=6),
legend.title=element_text(size=6), 
axis.title.x=element_text(size=7), 
axis.title.y=element_text(size=9))

#COMBINE ALL THE ABOVE:
message("Combining...")
t<-tracks(ou1, ou2,t1,direct1,trying2, xlab=paste0("Position along chromosome ",chromosome),heights=c(0.9,0.9,0.6,0.9,1), xlim=c(min(o1$start,o1$end,min(e1$SNPPos)), max(e1$SNPPos)))
ggsave(t,file=paste0(outdir,snp,"_sexadj/finemap_plus_coloc_plus_directions_sexadj_3.png"),
dpi=800, height=20, width=10, units="cm")
ggsave(t,file=paste0(outdir,snp,"_sexadj/finemap_plus_coloc_plus_directions_sexadj_3.pdf"),
height=12, width=7)
message("DONE!")
return(t)
}

ptprk<-getfigs("imm_6_128335625", gene="THEMIS")
ctsh<-getfigs("imm_15_77022012", gene="CTSH")
ikzf3<-getfigs("imm_17_35306733", gene="IKZF3")

