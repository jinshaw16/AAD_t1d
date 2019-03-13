#guessfm_readin_uk_only.R
## what files were created?

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

#define regions:
likelihoods<-read.table(file="/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_redo_n.txt",
header=T, sep="\t", as.is=T)
likelihoods<-likelihoods[likelihoods$p<(0.05/nrow(likelihoods)),]
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

#read in
readitallin<-function(snp){
load(file=paste0(mydir, snp,"/data.RData"))

mydir<-paste0(mydir,snp,"/")
list.files(mydir)

#create outdir if it doesn't exist:
if(!dir.exists(paste0(outdir,snp)))
dir.create(paste0(outdir,snp),recursive=T)


## read output with R2GUESS and run a basic qc plot
if (!file.exists(paste0(outdir,snp,"/model_diagnostics.png"))){
ess <- ess.read(f=paste0(mydir,"out_55000_55000"))
png(file=paste0(outdir,snp,"/model_diagnostics.png"),
res=200, height=20, width=20, units="cm")
k<-R2GUESS:::plot.ESS(ess)
dev.off()

png(file=paste0(outdir,snp,"/model_convergence.png"),
res=200, height=20, width=20, units="cm")
par(mfrow=c(1,2))
k1<-R2GUESS::check.convergence(ess)
dev.off()
}

## now read output with GUESSFM
if(file.exists(paste0(mydir,"out_55000_55000_sweeps_output_best_visited_models.txt"))){
system(paste0("mv ", mydir, "out_55000_55000_sweeps_output_best_visited_models.txt ",mydir,"out_55000_55000_output_best_visited_models.txt"))
system(paste0("mv ",mydir,"out_55000_55000_sweeps_output_marg_prob_incl.txt ", mydir,"out_55000_55000_output_marg_prob_incl.txt"))
}
dd <- read.snpmod(mydir)

#Examine posterior of number of SNPs in model:

png(file=paste0(outdir,snp,"/model_size.png"),
res=200, height=20, width=20, units="cm")
pp <- pp.nsnp(dd,plot=TRUE,expected=3, overdispersion = 1.00000001)
dev.off()

## examine the best models and SNPs with greatest marginal support within the tagged data.
best.models(dd)
best.snps(dd)

#Now expand these tag SNPs:
load(paste0(mydir,"/tags.RData"))
dx <- expand.tags(dd, tags)

#CTSH was behaving strangly, so handling this SNP slightly differently:
if(snp=="imm_15_77022012"){
best<-best.models(dx, pp.thr=0.001)
}
if(snp!="imm_15_77022012"){
best<-best.models(dx, cpp.thr = 0.99)
}

save(best, dx, file=paste0(mydir,"/best_models.RData"))
load(file=paste0(mydir,"/best_models.RData"))


colnames(DATA)<-gsub(":",".",colnames(DATA))
colnames(DATA)<-gsub("<",".",colnames(DATA))
colnames(DATA)<-gsub(">",".",colnames(DATA))

colnames(DATA)<-ifelse(substr(colnames(DATA),1,1) %in% c("1","2","3","4","5","6","7","8","9"),
paste0("X",colnames(DATA)),colnames(DATA))
Y<-as.matrix(Y)
#taking residuals from logistic regression including covriates as outcome... it's an approximation
abf <- abf.calc(y=Y,x=DATA,q=covariates,models=best$str,family="binomial", verbose=TRUE, approx.lm=T)

sm.all <- abf2snpmod(abf,expected=3, overdispersion = 1.00000001)
sp.all <-snp.picker(d=sm.all, data=DATA)

save(sm.all,sp.all,file=paste0(mydir,"/smsp_residuals.RData"))

load(file=paste0(mydir,"/smsp_residuals.RData"))
if (length(sp.all@groups)==0){
summx<-NULL
save(sm.all,sp.all,summx, file=paste0(mydir,"summx.RData"))
}
if(length(sp.all@groups)!=0){
#Grouping:

groups <- as(sp.all,"groups")

#for GLIS3, i think all should be part of the same group. just checking this here:
if (snp=="imm_9_4280823"){
test.groups<-function(tag1, tag2, groups){
w<-which(groups@tags %in% c(tag1, tag2))
groups@tags<-groups@tags[w]
groups@.Data<-list(groups@.Data[[w[1]]], groups@.Data[[w[2]]])
return(groups)
}
test<-test.groups(groups@tags[1], groups@tags[2], groups)
check.merge(sm.all, test)
#Looks like they might be the same signal as pp(all)<<pp(any)
g2<-groups.merge(groups, c(groups@tags[1], groups@tags[2]))
groups<-g2
}


#load the snp data:
chromosome<-substr(snp,5,6)
chromosome<-ifelse(substr(chromosome,2,2)=="_",substr(chromosome,1,1),chromosome)
l<-read.table(file=paste0("/well/1000G/WTCHG/1000GP_Phase3/1000GP_Phase3_chr",chromosome,".legend.gz"), header=T)
l$rsid<-gsub(":",".",l$id)
l$rsid<-ifelse((substr(l$rsid,1,1)==chromosome)|(substr(l$rsid,1,2)==chromosome),
paste0("X",l$rsid), l$rsid)
l$rsid<-gsub(">",".",l$rsid)
l$rsid<-gsub("<",".",l$rsid)
results<-l[l$rsid %in% colnames(DATA),]
rownames(results)<-results$rsid
summx <- guess.summ(sm.all,groups=groups,snps=results,position="position")
summx <- scalepos(summx,position="position")
summx$tag<-as.factor(summx$tag)

save(sm.all,sp.all,groups,results,summx, file=paste0(mydir,"summx.RData"))
load(file=paste0(mydir,"summx.RData"))


pat<-pattern.plot(sm.all,groups)
ggsave(pat,file=paste0(outdir,snp,"/cum_prob_plot.png"),
dpi=200, height=20, width=20, units="cm")
le<-length(unique(summx$tag))
summx$tag<-reorder(summx$tag,-summx$ppsum)
hues = c("red","blue","green","yellow","turquoise1","orange","purple","slategrey1","tan","pink","lightblue")
hues =hues[1:le]
#seq(from=10, to=(10+((le-1)*22)), by=22)
names(hues)<-levels(summx$tag)
signals<-signal.plot(summx) + scale_fill_manual(values=hues) + scale_colour_manual(values=hues)
chr<-ggchr(summx)  + scale_fill_manual(values=hues) + scale_colour_manual(values=hues)
lds<-ggld(DATA, summx)


write.table(summx, file=paste0(outdir,snp,"/credibles.txt"),col.names=T, row.names=F, quote=F, sep="\t")
#and bed files for Tony (in 37 and 38 builds):
s<-summx
s$chromosome=chromosome
s<-liftthem(s, snp.name="snp",chain="hg19ToHg38.over.chain",updateto="38")
s$posend<-s$position
s$posend38<-s$position38
s$chromosome<-paste0("chr",s$chromosome)
s1<-s[,c("chromosome","position","posend","snp")]
s2<-s[,c("chromosome","position38","posend38","snp")]
write.table(s1,file=paste0(outdir,snp,"/credibles_37.bed"),col.names=F, row.names=F, quote=F, sep="\t")
write.table(s2,file=paste0(outdir,snp,"/credibles_38.bed"),col.names=F, row.names=F, quote=F, sep="\t")
out<-tracks(chr,signals,lds, heights=c(1,2,1))
ggsave(out,file=paste0(outdir,snp,"/init_out.png"),
dpi=200, height=25, width=20, units="cm")

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
scale_colour_manual(name="Gene name",values=c(huesg)) + scale_y_continuous(limits=c(0, (m+0.2)))


data(genesymbol, package = "biovizBase")
g <- genesymbol[seqnames(genesymbol) == paste0('chr',chromosome)]
gr <- GRanges(
      seqnames = Rle(c(paste0("chr",chromosome)), c(1)),
      ranges = IRanges(min(summx$position), end = max(summx$position)))

o<-as.matrix(findOverlaps(g,gr))
g<-g[o[,1],]

if(length(g@seqnames)==0){
g <- genesymbol[seqnames(genesymbol) == paste0('chr',chromosome)]
gr <- GRanges(
      seqnames = Rle(c(paste0("chr",chromosome)), c(1)),
      ranges = IRanges(min(summx$position)-50000, end = max(summx$position)+50000))

o<-as.matrix(findOverlaps(g,gr))
g<-g[o[,1],]
}

t<-autoplot(Homo.sapiens, which = g)

if(snp!="imm_9_4281928"){
out1<-tracks(t,chr,signals,lds,heights=c(1,1,2,1))
}

e1<-e[e$SNPChr==chromosome & e$SNPPos>min(summx$position)-50000 & e$SNPPos< max(summx$position)+50000,]
if(snp=="imm_9_4281928"){
out1<-tracks(t,chr,signals,lds,heights=c(1,1,2,1),xlim=c(4280000,max(e1$SNPPos)))
}
ggsave(out1,file=paste0(outdir,snp,"/init_out_genes.png"),
dpi=200, height=30, width=20, units="cm")

e1$logzscore<-ifelse(e1$Zscore>0,log10(e1$Zscore),
ifelse(e1$Zscore<=0, log10(-e1$Zscore)*-1,NA))
if(snp %in% c("imm_17_35306733","imm_16_73809828")){
e1<-e1[e1$Pvalue<5*10^-150,]
}
if(snp %in% c("imm_1_170941164","imm_15_77022012","imm_1_205006527")){
e1<-e1[e1$Pvalue<5*10^-50,]
}
if(snp %in% c("imm_10_6170083")){
e1<-e1[e1$Pvalue<5*10^-25,]
}
if(snp %in% c("imm_9_4280823","imm_4_123335627","imm_6_128328079","imm_14_97557760")){
e1<-e1[e1$Pvalue<5*10^-8,]
}
summx$SNPPos<-summx$position

#run the disease comparison with this set of credible SNPs:
system(paste0("Rscript /users/todd/jinshaw/programs/aad/under_7/imputed_test_uk.R ",snp))

#and with the effect direction of the SNPs with disease for comparison:
direction<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/results/",snp,"_effects_uk.txt"),header=T,
sep="\t")
direction$se<-sqrt(direction$Var.beta)
direction$lb<-direction$beta-(qnorm(0.975)*direction$se)
direction$ub<-direction$beta+(qnorm(0.975)*direction$se)
direction<-merge(summx,direction,by="snp",all.x=T)
#export for the paper supplementary:
direction<-direction[direction$ppsum>0.9,]
ds<-direction[,c("snp","position","beta","ref","effect","pp","ppsum")]
ds$snp<-ifelse(substr(ds$snp,1,2)=="rs",gsub("\\..*","",ds$snp),ds$snp)
ds$ppsum<-round(ds$ppsum,digits=4)
ds$beta<-round(ds$beta,digits=4)
ds$pp<-format(ds$pp,scientific=T,digits=3)
write.table(ds,file=paste0("/well/todd/users/jinshaw/output/aad/under_7/guessfm/uk/",snp,"/supl_mat.txt"),col.names=T,row.names=F,
sep="\t",quote=F)

e1<-merge(e1,direction, by="SNPPos",all.x=T)

e1$logzscore<-ifelse(e1$OtherAllele==e1$effect & !is.na(e1$effect), e1$logzscore*-1,e1$logzscore)
trying<-ggplot(data=e1, aes(SNPPos, logzscore, colour=as.factor(GeneSymbol))) + geom_point() +
scale_y_continuous(name="Whole blood eQTL log z-score") +
scale_colour_hue(name="Gene") +
scale_x_continuous(name=paste0("Position along chromosome ",chromosome))

#make the credible SNPs the same colours as they are in the guessfm analysis (might be tricky) but think will make much clearer):
e1$tag1<-factor(e1$tag,levels=c(levels(e1$tag),"nogroup"))
e1[is.na(e1$tag),"tag1"]<-"nogroup"
e1$credible<-ifelse(!is.na(e1$a1),1,0)

hues1 = c(hues,"black")
names(hues1)<-levels(e1$tag1)

trying1<-ggplot(data=e1, aes(SNPPos, logzscore, colour=tag1,shape=as.factor(GeneSymbol))) + geom_point() +
scale_y_continuous(name="Whole blood eQTL log z-score") +
scale_shape_discrete(name="Gene") +
#scale_x_continuous(name=paste0("Position along chromosome ",chromosome)) +
theme(legend.position = c(0.5, 0.5)) + 
scale_colour_manual(name="Tag group",values = c(hues1))

if (snp!="imm_9_4280823"){
out2<-tracks(trying1, t, chr, signals,lds,heights=c(2,1,1,2,1))
}

if (snp=="imm_9_4280823"){
out2<-tracks(trying1, t, chr, signals,lds,heights=c(2,1,1,2,1), xlim=c(4280000,max(e1$SNPPos)))
}

ggsave(out2,file=paste0(outdir,snp,"/init_out_genes_eqtls.png"),
dpi=200, height=35, width=20, units="cm")

#and with the effect direction of the SNPs with disease for comparison:
direction<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/results/",snp,"_effects_uk.txt"),header=T,
sep="\t",as.is=T)
direction$se<-sqrt(direction$Var.beta)
direction$lb<-direction$beta-(qnorm(0.975)*direction$se)
direction$ub<-direction$beta+(qnorm(0.975)*direction$se)

direction<-merge(direction,summx,by="snp",all.x=T)
direction<-direction[direction$ppsum>0.9,]
direct<-ggplot(data=direction, aes(position, beta,colour=as.factor(tag))) + geom_point() +
geom_errorbar(data=direction, aes(ymin=lb, ymax=ub, x=position,colour=as.factor(tag))) +
scale_y_continuous(name="log-odds ratio (<7)") + 
scale_x_continuous(name=paste0("Position along chromosome ",chromosome)) +
geom_hline(yintercept=0, colour="red", linetype="dashed") +
theme(legend.position="none") +
scale_colour_manual(values = c(hues))

direct1<-ggplot(data=direction, aes(position, beta,colour=as.factor(tag))) + geom_point() +
geom_errorbar(data=direction, aes(ymin=lb, ymax=ub, x=position,colour=as.factor(tag)),size=0.1) +
scale_y_continuous(name="log-odds ratio (<7)") +
#scale_x_continuous(name=paste0("Position along chromosome ",chromosome)) +
geom_hline(yintercept=0, colour="red", linetype="dashed") +
scale_colour_manual(name="Tag group", values = c(hues))

direct3<-ggplot(data=direction, aes(position, beta,colour=as.factor(tag))) + geom_point() +
geom_errorbar(data=direction, aes(ymin=lb, ymax=ub, x=position,colour=as.factor(tag)),size=0.1) +
scale_y_continuous(name="log-odds ratio (<7)") +
#scale_x_continuous(name=paste0("Position along chromosome ",chromosome)) +
geom_hline(yintercept=0, colour="red", linetype="dashed") +
scale_colour_manual(name="Credible SNP group", values = c(hues))

trying2<-ggplot(data=e1, aes(SNPPos, logzscore, colour=tag1,shape=as.factor(GeneSymbol))) + geom_point() +
scale_y_continuous(name="Whole blood eQTL log z-score") +
scale_shape_discrete(name="Gene") +
#scale_x_continuous(name=paste0("Position along chromosome ",chromosome)) +
scale_colour_manual(name="Tag group",values = c(hues1)) +
geom_hline(yintercept=0, colour="red",linetype="dashed")

e2<-e1[order(e1$ppsum),]
e2$ppsum<-ifelse(is.na(e2$ppsum),0,e2$ppsum)
e2<-e2[!duplicated(e2$ppsum),]
e2<-e2[order(e2$ppsum),]
e2$ord<-c(1:nrow(e2))
e2<-e2[,c("tag1","ord")]
e1<-merge(e1,e2,by="tag1",all.x=T)
e1$tag1<-reorder(e1$tag1,e1$ord)
huesg1 = c("green","purple","turquoise1","orange","firebrick","tan","pink","lightblue")
allgenes<-c(names(huesg),unique(e1$GeneSymbol))
allgenes<-allgenes[!duplicated(allgenes)]
huesg1 =huesg1[1:length(allgenes)]
names(huesg1)<-allgenes

trying3<-ggplot(data=e1, aes(SNPPos, logzscore, alpha=as.factor(tag1),colour=as.factor(GeneSymbol))) + geom_point() +
scale_y_continuous(name="Whole blood eQTL log z-score") +
scale_alpha_discrete(name="Credible SNP group") +
#scale_x_continuous(name=paste0("Position along chromosome ",chromosome)) +
scale_colour_manual(name="Gene", values=huesg1) +
geom_hline(yintercept=0, colour="red",linetype="dashed")

#now trying a different graph without the wasted space.
lines<-ggplot(data=summx, aes(x=position, xend=position, y=0, yend=1,colour=as.factor(tag))) + geom_segment() +
scale_y_continuous(name="Credible SNPs") +
theme(axis.text.x=element_blank(),
axis.title.x=element_blank(),
axis.ticks.x=element_blank(),
axis.text.y=element_blank(),
axis.ticks.y=element_blank()) +
scale_colour_manual(name="Tag group",values=c(hues))

if(snp!="imm_9_4280823"){
out4<-tracks(t1,direct1,trying2, xlab=paste0("Position along chromosome ",chromosome),heights=c(1,1,3), xlim=c(min(o1$start,o1$end,min(e1$SNPPos)), max(e1$SNPPos)))
ggsave(out4,file=paste0(outdir,snp,"/new_eqtls_directions.png"),
dpi=800, height=20, width=15, units="cm")
}
if(snp=="imm_9_4280823"){
out4<-tracks(t1,direct1,trying2, xlab=paste0("Position along chromosome ",chromosome),heights=c(1,1,3), xlim=c(4150000, max(e1$SNPPos)))
ggsave(out4,file=paste0(outdir,snp,"/new_eqtls_directions.png"),
dpi=800, height=20, width=15, units="cm")
}

if(snp!="imm_9_4280823"){
out3<-tracks(direct, trying1,t,chr,signals,lds,heights=c(1,2,1,1,2,1),xlim=c(min(o1$start,o1$end), max(e1$SNPPos)))
ggsave(out3,file=paste0(outdir,snp,"/init_out_genes_eqtls_directions.png"),
dpi=200, height=40, width=20, units="cm")
}
if(snp=="imm_9_4280823"){
out3<-tracks(direct, trying1,t,chr,signals,lds,heights=c(1,2,1,1,2,1),xlim=c(4150000, max(e1$SNPPos)))
ggsave(out3,file=paste0(outdir,snp,"/init_out_genes_eqtls_directions.png"),
dpi=200, height=40, width=20, units="cm")
}

if(snp!="imm_9_4280823"){
out5<-tracks(t1,direct3,trying3, xlab=paste0("Position along chromosome ",chromosome),heights=c(1,1,3), xlim=c(min(o1$start,o1$end,min(e1$SNPPos)), max(e1$SNPPos)))
ggsave(out5,file=paste0(outdir,snp,"/new_eqtls_directions_tony.png"),
dpi=800, height=20, width=15, units="cm")
}
if(snp=="imm_9_4280823"){
out5<-tracks(t1,direct3,trying3, xlab=paste0("Position along chromosome ",chromosome),heights=c(1,1,3), xlim=c(4150000, max(e1$SNPPos)))
ggsave(out5,file=paste0(outdir,snp,"/new_eqtls_directions_tony.png"),
dpi=800, height=20, width=15, units="cm")
}
}
}

lapply(hits, readitallin)




