#hla_imputation.R
#imputing the HLAs of each individual in our analysis:
#need to use R 3.5 for this to work.
#library(snpStats,lib.loc="~/R/3.5.1/")
#library(GenomicRanges,lib.loc="~/R/3.5.1/")
#library(rtracklayer,lib.loc="~/R/3.5.1/")
#library(annotSnpStats, lib.loc="~/R/3.5.1/")

args = commandArgs(trailingOnly=TRUE)

#keep all individuals from multinomial analysis but don't drop any SNPs due to QC:

readit<-function(filename){
r<-read.plink(bed=paste0("/well/todd/users/jinshaw/t1d_risk/immunochip/",filename,".bed"),
bim=paste0("/well/todd/users/jinshaw/t1d_risk/immunochip/",filename,".bim"),
fam=paste0("/well/todd/users/jinshaw/t1d_risk/immunochip/",filename,".fam"))
r<-annot.plink(r)
return(r)
}
#genos<-lapply(c("sanger-controls-preqc", "t1d-cases-preqc", "uva-controls-preqc", "finn-preqc", "t1dgc-asp-preqc"),readit)
#rownames(genos[[2]])<-paste0(genos[[2]]@samples$pedigree,".",genos[[2]]@samples$member)
#rownames(genos[[5]])<-paste0(genos[[5]]@samples$pedigree,".",genos[[5]]@samples$member)

#geno<-rbind2(genos[[1]],genos[[2]])
#geno<-rbind2(geno,genos[[3]])
#genos[[4]]<-genos[[4]][,colnames(geno)]
#geno<-rbind2(geno,genos[[4]])
#geno<-rbind2(geno,genos[[5]])

#keep only individuals included in our analysis:
#inds<-read.table(file="/well/todd/users/jinshaw/aad/under_7/test_all_multinom_all.fam",header=F,as.is=T)

#try<-geno[rownames(geno) %in% inds$V2,]
#dim(try)
#try<-try[,colnames(try) %in% try@snps[try@snps$chromosome %in% c(1:22),"snp.name"]]
#write this out:
#write.plink(file.base="/well/todd/users/jinshaw/aad/under_7/nofilt",
#snps=as(try,"SnpMatrix"),
#pedigree=try@samples$pedigree,
#id=try@samples$member,
#mother=try@samples$mother,
#father=try@samples$father,
#sex=try@samples$sex,
#phenotype=try@samples$affected,
#chromosome=try@snps$chromosome,
#genetic.distance=try@snps$cM,
#position=try@snps$position,
#allele.1=try@snps$allele.1,
#allele.2=try@snps$allele.2)

#received this code from Dan.

dir="/well/todd/users/jinshaw/aad/under_7/"
fileStem="nofilt"

outFile="/well/todd/users/jinshaw/aad/under_7/imputation/HIBAG_new_l"

setwd(dir)

library(HIBAG,lib.loc="~/R/3.5.1/")

geno=hlaBED2Geno(bed.fn=paste(fileStem,".bed",sep=""),fam.fn=paste(fileStem,".fam",sep=""), bim.fn=paste(fileStem,".bim",sep=""), assembly="hg18")

fam=read.table(paste(fileStem,".fam",sep=""),header=F)

obj=hlaModelFiles("/gpfs2/well/todd/projects/HIBAG/ImmunoChip-European-HLA4-hg18.RData")

nLoci=length(obj)

i<-as.numeric(args[1])
     print(paste("Locus=",args[1]," of ",nLoci,sep=""))

    model=hlaModelFromObj(obj[[i]])

    predict=hlaPredict(model,geno, match.type="Position")

    res=predict[[4]]

write.table(res,file=paste0(outFile,i),quote=F,row.names=T,col.names=T,append=T)





