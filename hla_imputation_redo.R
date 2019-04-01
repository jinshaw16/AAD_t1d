#hla_imputation_redo.R
#imputing the HLAs of each individual in our analysis:
#need to use R 3.5 for this to work.
args = commandArgs(trailingOnly=TRUE)


#received this code from Dan.
dir="/well/todd/users/jinshaw/aad/"
fileStem="all_ages_n"

outFile="/well/todd/users/jinshaw/aad/under_7/imputation/HIBAG_new_redo_2"

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

write.table(res,file=paste0(outFile,i),quote=F,row.names=T,col.names=T)





