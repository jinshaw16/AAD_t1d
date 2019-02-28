#multinomial_split_iter_readin.R

#reading in the results from multinomial_split_iter.R
library(ggplot2)

getres<-function(args){
j<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_iter_",args,"_1_n.txt"),
header=T,as.is=T,sep="\t")
j$iter=args
j$locus=c("PTPN22", "TNFSF4","CAMSAP2",
"IL10", "AFF3", "ACOXL",
"IFIH1 (1)", "IFIH1 (2)", "IFIH1 (3)",
"CTLA4", "CCR5", "IL2/IL21 (1)", "IL2/IL21 (2)", "CPE",
"IL7R","PTPRK/THEMIS", "BACH2", "CENPW",
"IKZF1", "COBL",
"GLIS3", "IL2RA (1)", "IL2RA (2)", "IL2RA (3)","IL2RA (4)",
"PTEN","INS (1)","INS (2)", "CD69", "IKZF4",
"SH2B3", "GPR183", "LINC01550",
"MEG3","RASGRP1", "CTSH", "IL27", "DEXI (1)",
"DEXI (2)", "CTRB1","IKZF3",
"CCR7", "MAPT","PTPN2 (1)", "PTPN2 (2)", "CD226", "TYK2 (1)", "TYK2 (2)", "PRKD2",
"FUT2", "SIRPG",
"UBASH3A", "ICOSLG", "HORMAD2",
"C1QTNF6")
j$half=1
j1<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/results/inds_likelihoods_iter_",args,"_2_n.txt"),
header=T,as.is=T,sep="\t")
j1$iter=args
j1$locus=c("PTPN22", "TNFSF4","CAMSAP2",
"IL10", "AFF3", "ACOXL",
"IFIH1 (1)", "IFIH1 (2)", "IFIH1 (3)",
"CTLA4", "CCR5", "IL2/IL21 (1)", "IL2/IL21 (2)", "CPE",
"IL7R","PTPRK/THEMIS", "BACH2", "CENPW",
"IKZF1", "COBL",
"GLIS3", "IL2RA (1)", "IL2RA (2)", "IL2RA (3)","IL2RA (4)",
"PTEN","INS (1)","INS (2)", "CD69", "IKZF4",
"SH2B3", "GPR183", "LINC01550",
"MEG3","RASGRP1", "CTSH", "IL27", "DEXI (1)",
"DEXI (2)", "CTRB1","IKZF3",
"CCR7", "MAPT","PTPN2 (1)", "PTPN2 (2)", "CD226", "TYK2 (1)", "TYK2 (2)", "PRKD2",
"FUT2", "SIRPG",
"UBASH3A", "ICOSLG", "HORMAD2",
"C1QTNF6")
j1$half=2
j2<-rbind(j,j1)
return(j2)
}
g<-lapply(c(1:100),getres)
g<-do.call("rbind",g)
g$under<-g$p<0.05

g$prop<-NA
getprop<-function(locus){
p<-g[g$locus==locus,]
t<-table(p$p<0.05)
if(length(t)==2){
prop<-t[2]/sum(t)
}
if(length(t)==1 & names(t)=="TRUE"){
prop=1
}
if(length(t)==1	& names(t)=="FALSE"){
prop=0
}
g$prop<<-ifelse(g$locus==locus,prop,g$prop)
}
invisible(lapply(unique(g$locus),getprop))
test<-g[1:55,]
test<-test[order(-test$prop),]
test$order<-1:nrow(test)
test<-test[,c("locus","order")]
g<-merge(g,test,by="locus")
g$locus<-as.factor(g$locus)
g$locus<-reorder(g$locus, g$order)
png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/split_props_all_n.png",
width=25,height=20,units="cm",res=400)
ggplot(g,aes(as.factor(locus), fill=under)) + geom_bar(position="fill") +
theme(axis.text.x=element_text(angle=90)) +
scale_fill_hue(name="Nominal significance") +
scale_y_continuous(name="Proportion of samples") +
scale_x_discrete(name="Locus")
dev.off()

g1<-g[as.numeric(g$order)<25,]
png(file="/well/todd/users/jinshaw/output/aad/under_7/multinom/redo/split_props_top_n.png",
width=25,height=20,units="cm",res=400)
ggplot(g1,aes(as.factor(locus), fill=under)) + geom_bar(position="fill")	+
theme(axis.text.x=element_text(angle=90)) +
scale_fill_hue(name="Nominal significance") +
scale_y_continuous(name="Proportion of samples") +
scale_x_discrete(name="Locus")
dev.off()
