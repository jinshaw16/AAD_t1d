#hla_readin.R

#this reads in the HLA imputation result and tidies them up for analysis:
library(ggplot2)
library(epicalc)
library(gridExtra)
library(multinomRob)
d<-"/well/todd/users/jinshaw/t1d_risk/immunochip/"
aaddir<-"/well/todd/users/jinshaw/aad/"


#read imputation results into R:
read<-function(num,gene){
A<-read.table(file=paste0("/well/todd/users/jinshaw/aad/under_7/imputation/HIBAG_new_redo_3",num),header=T, as.is=T)
colnames(A)[2:5]<-c(paste0(gene,"1"),paste0(gene,"2"),paste0(gene,"prob"),paste0(gene,"matching"))
return(A)
}
j<-mapply(read, num=c(1:7), gene=c("A","B","C","DRB1","DQA1","DQB1","DPB1"), SIMPLIFY=F)
j<-do.call("cbind",j)
j<-j[,c(-6,-11,-16,-21,-26, -31)]


#get pheno and covar info to merge in:
load(file=paste0(d,"all_inds_unrel_postqc_3.RData"))
pheno<-all@samples

#names have changed over time for initial AAD publication so just re-mapping
j$sample.id<-gsub("sanger-controls-","sanger-controls.",j$sample.id)
j$sample.id<-gsub("t1d-cases-","t1d-cases.",j$sample.id)
j$sample.id<-gsub("uva-controls-","uva-controls.",j$sample.id)
j$sample.id<-gsub("ni-","ni.",j$sample.id)
j$sample.id<-gsub("finn-","finn.",j$sample.id)
j$sample.id<-ifelse(substr(j$sample.id,1,1) %in% c("1","2","4","5"),
gsub("-",".",j$sample.id),j$sample.id)
#j$fins1<-ifelse(substr(j$sample.id,1,3) %in% c("FIN"),sub("(^FIN-.*)(.*)(-)(.*)","\\1",j$sample.id),j$sample.id)
#j$fins2<-ifelse(substr(j$sample.id,1,3) %in% c("FIN"),sub("(^FIN-.*)(.*)(-)(.*)","\\2",j$sample.id),j$sample.id)
#j$fins3<-ifelse(substr(j$sample.id,1,3) %in% c("FIN"),sub("(^FIN-.*)(.*)(-)(.*)","\\3",j$sample.id),j$sample.id)
#j$fins4<-ifelse(substr(j$sample.id,1,3) %in% c("FIN"),sub("(^FIN-.*)(.*)(-)(.*)","\\4",j$sample.id),j$sample.id)
#j$sample.id<-ifelse(substr(j$sample.id,1,3)=="FIN",paste0(j$fins1,j$fins2,".",j$fins4),j$sample.id)

j<-j[j$sample.id %in% pheno$uniqueID,]
pheno<-pheno[j$sample.id,]

hla<-cbind(j, pheno)
hla$onset<-as.numeric(hla$onset)
hla$group<-ifelse(hla$affected==1,0,
ifelse(hla$onset<7 & hla$affected==2,1,
ifelse(hla$onset>=7 & hla$onset<13 & hla$affected==2,2,
ifelse(hla$onset>=13 & hla$affected==2,3,NA))))

#make those calls with <50% certainty be missing:
for( i in c("DRB1","DQA1","DQB1","A","B","DPB1")){
hla[,paste0(i,"1")]<-ifelse(hla[,paste0(i,"prob")]<0.50,NA,hla[,paste0(i,"1")])
hla[,paste0(i,"2")]<-ifelse(hla[,paste0(i,"prob")]<0.50,NA,hla[,paste0(i,"2")])
}

#########################################################################################################################
#define highest and lowest risk class 2 using DRB1 and DQB1 only (DQA1 less well imputed and included mainly due to LD):#
#########################################################################################################################

hla$risk1<-ifelse(is.na(hla$DRB11) | is.na(hla$DQB11),NA,"X")
hla$risk2<-ifelse(is.na(hla$DRB12) | is.na(hla$DQB12),NA,"X")

definehlarisk<-function(var,dr,dq, label){
hla[,paste0("risk",var)]<-ifelse(hla[,paste0("DRB1",var)]==dr & hla[,paste0("DQB1",var)]==dq,label,hla[,paste0("risk",var)])
hla <<- hla
}
invisible(mapply(definehlarisk, var=c(rep("1",21),rep("2",21)),
dr=c("03:01",
"04:01","04:02","04:04","04:05","04:08",
"04:01","04:02","04:04","04:05","04:08",
"04:01","04:02","04:04","04:05","04:08",
"13:03","11:04","15:01","07:01","14:01",
"03:01",
"04:01","04:02","04:04","04:05","04:08",
"04:01","04:02","04:04","04:05","04:08",
"04:01","04:02","04:04","04:05","04:08",
"13:03","11:04","15:01","07:01","14:01"),
dq=c("02:01",
"03:02","03:02","03:02","03:02","03:02",
"03:04","03:04","03:04","03:04","03:04",
"02:02","02:02", "02:02","02:02","02:02",
"03:01","03:01","06:02","03:03","05:03",
"02:01",
"03:02","03:02","03:02","03:02","03:02",
"03:04","03:04","03:04","03:04","03:04",
"02:02","02:02", "02:02","02:02","02:02",
"03:01","03:01","06:02","03:03","05:03"),
label=c("DR3",
"DR4","DR4","DR4","DR4","DR4",
"DR4","DR4","DR4","DR4","DR4",
"DR4","DR4","DR4","DR4","DR4",
"DR13","DR11","DR15","DR7","DR14",
"DR3",
"DR4","DR4","DR4","DR4","DR4",
"DR4","DR4","DR4","DR4","DR4",
"DR4","DR4","DR4","DR4","DR4",
"DR13","DR11","DR15","DR7","DR14")))


hla$t1d<-ifelse(hla$affected==1,0,ifelse(hla$affected==2,1,NA))

#define their combination:
hla$dr<-ifelse(is.na(hla$risk1) | is.na(hla$risk2),NA,"X/X")

defdr<-function(r1,r2, label){
hla$dr<-ifelse(hla$risk1==r1 & hla$risk2==r2 |
hla$risk1==r2 & hla$risk2==r1,label,hla$dr)
hla<<-hla
}
invisible(mapply(defdr,r1=c("DR3","DR3","DR4","DR13","DR11","DR15","DR7","DR14","DR3","DR4","DR13","DR11","DR15","DR7","DR14",
                  "DR3","DR14","DR3","DR3","DR4","DR7","DR3","DR4","DR11","DR7","DR11","DR3","DR4","DR13","DR11","DR13"),
             r2=c("DR4","X","X","X","X","X","X","X","DR3","DR4","DR13","DR11","DR15","DR7","DR14",
                  "DR7","DR15","DR14","DR15","DR15","DR14","DR11","DR7","DR15","DR15","DR13","DR13","DR14","DR15","DR14","DR14"),
          label=c("DR3/4","DR3/X","DR4/X","DR13/X","DR11/X","DR15/X","DR7/X","DR14/X","DR3/3","DR4/4",
                   "DR13/13","DR11/11","DR15/15","DR7/7","DR14/14","DR3/7","DR14/15","DR3/14","DR3/15",
                   "DR4/15","DR7/14","DR3/11","DR4/7","DR11/15","DR7/15","DR11/13","DR3/13","DR4/14","DR13/15","DR11/14","DR13/14")))


#redefine some of the smaller categories:
hla$dr<-ifelse(hla$dr=="DR11/11" |hla$dr=="DR11/13" |hla$dr=="DR11/15" |
hla$dr=="DR11/13" |hla$dr=="DR13/13" |hla$dr=="DR13/14" |hla$dr=="DR13/15" |
hla$dr=="DR14/14" |hla$dr=="DR11/13" |hla$dr=="DR11/14" |hla$dr=="DR7/7" |
hla$dr=="DR14/15"|hla$dr=="DR7/14"|hla$dr=="DR7/15","rareprot",hla$dr)

hla$dr<-ifelse(hla$dr=="DR3/7" |hla$dr=="DR3/13" |hla$dr=="DR3/15" |
hla$dr=="DR3/11" |hla$dr=="DR3/14"|hla$dr=="DR4/7" |hla$dr=="DR4/13" |hla$dr=="DR4/15" |
hla$dr=="DR4/11"|hla$dr=="DR4/14","Protective/Susceptible",hla$dr)

hla$ord<-ifelse(hla$dr=="X/X",0,
ifelse(hla$dr=="DR3/X",1,ifelse(hla$dr=="DR3/3",2,ifelse(hla$dr=="DR4/X",3,
ifelse(hla$dr=="DR4/4",4,ifelse(hla$dr=="DR3/4",5,ifelse(hla$dr=="DR13/X",7,
ifelse(hla$dr=="DR11/X",8,ifelse(hla$dr=="DR15/X",9,ifelse(hla$dr=="DR7/X",10,
ifelse(hla$dr=="DR14/X",11,ifelse(hla$dr=="DR15/15",12,ifelse(hla$dr=="protsus",13,
ifelse(hla$dr=="rareprot",14,NA))))))))))))))

hla$dr<-as.factor(hla$dr)
hla$dr<-reorder(hla$dr,hla$ord)


summary(glm(data=hla, t1d ~ as.factor(dr) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 , family="binomial"))
und<-hla[hla$group==0 | hla$group==1,]
midrange<-hla[hla$group==0 | hla$group==2,]
ov<-hla[hla$group==0 | hla$group==3,]
getresforgp<-function(frame,name){
underc2<-as.data.frame(summary(glm(data=frame, t1d ~ as.factor(dr) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 , family="binomial"))$coefficients)

underc2$hap<-rownames(underc2)
underc2$hap<-gsub(".*DR","",underc2$hap)
underc2<-underc2[2:14,]
underc2$type=name
underc2$ord<-c(nrow(underc2):1)
return(underc2)
}

both<-mapply(getresforgp,frame=list(und,midrange,ov),name=c("under","mid","over"),SIMPLIFY=FALSE)
both<-do.call("rbind",both)
both$lb<-both$Estimate-(qnorm(0.975)*both$`Std. Error`)
both$ub<-both$Estimate+(qnorm(0.975)*both$`Std. Error`)
both$hap<-as.factor(both$hap)
both$hap<-reorder(both$hap, both$ord)
g<-ggplot(data=both[both$type=="under",], aes(y=as.factor(hap),x=Estimate)) + geom_point(,colour="red") +
geom_errorbarh(data=both[both$type=="under",],aes(xmin=lb, xmax=ub, y=as.factor(hap)),colour="red",height=0) +
geom_point(data=both[both$type=="over",], aes(y=(as.numeric(hap)+0.4),x=Estimate),colour="blue") +
geom_errorbarh(data=both[both$type=="over",],aes(xmin=lb, xmax=ub, y=(as.numeric(hap)+0.4)),colour="blue",height=0) +
geom_point(data=both[both$type=="mid",], aes(y=(as.numeric(hap)+0.2),x=Estimate),colour="green") +
geom_errorbarh(data=both[both$type=="mid",],aes(xmin=lb, xmax=ub, y=(as.numeric(hap)+0.2)),colour="green",height=0) +
coord_cartesian(xlim=c(-5,5)) +
geom_vline(xintercept=0,colour="red",linetype="dashed") +
scale_y_discrete(name="Classical HLA allele (vs. DRX/X)") +
scale_x_continuous(name="T1D log-odds ratio")

#defining variables to use later in multinomial regression analyses:

definethem<-function(varname,hap){
allhap<-paste0(hap,"/",substr(hap,3,4))

#first under definition 'haplotpye vs. every other non interesting haplotype':
hla[,paste0(varname)]<-ifelse((hla$risk1==hap & hla$risk2!=hap)|
(hla$risk2==hap & hla$risk1!=hap),1,ifelse(hla$dr==allhap,2,ifelse(hla$dr=="X/X",0,NA)))

#second under definition 'haplotype vs. every other haplotype (including eg dr3 if examining dr4):
hla[,paste0(varname,"_1")]<-ifelse((hla$risk1==hap & hla$risk2!=hap)|
(hla$risk2==hap & hla$risk1!=hap),1,ifelse(hla$dr==allhap,2,ifelse(is.na(hla$dr),NA,0)))
h1<-hla[,c(paste0(varname),paste0(varname,"_1"))]
return(h1)
}
recoded<-mapply(definethem, hap=c("DR3","DR4","DR11","DR13","DR7","DR14","DR15"),
varname=c("dr3","dr4","dr11","dr13","dr7","dr14","dr15"), SIMPLIFY=FALSE)
names(recoded)<-NULL
recoded<-do.call("cbind",recoded)
hla<-cbind(hla, recoded)
#and define DR3/4:
hla$dr34<-ifelse(hla$risk1=="DR3" & hla$risk2=="DR4" |
hla$risk1=="DR4" & hla$risk2=="DR3", 1,
ifelse(hla$dr=="X/X",0,NA))
hla$dr34_1<-ifelse(hla$risk1=="DR3" & hla$risk2=="DR4" |
hla$risk1=="DR4" & hla$risk2=="DR3", 1,
ifelse(is.na(hla$dr),NA,0))

#CLASS 1:
#Seen evidence in literature of T1D association with DPB1*02:02 and DPB1*03:01 and DP*04:02:
hla$DPrisk<-ifelse(is.na(hla$DPB11) | is.na(hla$DPB12),NA,"X")

hla$DPB11_d<-ifelse(!hla$DPB11 %in% c("02:02","03:01","04:02") & !is.na(hla$DPB11),"X",ifelse(is.na(hla$DPB11),NA,hla$DPB11))
hla$DPB12_d<-ifelse(!hla$DPB12 %in% c("02:02","03:01","04:02") & !is.na(hla$DPB11),"X",ifelse(is.na(hla$DPB12),NA,hla$DPB12))
redefclass1<-function(gene,h1,h2,label,varname){
hla[,varname]<-ifelse((hla[,paste0(gene,"1_d")]==h1 & hla[,paste0(gene,"2_d")]==h2) |
(hla[,paste0(gene,"1_d")]==h2 & hla[,paste0(gene,"2_d")]==h1),label, hla[,varname])
hla<<-hla
}
invisible(mapply(redefclass1, gene=c(rep("DPB1",9)),
varname=c(rep("DPrisk",9)),
h1=c("02:02","02:02","03:01","03:01","04:02","04:02","02:02","02:02","03:01"),
h2=c("X","02:02","X","03:01","X","04:02","03:01","04:02","04:02"),
label=c("DPB1*02:02/X","DPB1*02:02/02:02","DPB1*03:01/X","DPB1*03:01/03:01","DPB1*04:02/X",
"DPB1*04:02/04:02","DPB1*02:02/03:01","DPB1*02:02/04:02","DPB1*03:01/04:02")))


hla$DPrisk<-ifelse(hla$DPrisk=="DPB1*02:02/02:02" | hla$DPrisk=="DPB1*02:02/X" |
hla$DPrisk=="DPB1*02:02/04:02"| hla$DPrisk=="DPB1*02:02/03:01", "DP*02:02",
ifelse(hla$DPrisk=="DPB1*03:01/04:02" | hla$DPrisk=="DPB1*02:02/04:02", "Susceptible/protective",hla$DPrisk))

hla$ordDP<-ifelse(hla$DPrisk=="X",0,
ifelse(hla$DPrisk=="DP*02:02",1,
ifelse(hla$DPrisk=="DPB1*03:01/X",2,
ifelse(hla$DPrisk=="DPB1*03:01/03:01",3,
ifelse(hla$DPrisk=="DPB1*04:02/X",4,
ifelse(hla$DPrisk=="DPB1*04:02/04:02",5,
ifelse(hla$DPrisk=="Susceptible/protective",6,NA)))))))

hla$DPrisk<-as.factor(hla$DPrisk)
hla$DPrisk<-reorder(hla$DPrisk,hla$ordDP)

#There are a number of potentially disease-associated A allele haplotypes:
#02:01, 11:01, 24:02, 66:01 and 32:01
hla$A1_d<-ifelse(!hla$A1 %in% c("02:01","11:01","24:02","66:01","32:01") & !is.na(hla$A1),"X", ifelse(is.na(hla$A1),NA,hla$A1))
hla$A2_d<-ifelse(!hla$A2 %in% c("02:01","11:01","24:02","66:01","32:01") & !is.na(hla$A2),"X", ifelse(is.na(hla$A2),NA,hla$A2))
hla$Arisk<-ifelse(is.na(hla$A1) | is.na(hla$A2),NA,"X")
invisible(mapply(redefclass1, gene=c(rep("A",20)),
varname=c(rep("Arisk",20)),
h1=c("02:01","02:01","11:01","11:01","24:02","24:02","66:01","66:01","32:01","32:01","02:01",
"02:01","02:01","02:01","24:02","11:01","11:01","11:01","32:01","24:02"),
h2=c("02:01","X","11:01","X","24:02","X","66:01","X","32:01","X","24:02","11:01","32:01","66:01",
"32:01","32:01","66:01","24:02","66:01","66:01"),
label=c("A*02:01/02:01","A*02:01/X","A*11:01/11:01","A*11:01/X","A*24:02/24:02","A*24:02/X",
"A*66:01/66:01","A*66:01/X","A*32:01/32:01","A*32:01/X","A*02:01/24:02",
"A*02:01/11:01", "A*02:01/32:01", "A*02:01/66:01","A*24:02/32:01",
"A*11:01/32:01","A*11:01/66:01","A*11:01/24:02","A*32:01/66:01",
"A*24:02/66:01")))


#combine some rare categories:
hla$Arisk<-ifelse(hla$Arisk=="A1*32:01/66:01" | hla$Arisk=="A1*11:01/66:01" | hla$Arisk=="A1*11:01/32:01", "Rareprotective",
ifelse(hla$Arisk=="A1*02:01/11:01"|hla$Arisk=="A1*02:01/66:01" | hla$Arisk=="A1*02:01/32:01" | hla$Arisk=="A1*11:01/24:02"|
hla$Arisk=="A1*24:01/66:01"|hla$Arisk=="A1*24:02/32:01","Susceptible/protective",hla$Arisk))

hla$ordA<-ifelse(hla$Arisk=="X",0,
ifelse(hla$Arisk=="A*02:01/X",1,
ifelse(hla$Arisk=="A*02:01/02:01",2,
ifelse(hla$Arisk=="A*24:02/X",3,
ifelse(hla$Arisk=="A*24:02/24:02",4,
ifelse(hla$Arisk=="A*02:01/24:02",5,
ifelse(hla$Arisk=="A*11:01/X",6,
ifelse(hla$Arisk=="A*11:01/11:01",7,
ifelse(hla$Arisk=="A*32:01/X",8,
ifelse(hla$Arisk=="A*32:01/32:01",9,
ifelse(hla$Arisk=="A*66:01/X",9,
ifelse(hla$Arisk=="A*66:01/66:01",10,
ifelse(hla$Arisk=="Susceptible/protective",11,
ifelse(hla$Arisk=="Rareprotective",12,NA))))))))))))))

hla$Arisk<-as.factor(hla$Arisk)
hla$Arisk<-reorder(hla$Arisk,hla$ordA)



#There are a number of potentially disease-associated B allele haplotypes:
#07:02, 18:01, 35:02, 39:06 and 44:03
hla$B1_d<-ifelse(!hla$B1 %in% c("07:02","18:01","35:02","39:06","44:03") & !is.na(hla$B1),"X",ifelse(is.na(hla$B1),NA,hla$B1))
hla$B2_d<-ifelse(!hla$B2 %in% c("07:02","18:01","35:02","39:06","44:03") & !is.na(hla$B2),"X",ifelse(is.na(hla$B2),NA,hla$B2))
hla$Brisk<-ifelse(is.na(hla$B1) | is.na(hla$B2),NA,"X")
invisible(mapply(redefclass1, gene=c(rep("B",17)),
varname=c(rep("Brisk",17)),
h1=c("07:02","07:02","18:01","18:01","35:02","35:02","39:06","39:06","44:03","44:03","39:06","18:01","18:01",
"35:02","18:01","35:02","35:02"),
h2=c("07:02","X","18:01","X","35:02","X","39:06","X","44:03","X","44:03","44:03","39:06","44:03","35:02","44:03","39:06"),
label=c("B*07:02/07:02","B*07:02/X","B*18:01/18:01","B*18:01/X",
"B*35:02/35:02","B*35:02/X","B*39:06/39:06","B*39:06/X","B*44:03/44:03","B*44:03/X",
"B*39:06/44:03","B*18:01/44:03","B*18:01/39:06","B*35:02/44:03","B*18:01/35:02",
"B*35:02/44:03","B*35:02/39:06")))


hla$Brisk<-ifelse(hla$Brisk=="B*35:02/39:06"|hla$Brisk=="B*18:01/35:02"|hla$Brisk=="B*39:06/44:03","Protective/susceptible",
ifelse(hla$Brisk=="B*35:02/44:03" | hla$Brisk=="B*35:02/35:02"|hla$Brisk=="B*44:03/44:03", "Rareprotective",
ifelse(hla$Brisk=="B*39:06/39:06"|hla$Brisk=="B*18:01/39:06"|hla$Brisk=="B*18:01/44:03"|hla$Brisk=="B*18:01/18:01","Raresusceptible",hla$Brisk)))

hla$ordB<-ifelse(hla$Brisk=="X",0,
ifelse(hla$Brisk=="B*18:01/X",1,
ifelse(hla$Brisk=="B*39:06/X",2,
ifelse(hla$Brisk=="B*07:02/X",3,
ifelse(hla$Brisk=="B*07:02/07:02",4,
ifelse(hla$Brisk=="B*35:02/X",5,
ifelse(hla$Brisk=="B*44:03/X",6,
ifelse(hla$Brisk=="Rareprotective",7,
ifelse(hla$Brisk=="Raresusceptible",8,
ifelse(hla$Brisk=="Protective/susceptible",9,NA))))))))))
hla$Brisk<-as.factor(hla$Brisk)
hla$Brisk<-reorder(hla$Brisk, hla$ordB)


#and one C disease-associated allele:C*03:03
hla$C1_d<-ifelse(!hla$C1 %in% c("03:03") & !is.na(hla$C1),"X",ifelse(is.na(hla$C1),NA,hla$C1))
hla$C2_d<-ifelse(!hla$C2 %in% c("03:03") & !is.na(hla$C2),"X",ifelse(is.na(hla$C2),NA,hla$C2))
hla$Crisk<-ifelse(is.na(hla$C1) | is.na(hla$C2),NA,"X")
invisible(mapply(redefclass1, gene=c(rep("C",2)),
varname=c(rep("Crisk",2)),
h1=c("03:03","03:03"),
h2=c("03:03","X"),
label=c("C*03:03/03:03","C*03:03/X")))


hla$ordC<-ifelse(hla$Crisk=="X",0,
ifelse(hla$Crisk=="C*03:03/X",1,
ifelse(hla$Crisk=="C*03:03/03/03",2,NA)))
hla$Crisk<-as.factor(hla$Crisk)
hla$Crisk<-reorder(hla$Crisk, hla$ordC)


hla$dpb10301<-ifelse(hla$DPrisk=="DPB1*03:01/X",1, ifelse(hla$DPrisk=="DPB1*03:01/03:01",2, ifelse(is.na(hla$DPrisk),NA,0)))
hla$dpb10402<-ifelse(hla$DPrisk=="DPB1*04:02/X",1, ifelse(hla$DPrisk=="DPB1*04:02/04:02",2, ifelse(is.na(hla$DPrisk),NA,0)))
hla$a0201<-ifelse(hla$Arisk=="A*02:01/X",1, ifelse(hla$Arisk=="A*02:01/02:01",2, ifelse(is.na(hla$Arisk),NA,0)))
hla$a2402<-ifelse(hla$Arisk=="A*24:02/X",1, ifelse(hla$Arisk=="A*24:02/24:02",2, ifelse(is.na(hla$Arisk),NA,0)))
hla$a0224<-ifelse(hla$Arisk=="A*02:01/24:02",1,ifelse(is.na(hla$Arisk),NA,0))
hla$a1101<-ifelse(hla$Arisk=="A*11:01/X",1, ifelse(hla$Arisk=="A*11:01/11:01",2, ifelse(is.na(hla$Arisk),NA,0)))
hla$a3201<-ifelse(hla$Arisk=="A*32:01/X",1, ifelse(hla$Arisk=="A*32:01/32:01",2, ifelse(is.na(hla$Arisk),NA,0)))
hla$b1801<-ifelse(hla$Brisk=="B*18:01/X",1,ifelse(is.na(hla$Brisk),NA,0))
hla$b3906<-ifelse(hla$Brisk=="B*39:06/X",1,ifelse(is.na(hla$Brisk),NA,0))
hla$b4403<-ifelse(hla$Brisk=="B*44:03/X",1,ifelse(is.na(hla$Brisk),NA,0))

#and defining in the other way (additive):
hla$dpb10301_1<-hla$dpb10301
hla$dpb10402_1<-hla$dpb10402
hla$a0201_1<-ifelse(hla$A1_d!="02:01" & hla$A2_d!="02:01",0,
ifelse((hla$A1_d=="02:01" & hla$A2_d!="02:01") | (hla$A1_d!="02:01" & hla$A2_d=="02:01"),1,
ifelse(hla$A1_d=="02:01" & hla$A2_d=="02:01",2,NA)))
hla$a2402_1<-ifelse(hla$A1_d!="24:02" & hla$A2_d!="24:02",0,
ifelse((hla$A1_d=="24:02" & hla$A2_d!="24:02") | (hla$A1_d!="24:02" & hla$A2_d=="24:02"),1,
ifelse(hla$A1_d=="24:02" & hla$A2_d=="24:02",2,NA)))
hla$a1101_1<-hla$a1101
hla$a3201_1<-hla$a3201
hla$b1801_1<-hla$b1801
hla$b3906_1<-hla$b3906
hla$b4403_1<-hla$b4403

save(hla, file="/well/todd/users/jinshaw/aad/under_7/hla_all_3.RData")
