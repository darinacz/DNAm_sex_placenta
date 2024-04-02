library(data.table)
library(GWASTools)
library(ggplot2)

#read in results from meta-analysis
meta<-fread("meta_BET_PREDO_ITU_random1.tbl",header=T,sep="\t")
#ARE: additive random effects, this is based on the DerSimonian & Laird algorith

length(which(meta$PvalueARE<9e-08)) #n=10,320 epigenome-wide hits

#load annotation, this is the annotation of the EPIC array as provided by the minfi package
load("annot_epic.Rdata")
annot<-annot_epic[,c(4,2,1,22,24)] #subset to columns of interest
res<-merge(meta,annot,by.x="MarkerName",by.y="Name")
res<-res[,-c(2,3)]

res$Chr<-substr(res$chr,start=4,stop=nchar(res$chr))
res$Chr<-as.numeric(res$Chr)
res<-res[order(res$Chr,res$pos),] #order by chromosome and position

write.table(res,"meta_results.txt",sep="\t",quote=F,row.names=F)
index<-which(res$PvalueARE<9e-08)
write.table(res[index,],"epigenome_wide_significant_hits.txt",sep="\t",quote=F,row.names=F)


#####Figure 1A####
#Manhattan plot
pdf("manhattan_meta.pdf")
manhattanPlot(p=res$PvalueARE,chr=res$Chr,signif=9e-08,ylim=c(0,250))
dev.off()

#plot top hit
#cg12691488, p_meta=2.9e-248
load("../../BET/rlm/BetaSet_BET_137.rda")
beta_bet<-exprs(Beta_BET_137)
pt_bet<-pData(Beta_BET_137)
load("../../PREDO/rlm/BetaSet_PREDO_139.rda")
beta_predo<-exprs(Beta_PREDO_139)
pt_predo<-pData(Beta_PREDO_139)
load("../../ITU/rlm/BetaSet_ITU_470.rda")
beta_itu<-exprs(Beta_ITU_470)
pt_itu<-pData(Beta_ITU_470)

bet<-rep("BET",137)
itu<-rep("ITU",470)
predo<-rep("PREDO",139)
test<-rbind(
  cbind(B2M(beta_bet["cg12691488",]),pt_bet$sex,bet),
  cbind(B2M(beta_itu["cg12691488",]),pt_itu$sex,itu),
  cbind(B2M(beta_predo["cg12691488",]),pt_predo$sex,predo))
pdf("cg12691488_meta.pdf")
boxplot(as.numeric(test[,1]) ~test[,2]*test[,3],col=c("blue","red"),ylab="M-value cg12691488",xlab="",cex.axis=1.4,cex.lab=1.4)
stripchart(as.numeric(test[,1])~test[,2]*test[,3],vertical = TRUE, method = "jitter", add = TRUE, pch = 16, cex=0.5, col = 'grey')
dev.off() 

#######Table S1########
#combine all results in one table
#need to merge res with individual results
bet<-read.table("../../BET/rlm/results_BET_sex.txt",header=T,sep="\t")
names(bet)[c(2:4)]<-c("beta_BET","se_BET","p_BET")
bet<-bet[,-c(5,6,7)]
all<-merge(res,bet,by.x="MarkerName",by.y="CpG",all.x=T,all.y=T)
itu<-read.table("../../ITU/rlm/results_ITU_sex.txt",header=T,sep="\t")
names(itu)[c(2:4)]<-c("beta_ITU","se_ITU","p_ITU")
itu<-itu[,-c(5,6,7)]
all<-merge(all,itu,by.x="MarkerName",by.y="CpG",all.x=T,all.y=T)
predo<-read.table("../../PREDO/rlm/results_PREDO_sex.txt",header=T,sep="\t")
names(predo)[c(2:4)]<-c("beta_PREDO","se_PREDO","p_PREDO")
predo<-predo[,-c(5,6,7)]
all<-merge(all,predo,by.x="MarkerName",by.y="CpG",all.x=T,all.y=T)
all<-all[,c(1,21,16,5,10,11,12,18,19,22:30)]
all<-all[,c(1,2,3,5,6,7,4,10:18,8,9)]
names(all)[c(1:7,17,18)]<-c("CpG","chr","pos_hg19","beta_meta","se_meta","p_meta","direction","gene","relation to gene")
#only keep first gene and first relation
all$gene <- sub(';.*$','',all$gene)
all$`relation to gene` <- sub(';.*$','',all$`relation to gene`)
names(all)[18]<-"position_in_gene"
meta_with_pos<-all
meta_with_pos<-meta_with_pos[order(meta_with_pos$chr,meta_with_pos$pos_hg19),]
write.table(meta_with_pos,"meta_final_with_pos.txt",sep="\t",quote=F,row.names=F)
index<-which(meta_with_pos$p_meta<9e-08) #n=10,320
write.table(meta_with_pos[index,],"meta_final_with_pos_epigenomewide.txt",sep="\t",quote=F,row.names=F) #This is Table S1

