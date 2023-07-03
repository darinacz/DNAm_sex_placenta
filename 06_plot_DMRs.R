library(GWASTools)

#DMRs were calculated using comb-p 
#CpG had to be at least nominally significant to start a region

dmrs<-read.table("dmr_seed_nominal.regions-t.bed",sep="\t",header=F)
dim(dmrs) #n=66,581

#how many have at least 2 CpGs?
index<-which(dmrs$V5>1) #n=11,765
dmr<-dmrs[index,]

#for plotting, take the middle of the DMR
dmr$pos<-dmr$V2+(dmr$V3-dmr$V2)/2
dmr$CHR<-gsub("chr","",dmr$V1)
dmr$CHR<-as.numeric(dmr$CHR)

#how many regions are epi-genome wide significant?
length(which(dmr$V6<9e-08)) #n=3,120
index<-which(dmr$V6<9e-08)
top<-dmr[index,]
top<-top[,c(1,2,3,5,6)]
names(top)<-c("chr","start","stop","number_of_CpGs","p_region")
write.table(top,"top_DMRs.txt",sep="\t",quote=F,row.names=F)

#manhattan plot
dmr<-dmr[order(dmr$CHR,dmr$pos),]
manhattanPlot(p=dmr$V6,chr=dmr$CHR,signif=9e-08,ylim=c(0,100))
dev.print("manhattan_meta.png",dev=png,width=800)



