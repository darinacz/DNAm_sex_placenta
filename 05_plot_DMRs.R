library(GWASTools)

#DMRs were calculated using comb-p 
#CpG had to be at least nominally significant to start a region

dmrs<-read.table("dmr_seed_nominal.regions-t.bed",sep="\t",header=F)
dim(dmrs) #n=23,650

index<-which(dmrs$n_probes>=2) #n=4,046
dmr<-dmrs[index,]
mean(dmr$n_probes) #2.64 CpGs
mean(dmr$end-dmr$start) #91.33 bps
#for plotting, take the middle of the DMR
dmr$pos<-dmr$start+(dmr$end-dmr$start)/2
dmr$CHR<-gsub("chr","",dmr$chrom)
dmr$CHR<-as.numeric(dmr$CHR)

#how many regions are epi-genome wide significant?
length(which(dmr$z_p<9e-08)) #n=1,199
index<-which(dmr$z_p<9e-08)

########Table S2###########
top<-dmr[index,]
top<-top[,c(1,2,3,5,6)]
names(top)<-c("chr","start","stop","number_of_CpGs","p_region")
write.table(top,"top_DMRs.txt",sep="\t",quote=F,row.names=F) #Table S2
summary(top$number_of_CpGs)


#########Figure S2###########
#manhattan plot
png("manhattan_meta.png",width=800)
manhattanPlot(p=dmr$z_p,chr=dmr$CHR,signif=9e-08,ylim=c(0,250))
dev.off()
pdf("manhattan_meta.pdf")
manhattanPlot(p=dmr$z_p,chr=dmr$CHR,signif=9e-08,ylim=c(0,250))
dev.off()

