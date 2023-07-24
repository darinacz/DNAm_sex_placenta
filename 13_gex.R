library(DESeq2)
library(ggplot2)
library(corrplot)

#check differential gene expression with sex in ITU placenta 
#load gex data
load("dds_sex.Rdata")
pt<-dds_sex@colData

counts<-counts(dds,norm=F)

dds_sex <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = pt,
                                   design = ~ SV1+Gestational_Age_Weeks+smoking+labor+Caesarian_Section+sex)
dds_sex<-estimateSizeFactors(dds_sex)
dds_sex<-DESeq(dds_sex)

save(dds_sex, file="dds_sex.Rdata")
results_sex<-results(dds_sex) 
index<-which(results_sex$padj<0.05) #n=50
results_sex[index,]
write.table(results_sex, "results_deseq2_sex_autosomes.txt",sep="\t",quote=F,row.names=F)

###volcano plot
res<-results_sex

res$expression = "stable"
res$expression[res$padj < 0.05 & res$log2FoldChange >= 0]<-"up"
res$expression[res$padj < 0.05 & res$log2FoldChange<0]<-"down"
p <- ggplot(data = as.data.frame(res), 
            aes(x = log2FoldChange, 
                y = -log10(res$pvalue), 
                colour=expression)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-0.4, 0.4)) +
  labs(x="effect size",
       y="-log10 (p-value)",
       title="Differential gene expression")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())
dev.print("volcano_plac.png",dev=png,width=800)
dev.print("volcano_plac.pdf",dev=pdf)

table(sign(results_sex$log2FoldChange))
#-1    1 
#4065 3890 

write.table(results_sex$gene,"7955_autosomal_genes.txt",sep="\t",quote=F,row.names=F)
index<-which(results_sex$padj<0.05) #n=50
write.table(results_sex$gene[index],"50_deg_sex.txt",sep="\t",quote=F,row.names=F)

