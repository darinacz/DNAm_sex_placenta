library(TwoSampleMR)
library(ieugwasr)
library(RColorBrewer)
library(lessR)

#######extract GWAS results#######
#SNPs to check for associations with GWAS
#per meQTLs, we take the SNP with the lowest p-value
load("consistent_best_SNP_dmps.Rdata")
unique(final_consistent$SNP) -> unique.exposure.snps

#GWAS to check
#this is the list of GWAS we went to check, see methods part of the manuscript
load("IEU_List_Selected_Outcomes.RData")

#extract GWAS results for each SNP
outcome.data.full<-rep(NA,16)
for (i in 1:7503)
{outcome.data <- extract_outcome_data(
  snps = unique.exposure.snps,
  outcomes = ao$id[i],
  proxies = FALSE
)
outcome.data.full<-rbind(outcome.data.full,outcome.data)
print(i)}

#first line is NAs
outcome.data.full<-outcome.data.full[-1,]
save(outcome.data.full, file="outcome.data.full.Rdata")


######perform pheWAS analysis#####
length(unique(outcome.data.full$outcome)) #n=7,415 outcomes
outcome.data.full$FDR<-p.adjust(outcome.data.full$pval.outcome, method="BH")

#save full results
save(outcome.data.full, file="outcome.data.full.Rdata")
length(which(outcome.data.full$FDR<0.05)) #n=14,192

index<-which(outcome.data.full$FDR<0.05)
top<-outcome.data.full[index,]
save(top,file="pheWAS_FDR_0.05.Rdata")
write.table(top,"pheWAS_FDR_0.05.txt",sep="\t",quote=F,row.names=F)
length(unique(top$outcome)) #1,861

#read in list of categories, manually created, see Suppl. Tables
hits<-read.table("categories_1861_traits.txt",sep="\t",header=T)

#sum up number of categories across all PheWAS hits
cat<-unique(hits$Category)
cat<-as.data.frame(cat)
cat$sum<-NA
for (i in 1:17)
{index<-which(hits$Category==cat[i,1])
tmp<-hits[index,]
cat$sum[i]<-sum(tmp$Hits)
}

#recode for piechart
test<-rep(NA,1)
for (i in 1:17)
{test<-c(test,rep(as.character(cat$cat[i]),cat$sum[i]))}
test<-na.omit(test)


#put low categories (below 1%) into extra categories: alcohol, hair, puberty, menopause, skin, sleep,
#smoking
hits1<-test
hits1[hits1=="alcohol"]<-"other"
hits1[hits1=="cognitive_performance"]<-"other"
hits1[hits1=="hair"]<-"other"
hits1[hits1=="menopause"]<-"other"
hits1[hits1=="puberty"]<-"other"
hits1[hits1=="skin"]<-"other"
hits1[hits1=="sleep"]<-"other"
hits1[hits1=="smoking"]<-"other"
length(unique(hits1)) # 9 unique categories

#piechart
color <- brewer.pal(9, "Paired") 
pdf("1861_traits_cat.pdf")
PieChart(hits1, hole = 0, values = "%", 
         col = color, main = "", fill=color, labels_cex = 0.6, main_cex=0.6)
dev.off()
dev.print("1861_traits_cat.pdf",dev=pdf)

######PheWAS based on sex-stratified results#######
##summary stats for the for SNPs to check for associations with GWAS with those traits found to 
#be genetically differentially regualated by sex
#from the Bernabeu et al. paper were downloaded from 
#https://datashare.ed.ac.uk/handle/10283/3917 (non-binary traits) and
#https://datashare.ed.ac.uk/handle/10283/3918 (binary traits)
load("results_nonbinary_traits_meQTLs_consistent.Rdata")
length(unique(results$SNP))
length(unique(results$exp))
results1<-results

load("results_binary_traits_meQTLs_consistent.Rdata")
results2<-results

#combine both
res1<-results1[,c(1,16,17)]
names(res1)<-c("SNP","p","exp")
res2<-results2[,c(1,23,25)]
names(res2)<-c("SNP","p","exp")
res<-rbind(res1,res2)
res<-na.omit(res)

res$FDR<-p.adjust(res$p, method="BH") #FDR-correction
length(which(res$FDR<0.05)) #n=85
index<-which(res$FDR<0.05)
res[index,]
length(unique(res$SNP[index])) #n=41
length(unique(res$exp[index])) #n=39
table(res$SNP[index],res$exp[index])
save(res,file="outcome.data.full.Rdata")

#need to check in trait list what these traits mean
table(res$exp[index])
#non-binary traits
#beta.sexcomparison__48_49_-0.0.tsv.gz: Waist circumference / Hip circumference, 9 hits
#beta.sexcomparison_1687-0.0.tsv.gz: Comparative body size at age 10, 1 hit
#beta.sexcomparison_23098-0.0.tsv.gz: Weight, 1 hit
#beta.sexcomparison_23101-0.0.tsv.gz: Whole body fat-free mass, 3 hits
#beta.sexcomparison_23102-0.0.tsv.gz: Whole body water mass, 3 hits
#beta.sexcomparison_23105-0.0.tsv.gz: Basal metabolic rate, 3 hits
#beta.sexcomparison_23106-0.0.tsv.gz: Impedance of whole body, 2 hits
#beta.sexcomparison_23107-0.0.tsv.gz: Impedance of leg (right), 2 hits 
#beta.sexcomparison_23110-0.0.tsv.gz: Impedance of arm (left), 1 hit
#beta.sexcomparison_23113-0.0.tsv.gz: Leg fat-free mass (right), 1 hit
#beta.sexcomparison_23114-0.0.tsv.gz: Leg predicted mass (right), 1 hit 
#beta.sexcomparison_23115-0.0.tsv.gz: Leg fat percentage (left),1 hit  
#beta.sexcomparison_23117-0.0.tsv.gz: Leg fat-free mass (left), 1 hit 
#beta.sexcomparison_23118-0.0.tsv.gz: Leg predicted mass (left), 1 hit 
#beta.sexcomparison_23119-0.0.tsv.gz: Arm fat percentage (right), 1 hit 
#beta.sexcomparison_23120-0.0.tsv.gz: Arm fat mass (right), 1 hit 
#beta.sexcomparison_23121-0.0.tsv.gz: Arm fat-free mass (right), 3 hits 
#beta.sexcomparison_23122-0.0.tsv.gz: Arm predicted mass (right), 4 hits 
#beta.sexcomparison_23125-0.0.tsv.gz: Arm fat-free mass (left), 3 hits
#beta.sexcomparison_23126-0.0.tsv.gz: Arm predicted mass (left), 3 hits
#beta.sexcomparison_23127-0.0.tsv.gz: Trunk fat percentage, 2 hits 
#beta.sexcomparison_23129-0.0.tsv.gz: Trunk fat-free mass, 3 hits
#beta.sexcomparison_23130-0.0.tsv.gz: Trunk predicted mass, 2 hits
#beta.sexcomparison_30100-0.0.tsv.gz: Mean platelet (thrombocyte) volume, 1 hit 
#beta.sexcomparison_30110-0.0.tsv.gz: Platelet distribution width, 2 hits 
#beta.sexcomparison_30180-0.0.tsv.gz: Lymphocyte percentage, 1 hit
#beta.sexcomparison_30210-0.0.tsv.gz: Eosinophill percentage, 2 hits 
#beta.sexcomparison_50-0.0.tsv.gz: Standing height, 3 hits
#binary traits
#beta.sexcomparison_REGENIE_clinical_c_Block_E00-E07.tsv.gz: disorders of thyroid gland, 2 hits 
#beta.sexcomparison_REGENIE_clinical_c_Block_I20-I25.tsv.gz: Ischaemic heart diseases ,1 hit
#beta.sexcomparison_REGENIE_clinical_c_Block_K40-K46.tsv.gz: hernia, 1 hit  
#beta.sexcomparison_REGENIE_clinical_c_E03.tsv.gz: Other hypothyroidism, 2 hits 
#beta.sexcomparison_REGENIE_clinical_c_K90.tsv.gz: Intestinal malabsorption, 2 hits 
#beta.sexcomparison_REGENIE_selfReported_n_1082.tsv.gz:  heart/cardiac problem, 1 hits
#beta.sexcomparison_REGENIE_selfReported_n_1112.tsv.gz:  deep venous thrombosis, 1 hit 
#beta.sexcomparison_REGENIE_selfReported_n_1249.tsv.gz: thyroid problem (not cancer), 4 hits  
#beta.sexcomparison_REGENIE_selfReported_n_1251.tsv.gz: hypothyroidism/myxoedema, 5 hits 
#beta.sexcomparison_REGENIE_selfReported_n_1516.tsv.gz: malabsorption/coeliac disease, 4 hits 
#beta.sexcomparison_REGENIE_selfReported_n_1526.tsv.gz: gout, 1 hit 
 
#read in manually created trait/domain-list, see Suppl.Tables
traits<-read.table("traits_domains.txt",sep="",header=T)

#link hits to domains
hits<-rbind(cbind("Waist circumference / Hip circumference", "anthropometric", 9),
            cbind("Comparative body size at age 10", "anthropometric", 1 ),
            cbind("Weight", "anthropometric",1),
            cbind("Whole body fat-free mass","anthropometric",3),
            cbind("Whole body water mass","anthropometric",3),
            cbind("Basal metabolic rate","metabolism",3),
            cbind("Impedance of whole body","anthropometric",2),
            cbind("Impedance of leg (right)","anthropometric",2),
            cbind("Impedance of arm (left)","anthropometric",1),
            cbind("Leg fat-free mass (right)","anthropometric",1),
            cbind("Leg predicted mass (right)","anthropometric",1),
            cbind("Leg fat percentage (left)","anthropometric", 1),
            cbind("Leg fat-free mass (left)","anthropometric",1),
            cbind("Leg predicted mass (left)","anthropometric", 1),
            cbind("Arm fat percentage (right)","anthropometric", 1),
            cbind("Arm fat mass (right)","anthropometric", 1),
            cbind("Arm fat-free mass (right)","anthropometric", 3),
            cbind("Arm predicted mass (right)","anthropometric",4),
            cbind("Arm fat-free mass (left)","anthropometric", 3),
            cbind("Arm predicted mass (left)","anthropometric",3),
            cbind("Trunk fat percentage","anthropometric",3),
            cbind("Trunk fat-free mass","anthropometric", 3),
            cbind("Trunk predicted mass","anthropometric", 2),
            cbind("Mean platelet (thrombocyte) volume","blood_composition",1),
            cbind("Platelet distribution width","blood_composition",2),
            cbind("Lymphocyte percentage","blood_composition",1),
            cbind("Eosinophill percentage","blood_composition",2),
            cbind("Standing height","anthropometric",3),
            cbind("Disorders of thyroid gland","disease",2),
            cbind("Ischaemic heart diseases","disease",1),
            cbind("Hernia","disease",1 ),
            cbind("Other hypothyroidism","disease",2 ),
            cbind("Intestinal malabsorption","disease",2),
            cbind("Heart/cardiac problem","disease", 1),
            cbind("Deep venous thrombosis","disease", 1),
            cbind("Thyroid problem (not cancer)","disease", 4),
            cbind("Hypothyroidism/myxoedema","disease", 5),
            cbind("Malabsorption/coeliac disease","disease", 4),
            cbind("Gout","disease",1))

hits<-as.data.frame(hits)
names(hits)<-c("trait","Category","hits")
hits$hits<-as.numeric(hits$hits)

cat<-unique(hits$Category)
cat<-as.data.frame(cat)
cat$sum<-NA
for (i in 1:4)
{index<-which(hits$Category==cat[i,1])
tmp<-hits[index,]
cat$sum[i]<-sum(tmp$hits)
}

test<-rep(NA,1)
for (i in 1:4)
{test<-c(test,rep(as.character(cat$cat[i]),cat$sum[i]))}
test<-na.omit(test)

pdf("39_traits_cat.pdf")
PieChart(test, hole = 0, values = "%",
         col=c("#A6CEE3","#1F78B4","#B2DF8A","#FDBF6F"), fill = c("#A6CEE3","#1F78B4","#33A02C","#FDBF6F"), labels_cex = 0.6, main="", main_cex=0.6)
dev.off()











