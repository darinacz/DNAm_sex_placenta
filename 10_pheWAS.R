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
length(which(outcome.data.full$FDR<0.05)) #n=12,035

index<-which(outcome.data.full$FDR<0.05)
top<-outcome.data.full[index,]
save(top,file="pheWAS_FDR_0.05.Rdata")
write.table(top,"pheWAS_FDR_0.05.txt",sep="\t",quote=F,row.names=F)
length(unique(top$outcome)) #1,1699

#read in list of categories, manually created, see Suppl. Table S11
hits<-read.table("categories_1698_traits.txt",sep="\t",header=T)

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
hits1[hits1=="skin"]<-"other"
hits1[hits1=="hair"]<-"other"
hits1[hits1=="puberty"]<-"other"
hits1[hits1=="menopause"]<-"other"
hits1[hits1=="cognition"]<-"other"
hits1[hits1=="exercise"]<-"other"
hits1[hits1=="urine_composition"]<-"other"
hits1[hits1=="sex"]<-"other"
length(unique(hits1)) # 10 unique categories

####Figure 2B#####
#piechart
color <- brewer.pal(10, "Paired") 
pdf("1698_traits_cat.pdf")
PieChart(hits1, hole = 0, values = "%", 
         col = color, main = "", fill=color, labels_cex = 1.2, main_cex=1.2, values_position="label")
dev.off()


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
length(which(res$FDR<0.05)) #n=87
index<-which(res$FDR<0.05)
res[index,]
length(unique(res$SNP[index])) #n=41
length(unique(res$exp[index])) #n=41
table(res$SNP[index],res$exp[index])
save(res,file="outcome.data.full.Rdata")

#need to check in trait list what these traits mean, Table S12
table(res$exp[index])
#non-binary traits
#beta.sexcomparison__48_49_-0.0.tsv.gz: Waist circumference / Hip circumference, 8 hits
#beta.sexcomparison_1687-0.0.tsv.gz: Comparative body size at age 10	Non-Binary, 1 hit
#beta.sexcomparison_23098-0.0.tsv.gz: Weight, 3 hits 
#beta.sexcomparison_23101-0.0.tsv.gz: Whole_body_fat-free_mass, 3 hits
#beta.sexcomparison_23102-0.0.tsv.gz: Whole body water mass, 2 hits
#beta.sexcomparison_23105-0.0.tsv.gz: Basal metabolic rate, 3 hits
#beta.sexcomparison_23106-0.0.tsv.gz: Impedance of whole body, 3 hits
#beta.sexcomparison_23107-0.0.tsv.gz: Impedance of leg (right), 1 hits 
#beta.sexcomparison_23108-0.0.tsv.gz: Impedance of leg (left), 1 hits 
#beta.sexcomparison_23110-0.0.tsv.gz: Impedance of arm (left), 1 hits 
#beta.sexcomparison_23113-0.0.tsv.gz: Leg fat-free mass (right), 1 hit
#beta.sexcomparison_23114-0.0.tsv.gz: Leg predicted mass (right), 1 hit 
#beta.sexcomparison_23115-0.0.tsv.gz: Leg fat percentage (left),1 hit  
#beta.sexcomparison_23117-0.0.tsv.gz: Leg_fat-free_mass_left,1 hit 
#beta.sexcomparison_23118-0.0.tsv.gz: Leg_predicted_mass_left,1 hit 
#beta.sexcomparison_23119-0.0.tsv.gz: Arm_fat_percentage_right,1 hit 
#beta.sexcomparison_23121-0.0.tsv.gz: Arm fat-free mass (right), 3 hits 
#beta.sexcomparison_23122-0.0.tsv.gz: Arm predicted mass (right), 4 hits 
#beta.sexcomparison_23125-0.0.tsv.gz: Arm fat-free mass (left), 3 hits
#beta.sexcomparison_23126-0.0.tsv.gz: Arm predicted mass (left), 3 hits
#beta.sexcomparison_23127-0.0.tsv.gz: Trunk_fat_percentage anthropometric, 2 hits
#beta.sexcomparison_23129-0.0.tsv.gz: Trunk fat-free mass, 4 hits
#beta.sexcomparison_23130-0.0.tsv.gz: Trunk predicted mass, 2 hits
#beta.sexcomparison_30100-0.0.tsv.gz: Mean_platelet_thrombocyte_volume, 1 hit
#beta.sexcomparison_30110-0.0.tsv.gz: Platelet_distribution_width blood_composition, 4 hits
#beta.sexcomparison_30180-0.0.tsv.gz: Lymphocyte_percentage blood_composition, 1 hit
#beta.sexcomparison_30210-0.0.tsv.gz: Eosinophill_percentage blood_composition, 2 hits
#beta.sexcomparison_50-0.0.tsv.gz: Standing height, 3 hits
#binary traits
#beta.sexcomparison_REGENIE_clinical_c_Block_E00-E07.tsv.gz: disorders of thyroid gland, 1 hits 
#beta.sexcomparison_REGENIE_clinical_c_Block_I20-I25.tsv.gz: Ischaemic_heart_diseases, 1 hit
#beta.sexcomparison_REGENIE_clinical_c_Block_K40-K46.tsv.gz: Hernia, 1 hit
#beta.sexcomparison_REGENIE_clinical_c_Block_K90-K93.tsv.gz: Other_diseases_of_the_digestive_system, 1 hit
#beta.sexcomparison_REGENIE_clinical_c_E03.tsv.gz: Other hypothyroidism, 1 hit
#beta.sexcomparison_REGENIE_clinical_c_K90.tsv.gz: Intestinal malabsorption, 1 hit
#beta.sexcomparison_REGENIE_selfReported_n_1082.tsv.gz: Heart_cardiac_problem, 1 hit
#beta.sexcomparison_REGENIE_selfReported_n_1112.tsv.gz: Deep_venous_thrombosis, 1 hit
#beta.sexcomparison_REGENIE_selfReported_n_1249.tsv.gz: thyroid problem (not cancer), 1 hit  
#beta.sexcomparison_REGENIE_selfReported_n_1251.tsv.gz: hypothyroidism/myxoedema, 4 hits
#beta.sexcomparison_REGENIE_selfReported_n_1284.tsv.gz: Chronic_degenerative_neurological_problem, 1 hit
#beta.sexcomparison_REGENIE_selfReported_n_1516.tsv.gz: Malabsorption_coeliac_disease, 4 hits
#beta.sexcomparison_REGENIE_selfReported_n_1526.tsv.gz: Gout, 1 hit

 










