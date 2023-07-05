library(data.table)
library(stringr)
library(VennDiagram)
library(ieugwasr)
library(TwoSampleMR)


########read in meta-analysis results from sex-stratifed meQTLs analysis########
#males
males<-fread("meta_BET_ITU_PREDO_males_random1.tbl", header=T)
males$P<-as.numeric(males$PvalueARE)
length(which(males$P<1e-08)) #n=2,150,477 combinations, this is based on the threshold from Min et al.
index<-which(males$P<1e-08)
male_sub<-males[index,]
table(male_sub$Direction)
#we only take these meQTLs were at least information from two studies is available
index<-which(male_sub$Direction=="-??" | male_sub$Direction=="?-?" |male_sub$Direction=="??-" |
               male_sub$Direction=="??+" | male_sub$Direction=="?+?" | male_sub$Direction=="+??")  
male_sub<-male_sub[-index,] #n=1,634,651 combinations
male_sub$CpG<-str_split_fixed(male_sub$MarkerName, "_", 2)[,1]
male_sub$SNP<-str_split_fixed(male_sub$MarkerName, "_", 2)[,2]
length(unique(male_sub$SNP))#n=721,621 SNPs
length(unique(male_sub$CpG))#n=42,751 CpGs
names(male_sub)[18]<-"P_random_effects"
male_sub<-male_sub[,-19]
save(male_sub,file="meta_males_below_e-08.Rdata")

#females
females<-fread("meta_BET_ITU_PREDO_females_random1.tbl", header=T)
females$P<-as.numeric(females$PvalueARE)
length(which(females$P<1e-08)) #n=1,945,122 combinations, this is based on the threshold from Min et al.
index<-which(females$P<1e-08)
female_sub<-females[index,]
table(female_sub$Direction)
#we only take these were at least information from two studies is available
index<-which(female_sub$Direction=="-??" | female_sub$Direction=="?-?" |female_sub$Direction=="??-" |
               female_sub$Direction=="??+" | female_sub$Direction=="?+?" | female_sub$Direction=="+??")  
female_sub<-female_sub[-index,] #n=1,461,207 combinations
female_sub$CpG<-str_split_fixed(female_sub$MarkerName, "_", 2)[,1]
female_sub$SNP<-str_split_fixed(female_sub$MarkerName, "_", 2)[,2]
length(unique(female_sub$SNP))#n=664,352 SNPs
#how many overlap with males?
length(unique(female_sub$CpG))#n=39,918 CpGs
#how many overlap with males?
names(female_sub)[18]<-"P_random_effects"
female_sub<-female_sub[,-19]
save(female_sub,file="meta_females_below_e-08.Rdata")

########read in meta-analysis results from combined meQTLs analysis########
males<-fread("meta_BET_ITU_PREDO_random1.tbl", header=T)
all$P<-as.numeric(all$PvalueARE)
length(which(all$P<1e-08)) #n=4,154,284 combinations, this is based on the threshold from Min et al.
index<-which(all$P<1e-08)
all_sub<-all[index,]
table(all_sub$Direction)
#we only take these were at least information from two studies is available
index<-which(all_sub$Direction=="-??" | all_sub$Direction=="?-?" |all_sub$Direction=="??-" |
               all_sub$Direction=="??+" | all_sub$Direction=="?+?" | all_sub$Direction=="+??")  
all_sub<-all_sub[-index,] #n=3,040,103 combinations
all_sub$CpG<-str_split_fixed(all_sub$MarkerName, "_", 2)[,1]
all_sub$SNP<-str_split_fixed(all_sub$MarkerName, "_", 2)[,2]
length(unique(all_sub$SNP))#n=1,087,708 SNPs
length(unique(all_sub$CpG))#n=70,855 CpGs
names(all_sub)[18]<-"P_random_effects"
all_sub<-all_sub[,-19]
save(all_sub,file="meta_below_e-08.Rdata")

########merge with results from males and females to check for consistency of effect directions across sexes######
all<-merge(all_sub,males, by="MarkerName", all.x=T)
all1<-all[,c(1,7,12,13,14,18,19,20,26,31,32,33,37)]
names(all1)<-c("MarkerName","Direction","EffectARE","StdErrARE","PvalueARE", "P_random_effects",
               "CpG","SNP","Direction_males","EffectARE_males","StdErrARE_males","PvalueARE_males","P_random_effects_males")

all1<-merge(all1,females, by="MarkerName", all.x=T)
all1<-all1[,c(1:13,19,24,25,26,30)]
names(all1)<-c("MarkerName","Direction","EffectARE","StdErrARE","PvalueARE", "P_random_effects",
               "CpG","SNP","Direction_males","EffectARE_males","StdErrARE_males","PvalueARE_males","P_random_effects_males",
               "Direction_females","EffectARE_females","StdErrARE_females","PvalueARE_females","P_random_effects_females")
all<-all1
save(all,file="meta_below_e-08_with_males_females.Rdata")


#####define sex-consistent and sex-inconsistent CpGs####
table(sign(all$EffectARE_males),sign(all$EffectARE_females))
#         -1       1
#-1 1551035      91
#1      135 1473604
#most are consistent across males and females, remove inconsistent ones
index<-which(sign(all$EffectARE_males)== 1 & sign(all$EffectARE_females)== -1 |
               sign(all$EffectARE_males)== -1 & sign(all$EffectARE_females)== 1) #n=226
consistent<-all[-index,]
               
#identify those which are different between the sexes
males_females<-merge(male_sub,females,by="MarkerName",all.x=T) 
males_females<-males_females[,c(1,7,12,13,14,18,19,20,26,31,32,33,37)]
names(males_females)<-c("MarkerName","Direction_males","EffectARE_males","StdErrARE_males","PvalueARE_males", "P_random_effects_males",
               "CpG","SNP","Direction_females","EffectARE_females","StdErrARE_females","PvalueARE_females","P_random_effects_females")

table(sign(males_females$EffectARE_males),sign(males_females$EffectARE_females))
#       -1      1
#-1 839904    263
#1     387 783071
#mostly consistent, get inconsistent ones
index<-which(sign(males_females$EffectARE_males)==1 & sign(males_females$EffectARE_females) == -1 |
               sign(males_females$EffectARE_males)== -1 & sign(males_females$EffectARE_females) == 1) #n=650
inconsistent_male_female<-males_females[index,]

females_males<-merge(female_sub,males,by="MarkerName",all.x=T) 
females_males<-females_males[,c(1,7,12,13,14,18,19,20,26,31,32,33,37)]
names(females_males)<-c("MarkerName","Direction_females","EffectARE_females","StdErrARE_females","PvalueARE_females", "P_random_effects_females",
                        "CpG","SNP","Direction_males","EffectARE_males","StdErrARE_males","PvalueARE_males","P_random_effects_males")

table(sign(females_males$EffectARE_females),sign(females_males$EffectARE_males))
#       -1      1
#-1 740864    234
#1     261 712302
#mostly consistent, get inconsistent ones
index<-which(sign(females_males$EffectARE_females)==1 & sign(females_males$EffectARE_males) == -1 |
               sign(females_males$EffectARE_females)== -1 & sign(females_males$EffectARE_males) == 1) #n=650
inconsistent_female_male<-females_males[index,]

#For those CpGs for which meQTLs in the sex-consisten and sex-inconsistent groups are present, 
#we use the meQTL with the lowest p-value to decide for one group 
inconsistent<-rbind(inconsistent_female_male[,c("CpG","SNP","P_random_effects_females")],inconsistent_male_female[,c("CpG","SNP","P_random_effects_males")],use.names=F)
names(inconsistent)[3]<-"P_random_effects"

index<-which(consistent$CpG %in% inconsistent$CpG)
both<-consistent[index,]
CpGs<-unique(both$CpG) #n=167 CpGs in both
groups<-cbind(NA,NA)
for (i in 1: length(CpGs))
  {index<-which(both$CpG %in% CpGs[i])
  p1<-min(both$P_random_effects[index])
  index1<-which(inconsistent$CpG %in% CpGs[i])
  p2<-min(inconsistent$P_random_effects[index1])
  if (p1<p2) {groups<-rbind(groups,cbind(CpGs[i],"cons"))}
  else (groups<-rbind(groups,cbind(CpGs[i],"incons")))
print(i)}

table(groups[,2])
#cons incons 
#160      7 

index<-which(groups[,2]=="incons") #n=7
cpgs_rem_consistent<-groups[index,1]
index<-which(consistent$CpG %in% cpgs_rem_consistent)
consistent<-consistent[-index,] #n=3,039,847 combinations
length(unique(consistent$CpG)) #n=70,848 CpGs
length(unique(consistent$SNP)) #n=1,087,627 SNPs

index<-which(groups[,2]=="cons") #n=160
cpgs_rem_inconsistent<-groups[index,1]
index<-which(inconsistent$CpG %in% cpgs_rem_inconsistent)
inconsistent<-inconsistent[-index,] #n=173 combinations
length(unique(inconsistent$CpG)) #n=22 CpGs
length(unique(inconsistent$SNP)) #n=169 SNPs


######check overlap with sex-differential DMPs####
dmps<-read.table("meta_BET_PREDO_ITU_random1.tbl",sep="\t",header=T)
index<-which(dmps$PvalueARE<9e-08) #n=10,737 DMPs
top<-dmps[index,] #epigenome-wide significant DMPs
other<-dmps[-index,] #non epigenome-wide significant CpGs
length(which(top$MarkerName %in% unique(consistent$CpG))) #n=2,262

#test for enrichment for CpGs involved in sex-consistent meQTLs: significant enrichment
a=length(which(top$MarkerName %in% unique(consistent$CpG))) #2,262
b=length(which(!(top$MarkerName %in% unique(consistent$CpG)))) #8,475
c=length(which(other$MarkerName %in% unique(consistent$CpG)))#68,586
d=length(which(!(other$MarkerName %in% unique(consistent$CpG))))#678,778
fisher.test(matrix(c(a,b,c,d),byrow=T,ncol=2)) #p=8.384236e-297

#test for enrichment for CpGs involved in sex-inconsistent meQTLs: no significant enrichment
a=length(which(top$MarkerName %in% unique(inconsistent$CpG))) #0
b=length(which(!(top$MarkerName %in% unique(inconsistent$CpG)))) #10,737
c=length(which(other$MarkerName %in% unique(inconsistent$CpG)))#22
d=length(which(!(other$MarkerName %in% unique(inconsistent$CpG))))#747,342
fisher.test(matrix(c(a,b,c,d),byrow=T,ncol=2)) #p=1

save(consistent,file="sex_consistent_meQTLs.Rdata")
save(inconsistent,file="sex_specific_meQTLs.Rdata")




 












