library(Biobase)
library(foreign)
library(MASS)
library(lmtest)
library(sandwich)

######Is DNAm significantly different between males and females in the BET cohort?#####

#load final, normalized and batch-corrected methylation data
load("BetaSet_BET_137.rda")
beta<-exprs(Beta_BET_137) 
pt<-pData(Beta_BET_137)

#run rlm using White's estimator and store results
results<-matrix(ncol=4,nrow=dim(beta)[1])
for (i in 1: dim(beta)[1])
  {
  results[i,1]<-rownames(beta)[i]
  #rlm per CpG methylation on sex using gestational age, maternal smoking and cell types as covariates
  #cell type proportions add up to 1, so one celltype is omitted
  mod<-rlm(beta[i,]~pt$gestage_days+pt$smoking+pt$Trophoblasts+pt$Stromal+pt$Endothelial+
             pt$nRBC+pt$Syncytiotrophoblast+pt$sex)
cf = try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))
if (class(cf)=="try-error") {
  results[i,2]<-"NA"
  results[i,3]<-"NA"
  results[i,4]<-"NA"
}
else{
  results[i,2]<-cf[9,1]
  results[i,3]<-cf[9,2]
  results[i,4]<-cf[9,4]
}
print(i)
}

results<-as.data.frame(results)
names(results)<-c("CpG","beta_sex","se_sex","p_sex")
for (i in 2:4)
{results[,i]<-as.numeric(as.character(results[,i]))}

#save results
results->results_BET_sex
save(results_BET_sex,file="results_BET_sex.Rdata")
write.table(results_BET_sex,"results_BET_sex.txt",sep="\t",quote=F,row.names=F)

#####Is DNAm significantly different between males and females in the ITU cohort?#####

#load final, normalized and batch-corrected methylation data
load("BetaSet_ITU_470.rda")
beta<-exprs(Beta_ITU_470) 
pt<-pData(Beta_ITU_470)

#run rlm using White's estimator
results<-matrix(ncol=4,nrow=dim(beta)[1])
for (i in 1: dim(beta)[1])
  {
  results[i,1]<-rownames(beta)[i]
  #rlm per CpG methylation on sex using gestational age, maternal smoking and cell types as covariates
  #cell type proportions add up to 1, so one celltype is omitted
  mod<-rlm(beta[i,]~pt$gestage_days+pt$smoking+pt$Trophoblasts+pt$Stromal+pt$Endothelial+
             pt$nRBC+pt$Syncytiotrophoblast+pt$sex)
cf = try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))
if (class(cf)=="try-error") {
  results[i,2]<-"NA"
  results[i,3]<-"NA"
  results[i,4]<-"NA"
}
else{
  results[i,2]<-cf[9,1]
  results[i,3]<-cf[9,2]
  results[i,4]<-cf[9,4]
}
print(i)
}

results<-as.data.frame(results)
names(results)<-c("CpG","beta_sex","se_sex","p_sex")
for (i in 2:4)
{results[,i]<-as.numeric(as.character(results[,i]))}

#save results
results->results_ITU_sex
save(results_ITU_sex,file="results_ITU_sex.Rdata")
write.table(results_ITU_sex,"results_ITU_sex.txt",sep="\t",quote=F,row.names=F)

#####Is DNAm significantly different between males and females in the PREDO cohort?#####

#load final, normalized and batch-corrected methylation data
load("BetaSet_PREDO_139.rda")
beta<-exprs(Beta_PREDO_139)
pt<-pData(Beta_PREDO_139)

#run rlm using White's estimator

results<-matrix(ncol=4,nrow=dim(beta)[1])
for (i in 1: dim(beta)[1])
  {
  results[i,1]<-rownames(beta)[i]
  #rlm per CpG methylation on sex using gestational age, maternal smoking and cell types as covariates
  #cell type proportions add up to 1, so one celltype is omitted
  #etimated Hofbauer cell proportions were all at 0 and hence only four cell types are used here
  mod<-rlm(beta[i,]~pt$gestage_days+pt$smoking+pt$Trophoblasts+pt$Stromal+pt$Endothelial+
             pt$nRBC+pt$sex)
cf = try(coeftest(mod, vcov=vcovHC(mod, type="HC0")))
if (class(cf)=="try-error") {
  results[i,2]<-"NA"
  results[i,3]<-"NA"
  results[i,4]<-"NA"
}
else{
  results[i,2]<-cf[8,1]
  results[i,3]<-cf[8,2]
  results[i,4]<-cf[8,4]
}
print(i)
}

results<-as.data.frame(results)
names(results)<-c("CpG","beta_sex","se_sex","p_sex")
for (i in 2:4)
{results[,i]<-as.numeric(as.character(results[,i]))}

#save results
results->results_PREDO_sex
save(results_PREDO_sex,file="results_PREDO_sex.Rdata")
write.table(results_PREDO_sex,"results_PREDO_sex.txt",sep="\t",quote=F,row.names=F)










