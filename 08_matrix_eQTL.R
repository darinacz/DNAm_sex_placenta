#matrix eQTL was run for each cohort separately (BET, ITU, PREDO) in females and males as well as in males and females only. 
#Afterwards, we meta-analysed the results (all, females, males) using METAL in the same set-up like is given in 02_meta-analysis.sh
#as example, we provide here the code for matrixeQTL in females in ITU, all other matrixEQTL analyses were set up accordingly

library(MatrixEQTL)

base_dir <- getwd()

#set thresholds and parameters
cis_threshold <- 5e-2 #cis-meqtlS cutoff=0.05
cis_dist <- 1.5e5 #Distance for local (cis) gene-SNP pairs: cis window of 100kb
useModel = modelLINEAR; #linear model to model the effect of the genotype as additive linear and test for its significance using t-statistic.
errorCovariance = numeric() # Error covariance matrix

#set file names and paths
SNP_file_name <-  "../genetic_data/geno_itu_females.txt";
methylation_file_name <- "../methylation_data/methylation_itu_females.txt";
covariates_file_name <- "../cov_data/cov_itu_females.txt";
snps_location_file_name <- "../../ITU_all/snplocation.txt";
cpg_location_file_name <- "../../ITU_all/cpglocation.txt";

SNP_file_path <- file.path(base_dir, SNP_file_name)
methylation_file_path <- file.path(base_dir, methylation_file_name)
covariates_file_path <- file.path(base_dir, covariates_file_name)
snps_location_file_path <- file.path(base_dir, snps_location_file_name)
cpg_location_file_path <- file.path(base_dir, cpg_location_file_name)

output_file_name_cis = "meqtl_placenta_cis_itu_females.txt"

#Load genotype data
#one column per ID, one row per SNP (coded as 0,1,2)
snps <- SlicedData$new()
snps$fileDelimiter = "\t" # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values
snps$fileSkipRows = 1 # one row of column labels
snps$fileSkipColumns = 1 # one column of row labels
snps$fileSliceSize = 2000 # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)


#Load gene methylation data
#one column per ID, one row per CpG (M-values)
#same ID order as in genotype data file
gene = SlicedData$new()
gene$fileDelimiter = "\t" # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1 # one row of column labels
gene$fileSkipColumns = 1 # one column of row labels
gene$fileSliceSize = 2000 # read file in slices of 2,000 rows
gene$LoadFile(methylation_file_name)

#Load covariates
#one column per ID, one row per covariate
#same ID order as in genotype data file and methylation data file
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values;
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1     # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name)
}

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
cpgpos = read.table(cpg_location_file_name, header = TRUE, stringsAsFactors = FALSE)

#run analysis
me_meqtl_placenta = Matrix_eQTL_main(
  snps = snps, # SlicedData object with genotype information. 
  gene = gene, # SlicedData object with gene methylation information. 
  cvrt = cvrt, # SlicedData object with additional covariates.
  snpspos = snpspos, #3 columns - SNP name, chromosome, and position
  genepos = cpgpos, #4 columns - CpG name, chromosome, and positions of the left and right ends
  output_file_name.cis = output_file_name_cis, #local associations are saved to this file
  pvOutputThreshold.cis = cis_threshold, #cis threshold for reporting
  pvOutputThreshold = 0, #trans threshold or 0 for only cis
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = FALSE,
  cisDist = cis_dist, 
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE, 
  noFDRsaveMemory = FALSE) 

head(me_meqtl_placenta$cis$eqtls)
save(me_meqtl_placenta, file = "meqtl_itu_females.Rdata")
