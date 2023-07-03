#-------------
#this was run via bash script
#METAL was used in the random effects setting as provided in https://github.com/explodecomputer/random-metal
SCHEME   STDERR
SEPARATOR  TAB

#set-up files
MARKER   CpG
EFFECT   beta_sex
STDERR   se_sex
PVAL     p_sex

PROCESS results_BET_sex.txt
PROCESS results_PREDO_sex.txt
PROCESS results_ITU_sex.txt

#Execute random-effects meta-analysis 
OUTFILE meta_BET_PREDO_ITU_random .tbl
ANALYZE RANDOM

QUIT
