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
GENOMICCONTROL ON     

GENOMICCONTROL 1.13
PROCESS ../results_BET_sex.txt

GENOMICCONTROL 1.37
PROCESS ../results_ITU_sex.txt

GENOMICCONTROL 1.17
PROCESS ../results_PREDO_sex.txt


#Execute random-effects meta-analysis 
OUTFILE meta_BET_PREDO_ITU_random .tbl
ANALYZE RANDOM

QUIT
