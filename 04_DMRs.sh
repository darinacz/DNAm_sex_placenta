#DMRs were calculated using a bash script and comb-p

#bed file needs "chr" and needs to be sorted by chr and pos
#header of bed file: #chrom start end P SNP
#-c 4: p-values are in 4th column
#--seed 1e-08 require a p-value of 1e-8 to start a region
#--dist 200 extend region if find another p-value within this dist
#--p prefix dmr

#CpG has to be at least nominally significant to start a region
comb-p pipeline -c 4 --seed 0.05 --dist 200 -p dmr_seed_nominal meta_compb.bed





