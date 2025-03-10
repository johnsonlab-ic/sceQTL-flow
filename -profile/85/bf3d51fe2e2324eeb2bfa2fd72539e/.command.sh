#!/usr/bin/env Rscript
library(dplyr)
source("genotype_functions.r")
generate_genotype_matrix(gds_file="final_geno_440samples.gds", 
chain_fpath="/usr/local/src/hg19ToHg38.over.chain")
