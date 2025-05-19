process run_matrixeQTL {

  
    tag "${expression_mat}"
    label "process_high_memory"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path source_R
    path genotype_mat
    path snp_locations
    path expression_mat
    path gene_locations
    path cov_file
    tuple val(celltype_optimal), path(optimal_pc_file), optional: true

    output:
    path "*_cis_MatrixEQTLout.rds", emit: eqtl_results
    path "*"

    script:
    """
    #!/usr/bin/env Rscript
    source("$source_R")
    library(data.table)
    library(dplyr)

    exp_mat = fread("$expression_mat") %>% tibble::column_to_rownames(var="geneid")
    geno_mat = fread("$genotype_mat") %>% tibble::column_to_rownames(var="snp")
    geno_loc = fread("$snp_locations")
    exp_loc = fread("$gene_locations")

    celltype = gsub("_pseudobulk_normalised.csv", "", "$expression_mat")
    common_samples = intersect(colnames(exp_mat), colnames(geno_mat))

    exp_mat = exp_mat %>% select(all_of(common_samples))
    geno_mat = geno_mat %>% select(all_of(common_samples))

    common_genes = intersect(exp_loc %>% pull(geneid), rownames(exp_mat))
    exp_mat = exp_mat %>% filter(rownames(exp_mat) %in% common_genes)
    exp_loc = exp_loc %>% filter(geneid %in% common_genes)

    geno_loc = geno_loc[, c("annot", "chrom", "position")] %>% tibble::column_to_rownames(var="annot")
    geno_mat = geno_mat[rownames(geno_loc), ]
    geno_mat = geno_mat[complete.cases(geno_mat), ]
    geno_loc = geno_loc[rownames(geno_mat), ]
    geno_loc = geno_loc %>% mutate(annot = rownames(geno_loc)) %>% select(annot, chrom, position)
    
    # Check if we have optimal PCs for this cell type
    celltype = gsub("_pseudobulk_normalised.csv", "", "$expression_mat")
    optimal_pc_file = "$optimal_pc_file"
    
    if(file.exists(optimal_pc_file) && nchar(optimal_pc_file) > 0) {
        message("Using optimal PCs file: ", optimal_pc_file)
        # Read the optimal PCs file and use as covariates
        optimal_pcs = read.table(optimal_pc_file, header=TRUE, row.names=1)
        
        # Combine with any existing covariates
        if(file.exists("$cov_file") && file.size("$cov_file") > 0) {
            existing_covs = read.table("$cov_file", header=TRUE, row.names=1)
            covmat = rbind(existing_covs, optimal_pcs)
        } else {
            covmat = optimal_pcs
        }
    } else {
        message("No optimal PCs file provided for this cell type")
        # Use original covariates if they exist
        if(file.exists("$cov_file") && file.size("$cov_file") > 0) {
            covmat = read.table("$cov_file", header=TRUE, row.names=1)
        } else {
            covmat = NULL
        }
    }

    message("calculating eQTLs")
  

    outs=calculate_ciseqtl(
        exp_mat = exp_mat,
        exp_loc = exp_loc,
        geno_mat = geno_mat,
        geno_loc = geno_loc,
        name = celltype,
        pvOutputThreshold =0,
        cisDist = as.numeric(${params.cis_distance}),
        optimize_pcs = as.logical("${params.optimize_pcs}")
    )

    save_eqtls <- function(eqtls, prefix) {
        if (!is.null(eqtls) && nrow(eqtls) > 0) {
          names(eqtls)[names(eqtls) == "statistic"] <- "t.stat"
          names(eqtls)[names(eqtls) == "pvalue"] <- "p.value"
          names(eqtls)[names(eqtls) == "snps"] <- "SNP"
          saveRDS(eqtls, paste0(celltype, "_cis_MatrixEQTLout.rds"))
        }
      }

    save_eqtls(outs)
    
    """
}

// docker run -it --rm -v /var/lib/docker/alex_tmp/data/:/mnt/data ah3918/expression_image:latest

// docker run -it --rm -v ./:/ ah3918/expression_image:latest