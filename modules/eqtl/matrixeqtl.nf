process run_matrixeQTL {
  
    tag "${expression_mat}"
    label "process_eqtl"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path source_R
    path genotype_mat
    path snp_locations
    path expression_mat
    path gene_locations

    output:
    path "*_cis_eqtl_sig.rds", emit: sig
    path "*_eqtl_summary.rds", emit: summary
    path "*_cis_MatrixEQTLout.rds", emit: full, optional: true

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

    # Update to handle both residuals and normalized file naming patterns
    if (grepl("_final_residuals.csv\$", basename("$expression_mat"))) {
        celltype = gsub("_final_residuals.csv", "", basename("$expression_mat"))
    } else if (grepl("_residuals.csv\$", basename("$expression_mat"))) {
        celltype = gsub("_residuals.csv", "", basename("$expression_mat"))
    } else if (grepl("_pseudobulk_normalised.csv\$", basename("$expression_mat"))) {
        celltype = gsub("_pseudobulk_normalised.csv", "", basename("$expression_mat"))
    } else {
        # Fallback: remove .csv extension
        celltype = gsub(".csv\$", "", basename("$expression_mat"))
    }
    common_samples = intersect(colnames(exp_mat), colnames(geno_mat))
    
    cat(sprintf("[eQTL] %s: Testing %d individuals\\n", celltype, length(common_samples)))
    if (length(common_samples) < 20) {
        cat(sprintf("[eQTL] ⚠️  WARNING: Only %d individuals for %s (eQTL may be underpowered)\\n", length(common_samples), celltype))
    }

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

    # Covariates and PCs have already been regressed out (and, if
    # --standardize_residuals is set, the residuals rescaled to unit
    # variance) upstream in finalize_residuals — nothing left to pass here.
    ##finally, re-order inputs to same column order
    exp_mat = exp_mat[, common_samples]
    geno_mat = geno_mat[, common_samples]

    message("Calculating eQTLs")
    outs=calculate_ciseqtl(
        exp_mat = exp_mat,
        exp_loc = exp_loc,
        geno_mat = geno_mat,
        geno_loc = geno_loc,
        name = celltype,
        covmat = NULL,
        pvOutputThreshold = 0,
        cisDist = as.numeric(${params.cis_distance})
    )

    # Filter + summarise at the source so combine_eqtls never holds the full table.
    outs = if (is.null(outs)) data.frame() else as.data.frame(outs)
    if (nrow(outs) > 0) {
        names(outs)[names(outs) == "statistic"] <- "t.stat"
        names(outs)[names(outs) == "pvalue"]    <- "p.value"
        names(outs)[names(outs) == "snps"]      <- "SNP"
        sig_mask <- outs[["FDR"]] <= as.numeric(${params.fdr_threshold})
    } else {
        sig_mask <- logical(0)
    }

    summ <- data.frame(
        celltype       = celltype,
        n_individuals  = length(common_samples),
        n_genes_tested = if (nrow(outs) > 0) length(unique(outs[["gene"]])) else 0L,
        n_snps_tested  = if (nrow(outs) > 0) length(unique(outs[["SNP"]]))  else 0L,
        n_tests        = nrow(outs),
        n_sig_pairs    = sum(sig_mask),
        n_egenes       = if (nrow(outs) > 0) length(unique(outs[["gene"]][sig_mask])) else 0L,
        stringsAsFactors = FALSE
    )
    saveRDS(summ, paste0(celltype, "_eqtl_summary.rds"))
    saveRDS(outs[sig_mask, , drop = FALSE], paste0(celltype, "_cis_eqtl_sig.rds"))

    # Full per-cell-type cis stats are opt-in (for coloc/mashr); never concatenated.
    if (${params.save_full_eqtl ? 'TRUE' : 'FALSE'}) {
        saveRDS(outs, paste0(celltype, "_cis_MatrixEQTLout.rds"))
    }
    """
}
