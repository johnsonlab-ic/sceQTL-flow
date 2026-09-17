process finalize_residuals {
    tag "${celltype}"
    label "process_low_memory"
    publishDir "${params.outdir}/final_residuals/", mode: 'copy'

    input:
    tuple val(celltype), path(expression_mat), path(pcs_file)

    output:
    tuple val(celltype), path("*_final_residuals.csv"), path(pcs_file), emit: final_residuals
    path "*_covs_used.rds", emit: covs_used
    path "*_sdY_check.csv", emit: sdy_check

    script:
    """
    #!/usr/bin/env Rscript
    # NOTE: PCs are explicitly regressed out of the expression matrix here
    # (below), AND the PCs file is still forwarded to run_matrixeQTL, which
    # passes it to MatrixEQTL as cvrt. This is not redundant: MatrixEQTL's
    # cvrt orthogonalizes covariates out of BOTH expression and genotype --
    # explicit regression here only ever touched expression. Re-passing PCs
    # to MatrixEQTL is what adjusts genotype for population-structure/PC
    # effects; it's a no-op on the expression side, since residuals of an
    # OLS fit are already exactly orthogonal to the predictors that produced
    # them (regressing an already-orthogonal vector on the same predictors
    # again returns it unchanged). Dropping the genotype-side adjustment
    # (i.e. passing covmat=NULL to MatrixEQTL) was tried and shown to
    # materially change beta/t.stat/p-value/FDR versus the pre-refactor
    # pipeline -- this two-part design is what keeps results identical.
    library(data.table)
    library(dplyr)

    celltype <- "${celltype}"
    exp_mat <- fread("$expression_mat") %>% tibble::column_to_rownames(var="geneid")

    # PCs file has 0 rows when no PCs were selected/used (see the n_pcs=0 fix
    # in optimize_pcs.nf / select_pcs.nf / generate_fixed_pcs.nf).
    pcs_size <- file.info("$pcs_file")\$size
    if (is.na(pcs_size) || pcs_size == 0) {
        n_pcs <- 0
        pcs <- NULL
    } else {
        # No explicit header=TRUE: the PCs file was written by write.table()
        # with row.names=TRUE, so its header line has one fewer field than
        # the data rows (no label for the PC-name column). fread's default
        # "auto" header detection catches this and shifts correctly,
        # naming that first (PC-name) column "V1" -- forcing header=TRUE
        # instead defeats that detection and misreads the whole file.
        pcs_raw <- fread("$pcs_file", data.table = FALSE)
        if (nrow(pcs_raw) == 0) {
            n_pcs <- 0
            pcs <- NULL
        } else {
            pcs <- pcs_raw %>% tibble::column_to_rownames(var = "V1")
            pcs <- as.matrix(pcs[, colnames(exp_mat), drop = FALSE])
            n_pcs <- nrow(pcs)
        }
    }

    saveRDS(pcs, paste0(celltype, "_covs_used.rds"))

    cat(sprintf("[%s] Removing %d PC(s) from residuals prior to final standardization\\n", celltype, n_pcs))

    if (n_pcs > 0) {
        pcs_df <- as.data.frame(t(pcs))
        exp_mat <- t(apply(exp_mat, 1, function(gene_exp) {
            lm_data <- data.frame(gene = gene_exp, pcs_df)
            resid(lm(gene ~ ., data = lm_data))
        }))
    }

    do_standardize <- ${params.standardize_residuals ? 'TRUE' : 'FALSE'}
    if (do_standardize) {
        cat(sprintf("[%s] Centering and scaling residuals to unit variance (SdY=1)\\n", celltype))
        exp_mat <- t(scale(t(exp_mat), center = TRUE, scale = TRUE))
    } else {
        # Always explicitly mean-center (covariate/PC lm() residuals are ~0
        # already, but the no-covariate/no-PC case starts from raw log2CPM+1
        # values, which are not centered at all) so the tested phenotype is
        # on a consistent, well-defined scale either way.
        exp_mat <- t(scale(t(exp_mat), center = TRUE, scale = FALSE))
    }

    # Drop zero-variance genes (defensive; can arise after PC removal)
    gene_sd <- apply(exp_mat, 1, sd, na.rm = TRUE)
    keep <- which(gene_sd > 0 & !is.na(gene_sd))
    n_dropped <- nrow(exp_mat) - length(keep)
    if (n_dropped > 0) {
        cat(sprintf("[%s] Dropping %d zero-variance genes after final residualization\\n", celltype, n_dropped))
        exp_mat <- exp_mat[keep, , drop = FALSE]
    }

    # Record achieved per-gene SD for QC/reporting (should be ~1 everywhere
    # when do_standardize is TRUE; informative either way).
    sd_check <- data.frame(geneid = rownames(exp_mat), sd = apply(exp_mat, 1, sd, na.rm = TRUE))
    fwrite(sd_check, paste0(celltype, "_sdY_check.csv"))

    exp_mat <- as.data.frame(exp_mat) %>% mutate(geneid = row.names(.))
    fwrite(exp_mat, paste0(celltype, "_final_residuals.csv"))
    """
}
