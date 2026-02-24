process get_residuals {
    tag "${expression_mat}"
    label "process_high"
    publishDir "${params.outdir}/residuals/", mode: 'copy'

    input:

    path expression_mat
    path cov_file
    path source_R
    val covariates_to_include

    output:
    path "*_residuals.csv", emit: residuals_results

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(dplyr)

    exp_mat = fread("$expression_mat") %>% tibble::column_to_rownames(var="geneid")
    celltype = gsub("_pseudobulk_normalised.csv", "", "$expression_mat")
    cat(sprintf("[%s] Loaded expression matrix: %d genes × %d individuals\\n", celltype, nrow(exp_mat), ncol(exp_mat)))
    
    cov_file="$cov_file"
    if(file.size(cov_file) > 0){
        cov_mat = fread(cov_file, head=T) %>% as.data.frame()
        cat(sprintf("[%s] Loaded covariate matrix: %d individuals\\n", celltype, nrow(cov_mat)))
        
        # Normalize sample IDs to handle different separators (dash, slash, period)
        normalize_ids <- function(ids) {
            gsub("[-/.]", "_", as.character(ids))
        }
        
        exp_samples_norm <- normalize_ids(colnames(exp_mat))
        cov_samples_norm <- normalize_ids(cov_mat[["Individual_ID"]])
        
        # Find common samples
        common_norm <- intersect(exp_samples_norm, cov_samples_norm)
        exp_idx <- which(exp_samples_norm %in% common_norm)
        cov_idx <- match(exp_samples_norm[exp_idx], cov_samples_norm)
        
        cat(sprintf("[%s] Common samples: %d / %d\\n", celltype, length(common_norm), ncol(exp_mat)))
        
        # Subset to common samples in matching order
        exp_mat <- exp_mat[, exp_idx]
        cov_mat <- cov_mat[cov_idx, ]
        
        # Parse covariates to include
        covs_to_include = unlist(strsplit("${covariates_to_include}", ","))
        if(length(covs_to_include) == 0 || (length(covs_to_include) == 1 && covs_to_include[1] == "all")) {
            covs_to_include = colnames(cov_mat)[!colnames(cov_mat) %in% c("Individual_ID")]
        }
        cat(sprintf("[%s] Regressing out: %s\\n", celltype, paste(covs_to_include, collapse=", ")))
        
        # Relevel factors to use most frequent level as reference
        for (col_name in covs_to_include) {
            if (is.character(cov_mat[[col_name]])) {
                cov_mat[[col_name]] <- as.factor(cov_mat[[col_name]])
            }
            if (is.factor(cov_mat[[col_name]])) {
                ref_level <- names(sort(table(cov_mat[[col_name]]), decreasing = TRUE))[1]
                cov_mat[[col_name]] <- relevel(cov_mat[[col_name]], ref = ref_level)
            }
        }
        
        # Remove covariates with zero variance (constant after subsetting)
        covs_to_remove <- c()
        for (col_name in covs_to_include) {
            if (is.numeric(cov_mat[[col_name]])) {
                if (length(unique(cov_mat[[col_name]])) == 1) {
                    covs_to_remove <- c(covs_to_remove, col_name)
                    cat(sprintf("[%s] WARNING: Removing '%s' - constant value (%.2f)\\n", celltype, col_name, unique(cov_mat[[col_name]])[1]))
                }
            } else if (is.factor(cov_mat[[col_name]])) {
                n_levels <- nlevels(cov_mat[[col_name]])
                if (n_levels == 1) {
                    covs_to_remove <- c(covs_to_remove, col_name)
                    cat(sprintf("[%s] WARNING: Removing '%s' - single level (%s)\\n", celltype, col_name, levels(cov_mat[[col_name]])[1]))
                }
            }
        }
        
        # Update covariates list
        if (length(covs_to_remove) > 0) {
            covs_to_include <- setdiff(covs_to_include, covs_to_remove)
            cat(sprintf("[%s] Covariates after removing constant variables: %s\\n", celltype, paste(covs_to_include, collapse=", ")))
        }
        
        # Skip residualization if no valid covariates remain
        if (length(covs_to_include) == 0) {
            cat(sprintf("[%s] WARNING: No valid covariates to regress out. Using normalized expression directly.\\n", celltype))
            exp_mat <- exp_mat %>% as.data.frame() %>% mutate(geneid=row.names(.))
        } else {
            # Calculate residuals for each gene
            lm_formula = paste0("gene ~ ", paste(covs_to_include, collapse = " + "))
            exp_mat <- t(apply(exp_mat, 1, function(gene_exp) {
                lm_data <- data.frame(gene = gene_exp)
                lm_data <- cbind(lm_data, cov_mat[, covs_to_include, drop=FALSE])
                fit <- lm(lm_formula, data = lm_data)
                resid(fit)
            }))
            
            cat(sprintf("[%s] Residuals calculated: %d genes × %d individuals\\n", celltype, nrow(exp_mat), ncol(exp_mat)))

            # Remove zero-variance genes after residualization
            gene_sd <- apply(exp_mat, 1, sd, na.rm = TRUE)
            keep_genes <- which(gene_sd > 0 & !is.na(gene_sd))
            removed_genes <- nrow(exp_mat) - length(keep_genes)
            if (removed_genes > 0) {
                cat(sprintf("[%s] Removing %d zero-variance genes after residuals\\n", celltype, removed_genes))
                exp_mat <- exp_mat[keep_genes, , drop = FALSE]
            }
            
            exp_mat <- exp_mat %>% as.data.frame() %>% mutate(geneid=row.names(.))
        }
    }else{
        exp_mat = exp_mat %>% mutate(geneid=row.names(.))
        cat(sprintf("[%s] No covariate file\\n", celltype))
    }

    data.table::fwrite(exp_mat, paste0(celltype, "_residuals.csv"))
    """

}
