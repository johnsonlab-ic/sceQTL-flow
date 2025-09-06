process generate_fixed_pcs {
    tag "${expression_mat}"
    label "process_low"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path expression_mat
    val n_pcs

    output:
    tuple val(celltype), path("*_pcs.txt"), emit: exp_pcs

    script:
    // Extract celltype from filename in Nextflow
    celltype = expression_mat.name.contains("_residuals.csv") ? 
               expression_mat.name.replace("_residuals.csv", "") : 
               expression_mat.name.replace("_pseudobulk_normalised.csv", "")
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(dplyr)

    # Read expression matrix
    exp_mat <- fread("${expression_mat}") %>% tibble::column_to_rownames(var="geneid")
    
    # Use celltype passed from Nextflow
    celltype <- "${celltype}"
    
    # Perform PCA and extract fixed number of PCs
    exp_pcs <- prcomp(t(exp_mat), scale. = TRUE)
    n_pcs_to_use <- min(${n_pcs}, ncol(exp_pcs\$x))  # Don't exceed available PCs
    exp_pcs <- exp_pcs\$x[, 1:n_pcs_to_use]
    exp_pcs <- as.data.frame(exp_pcs)
    colnames(exp_pcs) <- paste0("PC", 1:n_pcs_to_use)
    exp_pcs <- t(exp_pcs)

    # Write the PC covariate matrix to a file
    write.table(exp_pcs, file=paste0(celltype, "_", n_pcs_to_use, "_pcs.txt"), 
                quote=FALSE, sep="\\t", col.names=TRUE, row.names=TRUE)
    
    # Save information about how many PCs were used
    writeLines(paste("Cell type:", celltype, "\\nFixed number of PCs:", n_pcs_to_use), 
               paste0("pc_info_", celltype, ".txt"))
    
    cat("Generated", n_pcs_to_use, "PCs for", celltype, "\\n")
    """
}
