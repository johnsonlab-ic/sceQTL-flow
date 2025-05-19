process select_pcs {

    label "process_low_memory"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    tuple val(celltype), path(egenes_files)
    tuple val(celltype), path(exp_matrix)

    output:
    tuple val(celltype), path("optimal_pcs_${celltype}.txt"), emit: exp_pcs

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(ggplot2)
    library(dplyr)

    # Combine all egenes files into a single data frame
    file_list <- unlist(strsplit("${egenes_files}", " "))
    results <- rbindlist(lapply(file_list, fread))

    # Find the optimal n_pcs for the cell type
    optimal_pcs <- results[which.max(n_egenes), ]
    n_pcs <- optimal_pcs\$n_pcs

    # Perform PCA on the expression matrix
    exp_mat <- fread("${exp_matrix}") %>% tibble::column_to_rownames(var="geneid")
    exp_pcs <- prcomp(t(exp_mat), scale. = TRUE)
    exp_pcs <- exp_pcs\$x[, 1:n_pcs]
    exp_pcs <- as.data.frame(exp_pcs)
    colnames(exp_pcs) <- paste0("PC", 1:n_pcs)
    exp_pcs <- t(exp_pcs)

    # Write the PC covariate matrix to a file
    write.table(exp_pcs, file="optimal_pcs_${celltype}.txt", quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
    
    # Also save information about how many PCs were chosen
    writeLines(paste("Cell type:", "${celltype}", "\nOptimal number of PCs:", n_pcs), "pc_info_${celltype}.txt")
    """
}