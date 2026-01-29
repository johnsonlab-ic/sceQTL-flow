process select_pcs {
    tag "Selecting optimal PCs for ${celltype}"
    label "process_low_memory"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    tuple val(celltype), path(egenes_files)
    tuple val(celltype), path(exp_matrix)

    output:
    tuple val(celltype), path("*_pcs.txt"), emit: exp_pcs

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(ggplot2)
    library(dplyr)

    # Combine all egenes files into a single data frame
    file_list <- unlist(strsplit("${egenes_files}", " "))
    results <- rbindlist(lapply(file_list, fread))

    # Find the optimal n_pcs using elbow tolerance (smallest n_pcs within tol of max)
    results\$n_pcs <- as.integer(results\$n_pcs)
    results\$n_assoc <- as.numeric(results\$n_assoc)
    results <- results[order(results\$n_pcs)]

    max_assoc <- max(results\$n_assoc, na.rm = TRUE)
    threshold <- max_assoc * (1 - ${params.pc_elbow_tol})
    candidates <- results[results\$n_assoc >= threshold, ]

    if (nrow(candidates) == 0) {
        optimal_idx <- which.max(results\$n_assoc)
        n_pcs <- results[optimal_idx, ]\$n_pcs
    } else {
        n_pcs <- candidates\$n_pcs[1]
    }

    cat("Selected n_pcs:", n_pcs, "max_n_assoc:", max_assoc, "threshold:", threshold, "\n")

    # Perform PCA on the expression matrix
    exp_mat <- fread("${exp_matrix}") %>% tibble::column_to_rownames(var="geneid")
    exp_pcs <- prcomp(t(exp_mat), scale. = TRUE)
    exp_pcs <- exp_pcs\$x[, 1:n_pcs]
    exp_pcs <- as.data.frame(exp_pcs)
    colnames(exp_pcs) <- paste0("PC", 1:n_pcs)
    exp_pcs <- t(exp_pcs)

    # Write the PC covariate matrix to a file
    write.table(exp_pcs, file=paste0("${celltype}_",n_pcs,"_pcs.txt"), quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
    
    # Also save information about how many PCs were chosen
    writeLines(paste("Cell type:", "${celltype}", "\nOptimal number of PCs:", n_pcs), "pc_info_${celltype}.txt")
    """
}
