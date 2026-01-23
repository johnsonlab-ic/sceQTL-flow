process merge_seurat_objects {
    tag "merge_${seurat_files.size()}_objects"
    label "process_high_memory"
    publishDir "${params.outdir}/merged_objects/", mode: "copy"

    input: 
    path seurat_files
    path source_R

    output:
    path "merged_seurat.rds", emit: merged_seurat
    
    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    library(dplyr)

    source("$source_R")

    # Get all RDS files
    seurat_files <- list.files(pattern = "\\\\.rds\$", full.names = TRUE)
    
    cat("Found", length(seurat_files), "Seurat objects to merge\\n")
    cat("Files:", paste(seurat_files, collapse = ", "), "\\n")
    
    # Load all Seurat objects
    seurat_list <- lapply(seurat_files, readRDS)
    
    # Validate metadata columns exist in all objects
    for (i in seq_along(seurat_list)) {
        obj <- seurat_list[[i]]
        if (!("${params.celltype_column}" %in% colnames(obj@meta.data))) {
            stop(sprintf("Object %d missing celltype column: ${params.celltype_column}", i))
        }
        if (!("${params.individual_column}" %in% colnames(obj@meta.data))) {
            stop(sprintf("Object %d missing individual column: ${params.individual_column}", i))
        }
        cat(sprintf("Object %d: %d cells\\n", i, ncol(obj)))
    }
    
    # Merge all objects
    cat("Merging Seurat objects...\\n")
    if (length(seurat_list) == 1) {
        merged_seurat <- seurat_list[[1]]
    } else {
        merged_seurat <- merge(
            x = seurat_list[[1]], 
            y = seurat_list[-1],
            merge.data = TRUE
        )
    }
    
    cat("Merged object contains", ncol(merged_seurat), "cells\\n")
    cat("Cell types:", paste(unique(merged_seurat@meta.data[["${params.celltype_column}"]]), collapse = ", "), "\\n")
    cat("Individuals:", paste(unique(merged_seurat@meta.data[["${params.individual_column}"]]), collapse = ", "), "\\n")
    
    # Save merged object
    saveRDS(merged_seurat, "merged_seurat.rds")
    cat("Merge complete!\\n")
    """
}
