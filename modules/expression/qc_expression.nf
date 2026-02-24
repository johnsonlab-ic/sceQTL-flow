process qc_expression {
    tag "${pseudobulk_file}"
    label "process_single"
    publishDir "${params.outdir}/expression_matrices/", mode: 'copy'

    input:

    path pseudobulk_file

    output:

    path "*pseudobulk_normalised.csv" , emit: pseudobulk_normalised

    script:

    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(dplyr)
    pseudobulk_data <- fread("$pseudobulk_file")
    pseudobulk_data <- pseudobulk_data %>% tibble::column_to_rownames(var="geneid")
    
    # Log initial counts
    n_genes_initial <- nrow(pseudobulk_data)
    n_individuals_initial <- ncol(pseudobulk_data)
    cat(sprintf("[QC] Initial: %d genes × %d individuals\\n", n_genes_initial, n_individuals_initial))
    
    # Filter genes by expression threshold
    min_percentage <- as.numeric(${params.min_expression})
    min_individuals <- min_percentage * ncol(pseudobulk_data)
    pseudobulk_data <- pseudobulk_data[rowSums(pseudobulk_data > 0) >= min_individuals, ]
    
    # Log filtered counts
    n_genes_final <- nrow(pseudobulk_data)
    n_genes_removed <- n_genes_initial - n_genes_final
    n_individuals_final <- ncol(pseudobulk_data)
    cat(sprintf("[QC] Removed: %d genes (%.1f%%) | Kept: %d individuals (threshold: ≥%.1f%% expressing)\\n", 
                n_genes_removed, 100*n_genes_removed/n_genes_initial, n_individuals_final, 100*min_percentage))
    cat(sprintf("[QC] Final: %d genes × %d individuals\\n", n_genes_final, n_individuals_final))
    
    cell_type_name <- gsub("_pseudobulk.csv", "", "$pseudobulk_file")
    pseudobulk_data=log2(edgeR::cpm(pseudobulk_data)+1) %>% as.data.frame()
    pseudobulk_data = pseudobulk_data %>% mutate(geneid=row.names(.))
    fwrite(pseudobulk_data, paste0(cell_type_name, "_pseudobulk_normalised.csv"))
    """
}
