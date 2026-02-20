process qc_atac {
    tag "${pseudobulk_file}"
    label "process_single"
    publishDir "${params.outdir}/expression_matrices/", mode: 'copy'

    input:
    path pseudobulk_file

    output:
    path "*pseudobulk_normalised.csv", emit: pseudobulk_normalised

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(dplyr)
    library(edgeR)

    pseudobulk_data <- fread("$pseudobulk_file")
    pseudobulk_data <- pseudobulk_data %>% tibble::column_to_rownames(var="geneid")

    # Log initial counts
    n_peaks_initial <- nrow(pseudobulk_data)
    n_individuals_initial <- ncol(pseudobulk_data)
    cat(sprintf("[QC] Initial: %d peaks × %d individuals\n", n_peaks_initial, n_individuals_initial))

    # Filter peaks by expression threshold
    min_percentage <- as.numeric(${params.min_expression})
    min_individuals <- min_percentage * ncol(pseudobulk_data)
    pseudobulk_data <- pseudobulk_data[rowSums(pseudobulk_data > 0) >= min_individuals, ]

    # Log filtered counts
    n_peaks_final <- nrow(pseudobulk_data)
    n_peaks_removed <- n_peaks_initial - n_peaks_final
    n_individuals_final <- ncol(pseudobulk_data)
    cat(sprintf("[QC] Removed: %d peaks (%.1f%%) | Kept: %d individuals (threshold: ≥%.1f%% expressing)\n",
                n_peaks_removed, 100*n_peaks_removed/n_peaks_initial, n_individuals_final, 100*min_percentage))
    cat(sprintf("[QC] Final: %d peaks × %d individuals\n", n_peaks_final, n_individuals_final))

    # edgeR TMM normalization
    dge <- edgeR::DGEList(counts = pseudobulk_data)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    norm_mat <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

    cell_type_name <- gsub("_pseudobulk.csv", "", basename("$pseudobulk_file"))
    norm_mat <- as.data.frame(norm_mat) %>% mutate(geneid = row.names(.))
    fwrite(norm_mat, paste0(cell_type_name, "_pseudobulk_normalised.csv"))
    """
}