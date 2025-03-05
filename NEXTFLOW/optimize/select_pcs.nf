process select_pcs {

    label "process_low_memory"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path combined_results

    output:
    path "optimal_pcs_per_celltype.txt", emit: optimal_pcs

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)

    # Read the combined results
    results <- fread("$combined_results")

    # Find the optimal n_pcs for each celltype
    optimal_pcs <- results[, .SD[which.max(n_genes)], by = celltype]

    # Write the optimal n_pcs to a file
    fwrite(optimal_pcs, "optimal_pcs_per_celltype.txt")

    """
}
