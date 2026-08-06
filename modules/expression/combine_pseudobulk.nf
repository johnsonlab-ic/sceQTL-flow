process combine_pseudobulk {
    label "process_medium"
    publishDir "${params.outdir}/expression_matrices/", mode: "copy"

    input:
    path partials
    path ncells
    path source_R
    path combine_script

    output:
    path "*_pseudobulk.csv", emit: pseudobulk_counts
    path "gene_locations.csv", emit: gene_locations
    path "cells_per_individual.csv", emit: cells_per_individual

    script:
    """
    Rscript ${combine_script} \\
        --indir . --outdir . \\
        --min-cells ${params.min_cells} \\
        --source ${source_R}
    """
}
