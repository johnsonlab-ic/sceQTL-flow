process pseudobulk_seurat {
    tag "${tag}"
    label "process_high_memory"
    publishDir "${params.outdir}/expression_matrices/partials/", mode: "copy"

    input:
    tuple val(tag), path(rds)
    path script
    path cell_metadata

    output:
    path "*_partial.csv", emit: partials
    path "*_ncells.csv",  emit: ncells

    script:
    def meta_arg = cell_metadata.name != 'NO_FILE' ?
        "--cell-metadata ${cell_metadata} --id-col ${params.metadata_id_col}" : ""
    """
    Rscript ${script} \\
        --rds ${rds} \\
        --celltype-col ${params.celltype_column} \\
        --indiv-col ${params.individual_column} \\
        --counts-assay ${params.counts_assay} \\
        --counts-slot ${params.counts_slot} \\
        ${meta_arg} \\
        --outdir . --tag ${tag}
    """
}
