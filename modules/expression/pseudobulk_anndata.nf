process pseudobulk_anndata {
    tag "${tag}"
    label "process_high_memory"
    publishDir "${params.outdir}/expression_matrices/partials/", mode: "copy"

    input:
    tuple val(tag), path(h5ad)
    path script
    path cell_metadata

    output:
    path "*_partial.csv", emit: partials
    path "*_ncells.csv",  emit: ncells

    script:
    def meta_arg = cell_metadata.name != 'NO_FILE' ?
        "--cell-metadata ${cell_metadata} --id-col ${params.metadata_id_col}" : ""
    """
    python3 ${script} \\
        --h5ad ${h5ad} \\
        --celltype-col ${params.celltype_column} \\
        --indiv-col ${params.individual_column} \\
        --counts-layer ${params.counts_slot} \\
        ${meta_arg} \\
        --outdir . --tag ${tag}
    """
}
