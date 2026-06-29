process check_overlap {
    tag "genotype <-> single-cell overlap"
    label "process_single"
    publishDir "${params.outdir}/QC/", mode: "copy"

    input:
    path genotype_mat
    path cell_metadata
    path script

    output:
    path "overlap_report.txt", emit: report

    script:
    """
    set -o pipefail
    python3 ${script} \\
        --genotype-matrix ${genotype_mat} \\
        --cell-metadata ${cell_metadata} \\
        --indiv-col ${params.individual_column} \\
        --warn-frac ${params.overlap_warn_frac} | tee overlap_report.txt
    """
}
