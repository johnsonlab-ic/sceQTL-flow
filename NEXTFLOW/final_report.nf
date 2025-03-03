process final_report {
    label "process_single"
    publishDir "${params.outdir}", mode: 'copy'
    input: 
    path eqtl_results_filtered
    path eqtl_results
    path report_file
    output: 
    path "report.html"
    script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(input = "$report_file", output_file = "report.html", params = list(
        eqtl_results_filtered = "$eqtl_results_filtered",
        eqtl_results = "$eqtl_results",
        outdir = "${params.outdir}",
        gds_file = "${params.gds_file}",
        single_cell_file = "${params.single_cell_file}",
        counts_assay = "${params.counts_assay}",
        counts_slot = "${params.counts_slot}",
        celltype_column = "${params.celltype_column}",
        individual_column = "${params.individual_column}",
        min_cells = ${params.min_cells},
        min_expression = ${params.min_expression},
        cis_distance = ${params.cis_distance},
        fdr_threshold = ${params.fdr_threshold},
        optimize_pcs = ${params.optimize_pcs ? 'TRUE' : 'FALSE'}
    ))
    """
}
