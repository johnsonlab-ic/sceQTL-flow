process final_report {
    label "process_high_memory"
    publishDir "${params.outdir}", mode: 'copy'

    input: 
    path eqtl_results_filtered
    path eqtl_results
    path report_file
    path optimization_results

    output: 
    path "eqtl_report.html"

    script:
    """
    #!/usr/bin/env Rscript
     if (length("$optimization_results") == 0) {
        optimization_results <- NULL
    } else {
        optimization_results <- "$optimization_results"
    }
    rmarkdown::render(input = "$report_file", output_file = "eqtl_report.html", params = list(
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
        filter_chr="${params.filter_chr}",
        fdr_threshold = ${params.fdr_threshold},
        optimize_pcs = TRUE,
        optimization_results = "$optimization_results"
    ))
    """
}
