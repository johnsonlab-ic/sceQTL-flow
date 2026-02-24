process final_report {
    label "process_high"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    tuple path(eqtl_results_filtered), path(eqtl_results), path(report_file), path(coarse_summaries, stageAs: "coarse/*"), path(fine_summaries, stageAs: "fine/*"), path(covs_used, stageAs: "covs/*")

    output:
    path "eqtl_report.html"

    script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(input = "$report_file", output_file = "eqtl_report.html", params = list(
        eqtl_results_filtered = "$eqtl_results_filtered",
        eqtl_results = "$eqtl_results",
        coarse_summaries = list.files("coarse", full.names = TRUE),
        fine_summaries = list.files("fine", full.names = TRUE),
        covs_used = list.files("covs", full.names = TRUE),
        outdir = "${params.outdir}",
        gds_file = "${params.gds_file}",
        single_cell_file = "${params.single_cell_file}",
        single_cell_file_list = "${params.single_cell_file_list}",
        counts_assay = "${params.counts_assay}",
        counts_slot = "${params.counts_slot}",
        celltype_column = "${params.celltype_column}",
        individual_column = "${params.individual_column}",
        min_cells = ${params.min_cells},
        min_expression = ${params.min_expression},
        cis_distance = ${params.cis_distance},
        filter_chr = "${params.filter_chr}",
        fdr_threshold = ${params.fdr_threshold},
        cov_file = "${params.cov_file}",
        covariates_to_include = "${params.covariates_to_include}",
        subset_column = "${params.subset_column}",
        subset_values = "${params.subset_values}",
        optimize_pcs = ${params.optimize_pcs ? 'TRUE' : 'FALSE'},
        fixed_pcs = ${params.fixed_pcs},
        pc_coarse_step = ${params.pc_coarse_step},
        pc_fine_step = ${params.pc_fine_step},
        pc_fine_window = ${params.pc_fine_window},
        pc_elbow_tol = ${params.pc_elbow_tol},
        pc_early_stop_tol = ${params.pc_early_stop_tol},
        pc_early_stop_patience = ${params.pc_early_stop_patience},
        workflow = "${params.workflow}",
        profile = "${workflow.profile}"
    ))
    """
}
