include { create_genotype; qc_genotype } from '../modules/genotype/genotype.nf'
include { merge_seurat_objects } from '../modules/expression/merge_seurat.nf'
include { pseudobulk_singlecell } from '../modules/expression/pseudobulk.nf'
include { qc_expression } from '../modules/expression/qc_expression.nf'
include { pseudobulk_atac } from '../modules/atac/pseudobulk_atac.nf'
include { qc_atac } from '../modules/atac/qc_atac.nf'
include { preflight_check } from '../modules/qc/preflight_check.nf'
include { subset_samples } from '../modules/qc/subset_samples.nf'
include { get_residuals } from '../modules/residuals/get_residuals.nf'
include { generate_fixed_pcs } from '../modules/optimize/generate_fixed_pcs.nf'
include { run_tensorqtl } from '../modules/eqtl/tensorqtl_cis.nf'

workflow tensorqtl {
    def data_type = params.data_type?.toString()?.toUpperCase() ?: 'RNA'
    if (!(data_type in ['RNA', 'ATAC'])) {
        error "Unknown data_type '${params.data_type}'. Choose 'RNA' or 'ATAC'."
    }

    if (params.optimize_pcs) {
        log.warn "tensorQTL workflow currently supports fixed PCs only. Ignoring --optimize_pcs true and using --fixed_pcs ${params.fixed_pcs}."
    }

    def seurat_input_msg = params.single_cell_file_list != "none" && params.single_cell_file_list != "" ?
        "Multiple Seurat Files (will be merged): ${params.single_cell_file_list}" :
        "Single Seurat File: ${params.single_cell_file}"

    println """
    ========================================
    Running tensorQTL workflow (fixed PCs)
    ========================================

    Output Directory: ${params.outdir}
    GDS File: ${params.gds_file}
    ${seurat_input_msg}
    WorkDir: ${workflow.workDir}
    Using profile: ${workflow.profile}
    Covariate file: ${params.cov_file}
    Data type: ${data_type}

    Fixed PCs: ${params.fixed_pcs}
    Cis distance: ${params.cis_distance}
    Chromosomes used: ${params.filter_chr}
    Reporting: DISABLED for tensorQTL currently
    PC optimization: DISABLED for tensorQTL currently

    ========================================
    """

    create_genotype(
        gds_file = params.gds_file,
        params.genotype_source_functions
    )

    qc_genotype(
        genotype_mat = create_genotype.out.genotype_mat,
        snp_chromlocations = create_genotype.out.snp_chromlocations,
        filter_chr = params.filter_chr
    )

    if (params.single_cell_file_list != "none" && params.single_cell_file_list != "") {
        seurat_files_ch = nextflow.Channel.fromPath(params.single_cell_file_list.split(',') as List)
        merge_seurat_objects(seurat_files_ch.collect(), params.pseudobulk_source_functions)
        seurat_input = merge_seurat_objects.out.merged_seurat
    } else {
        seurat_input = params.single_cell_file
    }

    def pseudobulk_source_functions = data_type == 'ATAC' ? params.atac_source_functions : params.pseudobulk_source_functions

    if (data_type == 'ATAC') {
        pseudobulk_atac(seurat_input, pseudobulk_source_functions)
        pseudobulk_ch = pseudobulk_atac.out.pseudobulk_counts.flatten()
        gene_locations_ch = pseudobulk_atac.out.gene_locations
        qc_atac(pseudobulk_file = pseudobulk_ch)
        qc_output_ch = qc_atac.out.pseudobulk_normalised
    } else {
        pseudobulk_singlecell(seurat_input, pseudobulk_source_functions)
        pseudobulk_ch = pseudobulk_singlecell.out.pseudobulk_counts.flatten()
        gene_locations_ch = pseudobulk_singlecell.out.gene_locations
        qc_expression(pseudobulk_file = pseudobulk_ch)
        qc_output_ch = qc_expression.out.pseudobulk_normalised
    }

    has_cov = params.cov_file != "none" && params.cov_file != ""
    preflight_cov_file = has_cov ? params.cov_file : "${baseDir}/R/cov.txt"
    preflight_check(
        genotype_mat = qc_genotype.out.qc_genotype_mat,
        cov_file = preflight_cov_file,
        pseudobulk_files = qc_output_ch.collect(),
        has_cov_file = has_cov
    )

    def cov_file_to_use = params.cov_file
    if (params.subset_column != "none" && params.subset_values != "none" && params.cov_file != "none" && params.cov_file != "") {
        subset_samples(
            params.cov_file,
            params.subset_column,
            params.subset_values
        )
        cov_file_to_use = subset_samples.out.filtered_cov
        log.info "Subsetting samples by ${params.subset_column} = ${params.subset_values}"
    }

    if (params.cov_file != "none" && params.cov_file != "") {
        get_residuals(
            qc_output_ch.flatten(),
            cov_file_to_use,
            params.pseudobulk_source_functions,
            params.covariates_to_include
        )
        expression_files_ch = get_residuals.out.residuals_results.flatten()
    } else {
        expression_files_ch = qc_output_ch.flatten()
    }

    generate_fixed_pcs(expression_files_ch, params.fixed_pcs)

    expression_files_ch
        .map { file ->
            def celltype
            if (file.name.contains("_residuals.csv")) {
                celltype = file.name.replace("_residuals.csv", "")
            } else if (file.name.contains("_pseudobulk_normalised.csv")) {
                celltype = file.name.replace("_pseudobulk_normalised.csv", "")
            } else {
                celltype = file.getBaseName()
            }
            [celltype, file]
        }
        .set { ch_expression_matrices }

    ch_expression_matrices
        .join(generate_fixed_pcs.out.exp_pcs, failOnDuplicate: true, failOnMismatch: true)
        .set { expression_with_pcs }

    run_tensorqtl(
        qc_genotype.out.qc_genotype_mat,
        qc_genotype.out.qc_snp_chromlocations,
        expression_with_pcs.map { it[1] },
        gene_locations_ch,
        expression_with_pcs.map { it[2] }
    )
}

workflow.onComplete {
    def active_wf = params.workflow?.toLowerCase() ?: 'matrixeqtl'
    if (active_wf != 'tensorqtl') {
        return
    }

    println """
    ========================================
    tensorQTL workflow completed.
    ========================================

    Outputs are in: ${params.outdir}/eQTL_outputs/
    (parquet files: *_tensorqtl.cis_qtl_pairs.<chr>.parquet)

    """
}
