// MatrixEQTL workflow (extracted from main.nf)

include { create_genotype; qc_genotype } from '../modules/genotype/genotype.nf'
include { merge_seurat_objects } from '../modules/expression/merge_seurat.nf'
include { pseudobulk_singlecell } from '../modules/expression/pseudobulk.nf'
include { qc_expression } from '../modules/expression/qc_expression.nf'
include { pseudobulk_atac } from '../modules/atac/pseudobulk_atac.nf'
include { qc_atac } from '../modules/atac/qc_atac.nf'
include { preflight_check } from '../modules/qc/preflight_check.nf'
include { subset_samples } from '../modules/qc/subset_samples.nf'
include { get_residuals } from '../modules/residuals/get_residuals.nf'
include { run_matrixeQTL } from '../modules/eqtl/matrixeqtl.nf'
include { combine_eqtls } from '../modules/eqtl/combine_eqtls.nf'
include { final_report } from '../modules/reports/final_report.nf'
include { organize_pc_optimization } from '../modules/reports/organize_pc_optimization.nf'
include { optimize_pcs as optimize_pcs_coarse } from '../modules/optimize/optimize_pcs.nf'
include { optimize_pcs as optimize_pcs_fine } from '../modules/optimize/optimize_pcs.nf'
include { select_pcs } from '../modules/optimize/select_pcs.nf'
include { select_pcs_coarse } from '../modules/optimize/select_pcs_coarse.nf'
include { count_individuals } from '../modules/optimize/count_individuals.nf'
include { generate_fixed_pcs } from '../modules/optimize/generate_fixed_pcs.nf'

workflow matrixeqtl {
    def data_type = params.data_type?.toString()?.toUpperCase() ?: 'RNA'
    if (!(data_type in ['RNA', 'ATAC'])) {
        error "Unknown data_type '${params.data_type}'. Choose 'RNA' or 'ATAC'."
    }

    // Determine which input mode is being used
    def seurat_input_msg = params.single_cell_file_list != "none" && params.single_cell_file_list != "" ? 
        "Multiple Seurat Files (will be merged): ${params.single_cell_file_list}" : 
        "Single Seurat File: ${params.single_cell_file}"
    
    println """
    ========================================
    Welcome to the Nextflow eQTL pipeline
    ========================================
    
    File inputs:

    Output Directory: ${params.outdir}
    GDS File: ${params.gds_file}
    ${seurat_input_msg}
    WorkDir: ${workflow.workDir}
    Using profile: ${workflow.profile}
    Covariate file: ${params.cov_file}

    ==============================================

                RUN PARAMETERS

    Expression/Pseudobulking parameters:

    Min cells for pseudobulking: ${params.min_cells}
    Min percentage for genes: ${params.min_expression}
    Cell-type column: ${params.celltype_column}
    Individual column: ${params.individual_column}
    Assay used: ${params.counts_assay}
    Slot used: ${params.counts_slot}
    Data type: ${data_type}

    eQTL parameters:

    Cis distance: ${params.cis_distance}
    FDR threshold: ${params.fdr_threshold}
    Chromosomes used: ${params.filter_chr}
    Calculating residuals: ${params.cov_file != "none" && params.cov_file != "" ? "YES" : "NO"}
    Optimize PCs: ${params.optimize_pcs ? "YES" : "NO (using " + params.fixed_pcs + " PCs)"}
    PC optimization strategy: coarse_step=${params.pc_coarse_step}, fine_step=${params.pc_fine_step}, fine_window=${params.pc_fine_window}, elbow_tol=${params.pc_elbow_tol}, early_stop_tol=${params.pc_early_stop_tol}, early_stop_patience=${params.pc_early_stop_patience}
    Sample subsetting: ${params.subset_column != "none" && params.subset_values != "none" ? "YES (" + params.subset_column + " = " + params.subset_values + ")" : "NO"}

    If "optimize PCs" is set to TRUE, the pipeline will run longer.

    Generating reports via markdown too now!

    ==============================================

    """
    create_genotype(gds_file= params.gds_file, 
    params.genotype_source_functions)
    qc_genotype(
        genotype_mat=create_genotype.out.genotype_mat,
        snp_chromlocations=create_genotype.out.snp_chromlocations,
        filter_chr=params.filter_chr  // Pass the optional parameter
    )

    // Handle single file vs multiple files
    if (params.single_cell_file_list != "none" && params.single_cell_file_list != "") {
        // Multiple Seurat objects - need to merge first
        seurat_files_ch = nextflow.Channel.fromPath(params.single_cell_file_list.split(',') as List)
        merge_seurat_objects(seurat_files_ch.collect(), params.pseudobulk_source_functions)
        seurat_input = merge_seurat_objects.out.merged_seurat
    } else {
        // Single Seurat object - use directly
        seurat_input = params.single_cell_file
    }

    def pseudobulk_source_functions = data_type == 'ATAC' ? params.atac_source_functions : params.pseudobulk_source_functions

    if (data_type == 'ATAC') {
        pseudobulk_atac(seurat_input, pseudobulk_source_functions)
        pseudobulk_ch = pseudobulk_atac.out.pseudobulk_counts.flatten()
        gene_locations_ch = pseudobulk_atac.out.gene_locations
        qc_atac(pseudobulk_file= pseudobulk_ch)
        qc_output_ch = qc_atac.out.pseudobulk_normalised
    } else {
        pseudobulk_singlecell(seurat_input, pseudobulk_source_functions)
        pseudobulk_ch = pseudobulk_singlecell.out.pseudobulk_counts.flatten()
        gene_locations_ch = pseudobulk_singlecell.out.gene_locations
        qc_expression(pseudobulk_file= pseudobulk_ch)
        qc_output_ch = qc_expression.out.pseudobulk_normalised
    }
    
    // =============================================
    // RUN PREFLIGHT CHECK
    // =============================================
    has_cov = params.cov_file != "none" && params.cov_file != ""
    preflight_cov_file = has_cov ? params.cov_file : "${baseDir}/R/cov.txt"
    preflight_check(
        genotype_mat = qc_genotype.out.qc_genotype_mat,
        cov_file = preflight_cov_file,
        pseudobulk_files = qc_output_ch.collect(),
        has_cov_file = has_cov
    )
    
    // =============================================
    // CONDITIONALLY SUBSET SAMPLES
    // =============================================
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
    
    // Create a channel for expression data - either residuals or normalized data
    
    
    // Conditionally run get_residuals only if a covariate file is provided
    if (params.cov_file != "none" && params.cov_file != "") {
        // Run get_residuals when covariates are provided
        get_residuals(
            qc_output_ch.flatten(),
            cov_file_to_use,
            params.pseudobulk_source_functions,
            params.covariates_to_include
        )
        // Use residuals for downstream analysis
        expression_files_ch = get_residuals.out.residuals_results.flatten()
    } else {
        // Use normalized expression data directly if no covariates
        expression_files_ch = qc_output_ch.flatten()
    }

    // Conditional PC handling: optimize or use fixed number
    if (params.optimize_pcs) {
        // Count individuals in each expression file and determine coarse PC values to test
        count_individuals(expression_files_ch)

        // Create a channel for coarse PC values
        dynamic_pcs_coarse_ch = count_individuals.out.residuals_with_pcs
            .flatMap { file, pc_file ->
                def pc_values = pc_file.text.trim().split('\n')
                pc_values.collect { pc -> [file, pc.toInteger()] }
            }

        // Run the optimize_pcs process on the coarse grid
        optimize_pcs_coarse(
            params.eqtl_source_functions,
            qc_genotype.out.qc_genotype_mat,
            qc_genotype.out.qc_snp_chromlocations,
            dynamic_pcs_coarse_ch.map { it[0] },
            gene_locations_ch,
            dynamic_pcs_coarse_ch.map { it[1] },
            "coarse"
        )

        optimize_pcs_coarse.out.egenes_results
            .map { file ->
                celltype = file.name.replaceAll(/_egenes_vs_.+_coarse\.txt$/, "")
                [celltype, file]
            }
            .groupTuple(by: 0)
            .set { grouped_results_coarse }

        // Map the expression files by celltype
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
            .set { ch_residual_matrices }

        grouped_results_coarse.set { grouped_results_coarse_debug }

        ch_residual_matrices.set { ch_residual_matrices_debug }

        // Select coarse best and generate fine PC grid
        select_pcs_coarse(grouped_results_coarse_debug)

        // Collect coarse summary CSVs for reporting (only file paths)
        select_pcs_coarse.out.coarse_summary
            .map { tuple -> tuple[1] ?: tuple } // If tuple, get file path; else pass as is
            .collect()
            .set { collected_coarse_summaries }

        // Build fine PC values per celltype
        ch_residual_matrices_debug
            .join(select_pcs_coarse.out.fine_pc_values, failOnDuplicate: true, failOnMismatch: true)
            .flatMap { celltype, exp_file, pc_file ->
                def pc_values = pc_file.text.trim().split('\n')
                pc_values.collect { pc -> [exp_file, pc.toInteger()] }
            }
            .set { dynamic_pcs_fine_ch }

        // Run the optimize_pcs process on the fine grid (second invocation with alias)
        optimize_pcs_fine(
            params.eqtl_source_functions,
            qc_genotype.out.qc_genotype_mat,
            qc_genotype.out.qc_snp_chromlocations,
            dynamic_pcs_fine_ch.map { it[0] },
            gene_locations_ch,
            dynamic_pcs_fine_ch.map { it[1] },
            "fine"
        )

        optimize_pcs_fine.out.egenes_results
            .map { file ->
                celltype = file.name.replaceAll(/_egenes_vs_.+_fine\.txt$/, "")
                [celltype, file]
            }
            .groupTuple(by: 0)
            .set { grouped_results_fine }

        // For the markdown file
        optimize_pcs_fine.out.egenes_results
            .collect()
            .set { collected_results }

        // Run the select_pcs process on fine grid results
        grouped_results_fine
            .join(ch_residual_matrices_debug, failOnDuplicate: true, failOnMismatch: true)
            .set { paired_data_fine }

        select_pcs(paired_data_fine.map { celltype, egenes_files, exp_file -> [celltype, egenes_files] },
                  paired_data_fine.map { celltype, egenes_files, exp_file -> [celltype, exp_file] })

        // Collect fine summary CSVs for reporting (only file paths)
        select_pcs.out.fine_summary
            .map { tuple -> tuple[1] ?: tuple } // If tuple, get file path; else pass as is
            .collect()
            .set { collected_fine_summaries }

        // Join the residual matrices with their corresponding optimal PCs
        ch_residual_matrices
            .join(select_pcs.out.exp_pcs, failOnDuplicate: true, failOnMismatch: true)
            .set { residuals_with_pcs }

    } else {
        // Use fixed number of PCs - skip optimization
        println "Skipping PC optimization, using ${params.fixed_pcs} PCs for all cell types"

        // Generate fixed PCs for each expression file
        generate_fixed_pcs(expression_files_ch, params.fixed_pcs)

        // Map expression files by celltype for joining
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
            .set { ch_residual_matrices }

        // Join the residual matrices with their corresponding fixed PCs
        ch_residual_matrices
            .join(generate_fixed_pcs.out.exp_pcs, failOnDuplicate: true, failOnMismatch: true)
            .set { residuals_with_pcs }

        // Provide empty channels for coarse/fine summaries since optimization was skipped
        nextflow.Channel.empty().set { collected_coarse_summaries }
        nextflow.Channel.empty().set { collected_fine_summaries }
    }

    // Run matrixeQTL with PCs (either optimized or fixed)
    run_matrixeQTL(
        params.eqtl_source_functions,
        qc_genotype.out.qc_genotype_mat,
        qc_genotype.out.qc_snp_chromlocations,
        residuals_with_pcs.map { it[1] },  // expression file
        gene_locations_ch,
        residuals_with_pcs.map { it[2] }   // PCs file (optimized or fixed)
    )

    // Collect covariate matrices used per cell type for reporting
    run_matrixeQTL.out.covs_used
        .collect()
        .set { collected_covs_used }

    combine_eqtls(eqtls= run_matrixeQTL.out.eqtl_results.collect())

    if(params.report){
        def unified_report_file = data_type == 'ATAC' ? params.quarto_report_atac : params.quarto_report
        def report_inputs = combine_eqtls.out.mateqtlouts_FDR_filtered
            .combine(combine_eqtls.out.mateqtlouts)
            .combine(nextflow.Channel.value(unified_report_file))
            .combine(collected_coarse_summaries.map { files -> [files] })
            .combine(collected_fine_summaries.map { files -> [files] })
            .combine(collected_covs_used.map { files -> [files] })
        report_inputs | final_report
    }
}

workflow.onComplete {
    def active_wf = params.workflow?.toLowerCase() ?: 'matrixeqtl'
    if (active_wf != 'matrixeqtl') {
        return
    }

    println """
    ========================================
    Pipeline Completed!
    ========================================

    eQTL report generated at: ${params.outdir}/eQTL_outputs/eqtl_report.html

    """
}
