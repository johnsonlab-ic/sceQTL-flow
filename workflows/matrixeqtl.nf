// MatrixEQTL workflow (extracted from main.nf)

include { create_genotype; qc_genotype } from '../modules/genotype/genotype.nf'
include { merge_seurat_objects } from '../modules/expression/merge_seurat.nf'
include { pseudobulk_singlecell } from '../modules/expression/pseudobulk.nf'
include { qc_expression } from '../modules/expression/qc_expression.nf'
include { preflight_check } from '../modules/qc/preflight_check.nf'
include { subset_samples } from '../modules/qc/subset_samples.nf'
include { get_residuals } from '../modules/residuals/get_residuals.nf'
include { run_matrixeQTL } from '../modules/eqtl/matrixeqtl.nf'
include { combine_eqtls } from '../modules/eqtl/combine_eqtls.nf'
include { final_report } from '../modules/reports/final_report.nf'
include { optimize_pcs } from '../modules/optimize/optimize_pcs.nf'
include { select_pcs } from '../modules/optimize/select_pcs.nf'
include { count_individuals } from '../modules/optimize/count_individuals.nf'
include { generate_fixed_pcs } from '../modules/optimize/generate_fixed_pcs.nf'

workflow matrixeqtl {
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

    eQTL parameters:

    Cis distance: ${params.cis_distance}
    FDR threshold: ${params.fdr_threshold}
    Chromosomes used: ${params.filter_chr}
    Calculating residuals: ${params.cov_file != "none" && params.cov_file != "" ? "YES" : "NO"}
    Optimize PCs: ${params.optimize_pcs ? "YES" : "NO (using " + params.fixed_pcs + " PCs)"}
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
        seurat_files_ch = Channel.fromPath(params.single_cell_file_list.split(',') as List)
        merge_seurat_objects(seurat_files_ch.collect(), params.pseudobulk_source_functions)
        seurat_input = merge_seurat_objects.out.merged_seurat
    } else {
        // Single Seurat object - use directly
        seurat_input = params.single_cell_file
    }

    pseudobulk_singlecell(seurat_input, params.pseudobulk_source_functions)
    pseudobulk_ch = pseudobulk_singlecell.out.pseudobulk_counts.flatten()
    qc_expression(pseudobulk_file= pseudobulk_ch)
    
    // =============================================
    // RUN PREFLIGHT CHECK
    // =============================================
    has_cov = params.cov_file != "none" && params.cov_file != ""
    preflight_check(
        genotype_mat = qc_genotype.out.qc_genotype_mat,
        cov_file = has_cov ? params.cov_file : "",
        pseudobulk_files = qc_expression.out.pseudobulk_normalised.collect(),
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
            qc_expression.out.pseudobulk_normalised.flatten(),
            cov_file_to_use,
            params.pseudobulk_source_functions,
            params.covariates_to_include
        )
        // Use residuals for downstream analysis
        expression_files_ch = get_residuals.out.residuals_results.flatten()
    } else {
        // Use normalized expression data directly if no covariates
        expression_files_ch = qc_expression.out.pseudobulk_normalised.flatten()
    }

    // Conditional PC handling: optimize or use fixed number
    if (params.optimize_pcs) {
        // Count individuals in each expression file and determine PC values to test
        count_individuals(expression_files_ch)
        
        // Create a channel for dynamic PC values
        dynamic_pcs_ch = count_individuals.out.residuals_with_pcs
            .flatMap { file, pc_file -> 
                // Read PC values from the file
                def pc_values = pc_file.text.trim().split('\n')
                pc_values.collect { pc -> [file, pc.toInteger()] }
            }
            
        // Run the optimize_pcs process with dynamic PC values based on sample count
        optimize_pcs(
            params.eqtl_source_functions,
            qc_genotype.out.qc_genotype_mat,
            qc_genotype.out.qc_snp_chromlocations,
            dynamic_pcs_ch.map { it[0] },  // expression file
            pseudobulk_singlecell.out.gene_locations,
            dynamic_pcs_ch.map { it[1] }   // pc value
        )

        // Collect all egenes_results into a single list
        optimize_pcs.out.egenes_results
            .map { file -> 
                // Extract cell type from the file name - remove the egenes_vs_XX.txt suffix
                celltype = file.name.replaceAll(/_egenes_vs_.+\.txt$/, "")
                [celltype, file]
            }
            .groupTuple(by: 0)  // Group by cell type
            .set { grouped_results }
            
        // Map the expression files by celltype
        expression_files_ch
            .map { file ->
                // Handle different file naming patterns based on whether it's from residuals or normalized expression
                def celltype
                if (file.name.contains("_residuals.csv")) {
                    celltype = file.name.replace("_residuals.csv", "")
                } else if (file.name.contains("_pseudobulk_normalised.csv")) {
                    celltype = file.name.replace("_pseudobulk_normalised.csv", "")
                } else {
                    // Fallback: use basename without extension
                    celltype = file.getBaseName()
                }
                [celltype, file]
            }
            .set { ch_residual_matrices }

        // For the markdown file
        optimize_pcs.out.egenes_results
            .collect()
            .set { collected_results }

        // Debug: Print the cell types in each channel
        grouped_results.map { celltype, files -> 
            println "Grouped results for cell type: $celltype with ${files.size()} files"
            return [celltype, files]
        }.set { grouped_results_debug }
        
        ch_residual_matrices.map { celltype, file -> 
            println "Expression file for cell type: $celltype - ${file.name}"
            return [celltype, file]
        }.set { ch_residual_matrices_debug }
        
        // Run the select_pcs process
        grouped_results_debug
            .join(ch_residual_matrices_debug, failOnDuplicate: true, failOnMismatch: true)
            .map { celltype, egenes_files, exp_file ->
                println "Successfully paired: $celltype with ${egenes_files.size()} egenes files and ${exp_file.name}"
                return [celltype, egenes_files, exp_file]
            }
            .set { paired_data }
        
        select_pcs(paired_data.map { celltype, egenes_files, exp_file -> [celltype, egenes_files] }, 
                  paired_data.map { celltype, egenes_files, exp_file -> [celltype, exp_file] })
        
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
            
        // Provide a placeholder value so downstream reporting still runs
        Channel.value("").set { collected_results }
    }

    // Run matrixeQTL with PCs (either optimized or fixed)
    run_matrixeQTL(
        params.eqtl_source_functions,
        qc_genotype.out.qc_genotype_mat,
        qc_genotype.out.qc_snp_chromlocations,
        residuals_with_pcs.map { it[1] },  // expression file
        pseudobulk_singlecell.out.gene_locations,
        residuals_with_pcs.map { it[2] }   // PCs file (optimized or fixed)
    )

    combine_eqtls(eqtls= run_matrixeQTL.out.eqtl_results.collect())

    if(params.report){
        final_report(
            eqtl_results_filtered = combine_eqtls.out.mateqtlouts_FDR_filtered,
            eqtl_results = combine_eqtls.out.mateqtlouts,
            report_file = params.quarto_report,
            optimization_results = collected_results
        )
    }
}

workflow.onComplete {
    println """
    ========================================
    Pipeline Completed!
    ========================================

    eQTL report generated at: ${params.outdir}/eQTL_outputs/eqtl_report.html

    """
}
