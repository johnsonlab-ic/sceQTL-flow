nextflow.enable.dsl=2

//default inputfiles
params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds"
params.single_cell_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/roche_ms_decontx.rds"
params.cov_file="none"


// source functions for easy troubleshooting
params.genotype_source_functions="${baseDir}/R/genotype_functions/genotype_functions.r"
params.pseudobulk_source_functions="${baseDir}/R/expression_functions/pseudobulk_functions.r"
params.eqtl_source_functions="${baseDir}/R/MatrixEQTL_functions/matrixeqtl_source.r"
params.quarto_report="${baseDir}/R/rmarkdown_reports/final_report.Rmd"

params.min_cells=5
params.min_expression=0.1
params.celltype_column="celltype"
params.individual_column="individual"
params.counts_assay="RNA"
params.counts_slot="counts"
params.cis_distance=1000000
params.fdr_threshold=0.05
params.filter_chr = "all" // Optional parameter for filtering by chromosome. use "chr6"


include { create_genotype } from './NEXTFLOW/genotype.nf'
include { qc_genotype } from './NEXTFLOW/genotype.nf'
include { pseudobulk_singlecell } from './NEXTFLOW/pseudobulk.nf'
include { qc_expression } from './NEXTFLOW/qc_expression.nf'
include { get_residuals } from './NEXTFLOW/get_residuals.nf'
include { run_matrixeQTL } from './NEXTFLOW/matrixeqtl.nf'
include { combine_eqtls } from './NEXTFLOW/combine_eqtls.nf'
include { final_report } from './NEXTFLOW/final_report.nf'
include { optimize_pcs } from './NEXTFLOW/optimize/optimize_pcs.nf'
include { select_pcs } from './NEXTFLOW/optimize/select_pcs.nf'
include { count_individuals } from './NEXTFLOW/optimize/count_individuals.nf'


// default parameters 
// Add a new parameter for specifying which covariates to include
params.covariates_to_include = "all" // Default to include all covariates

workflow {
    println """
    ========================================
    Welcome to the Nextflow eQTL pipeline
    ========================================
    
    File inputs:

    Output Directory: ${params.outdir}
    GDS File: ${params.gds_file}
    Input Seurat File: ${params.single_cell_file}
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

    pseudobulk_singlecell(params.single_cell_file, params.pseudobulk_source_functions)
    pseudobulk_ch = pseudobulk_singlecell.out.pseudobulk_counts.flatten()
    qc_expression(pseudobulk_file= pseudobulk_ch)
    
    // Create a channel for expression data - either residuals or normalized data
    
    
    // Conditionally run get_residuals only if a covariate file is provided
    if (params.cov_file != "none" && params.cov_file != "") {
        // Run get_residuals when covariates are provided
        get_residuals(
            qc_expression.out.pseudobulk_normalised.flatten(),
            params.cov_file,
            params.pseudobulk_source_functions,
            params.covariates_to_include
        )
        // Use residuals for downstream analysis
        expression_files_ch = get_residuals.out.residuals_results.flatten()
    } else {
        // Use normalized expression data directly if no covariates
        expression_files_ch = qc_expression.out.pseudobulk_normalised.flatten()
    }

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
            // Extract cell type from the file name
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
            if (file.name.contains("_residuals")) {
                celltype = file.getBaseName().replace("_residuals", "")
            } else {
                celltype = file.getBaseName().replace("_pseudobulk_normalised", "")
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

    // Run matrixeQTL with optimized PCs for each cell type
    run_matrixeQTL(
        params.eqtl_source_functions,
        qc_genotype.out.qc_genotype_mat,
        qc_genotype.out.qc_snp_chromlocations,
        residuals_with_pcs.map { it[1] },  // expression file
        pseudobulk_singlecell.out.gene_locations,
        residuals_with_pcs.map { it[2] }   // optimized PCs file
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

