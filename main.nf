nextflow.enable.dsl=2

//default inputfiles
params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds"
params.single_cell_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/roche_ms_decontx.rds"
params.cov_file="${baseDir}/R/cov.txt"


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

    If "optimize PCs" is set to TRUE, the pipeline will run longer.

    Generating reports via markdown too now!

    ==============================================

    This is the DEV branch. 

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
    
    // Add get_residuals step with the new covariates parameter
    get_residuals(
        qc_expression.out.pseudobulk_normalised.flatten(),
        params.cov_file,
        params.pseudobulk_source_functions,
        params.covariates_to_include
    )

    // Define the n_pcs channel - always runs now
    n_pcs_ch = Channel.from(1..20)

    // Combine with get_residuals output
    get_residuals.out.residuals_results.flatten()
        .combine(n_pcs_ch)
        .set { combined_ch }

    // Run the optimize_pcs process
    optimize_pcs(
        params.eqtl_source_functions,
        qc_genotype.out.qc_genotype_mat,
        qc_genotype.out.qc_snp_chromlocations,
        combined_ch.map { it[0] },
        pseudobulk_singlecell.out.gene_locations,
        params.cov_file,
        combined_ch.map { it[1] }
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
        
    // Map the residual files by celltype
    get_residuals.out.residuals_results.flatten()
        .map { file ->
            def celltype = file.getBaseName().replace("_residuals", "")
            [celltype, file]
        }
        .set { ch_residual_matrices }

    // For the markdown file
    optimize_pcs.out.egenes_results
        .collect()
        .set { collected_results }

    // Run the select_pcs process
    select_pcs(grouped_results, ch_residual_matrices)
    
    // Join the residual matrices with their corresponding optimal PCs
    ch_residual_matrices
        .join(select_pcs.out.exp_pcs)
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

