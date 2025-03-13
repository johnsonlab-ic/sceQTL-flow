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
params.optimize_pcs = true

include { create_genotype } from './NEXTFLOW/genotype.nf'
include { qc_genotype } from './NEXTFLOW/genotype.nf'
include { pseudobulk_singlecell } from './NEXTFLOW/pseudobulk.nf'
include { qc_expression } from './NEXTFLOW/qc_expression.nf'
include { run_matrixeQTL } from './NEXTFLOW/matrixeqtl.nf'
include { combine_eqtls } from './NEXTFLOW/combine_eqtls.nf'
include { final_report } from './NEXTFLOW/final_report.nf'
include { optimize_pcs } from './NEXTFLOW/optimize/optimize_pcs.nf'
include { select_pcs } from './NEXTFLOW/optimize/select_pcs.nf'

// default parameters 



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
    Optimize PCs: ${params.optimize_pcs}
    Chromosomes used: ${params.filter_chr}

    If "optimize PCs" is set to TRUE, the pipeline will run longer.

    Generating reports via markdown too now!

    ==============================================

    This is the stable release.

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

    if (params.optimize_pcs) {
        // Define the n_pcs channel
        n_pcs_ch = Channel.from(1..5)

        // Combine pseudobulk_ch with n_pcs_ch
        qc_expression.out.pseudobulk_normalised.flatten()
            .combine(n_pcs_ch)
            .set { combined_ch }

        // Run the optimize_pcs process
        optimize_pcs(
            params.eqtl_source_functions,  // source_R
            qc_genotype.out.qc_genotype_mat,  // genotype_mat
            qc_genotype.out.qc_snp_chromlocations,  // snp_locations
            combined_ch.map { it[0] },  // expression_mat (file from pseudobulk_ch)
            pseudobulk_singlecell.out.gene_locations,  // gene_locations
            combined_ch.map { it[1] }  // n_pcs (value from n_pcs_ch)
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
        // Pass the collected results to select_pcs
        select_pcs(grouped_results)
        
        // this is for the markdown file
        optimize_pcs.out.egenes_results
            .collect()
            .set { collected_results }

            
    }else{
            // Define a dummy file for optimization_results
        dummy_file = file("dummy_optimization_results.txt")

        // Initialize collected_results with the dummy file
        collected_results = Channel.from(dummy_file)
    }

    run_matrixeQTL(
        params.eqtl_source_functions,
        qc_genotype.out.qc_genotype_mat,
        qc_genotype.out.qc_snp_chromlocations,
        qc_expression.out.pseudobulk_normalised.flatten(),
        pseudobulk_singlecell.out.gene_locations,
        params.cov_file
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
    // ...existing code...
}

workflow.onComplete {
    println """
    ========================================
    Pipeline Completed!
    ========================================

    """
}

