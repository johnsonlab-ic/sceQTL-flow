nextflow.enable.dsl=2

params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds"

params.genotype_source_functions="${baseDir}/R/genotype_functions/genotype_functions.r"
params.pseudobulk_source_functions="${baseDir}/R/expression_functions/pseudobulk_functions.r"
params.eqtl_source_functions="${baseDir}/R/MatrixEQTL_functions/matrixeqtl_source.r"
params.quarto_report="${baseDir}/R/rmarkdown_reports/final_report.Rmd"

params.single_cell_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/roche_ms_decontx.rds"

// ...existing code...

include { create_genotype } from './NEXTFLOW/genotype.nf'
include { pseudobulk_singlecell } from './NEXTFLOW/pseudobulk.nf'
include { qc_expression } from './NEXTFLOW/qc_expression.nf'
include { run_matrixeQTL } from './NEXTFLOW/matrixeqtl.nf'
include { combine_eqtls } from './NEXTFLOW/combine_eqtls.nf'
include { final_report } from './NEXTFLOW/final_report.nf'

// ...existing code...

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

    If "optimize PCs" is set to TRUE, the pipeline will run longer.

    Generating reports via markdown too now!

    ==============================================

    This is the stable release.

    """
    create_genotype(gds_file= params.gds_file)
    pseudobulk_singlecell(single_cell_file= params.single_cell_file)
    pseudobulk_ch=pseudobulk_singlecell.out.pseudobulk_counts.flatten()
    qc_expression(pseudobulk_file= pseudobulk_ch)
    run_matrixeQTL(
        genotype_mat= create_genotype.out.genotype_mat,
        snp_locations= create_genotype.out.snp_chromlocations,
        expression_mat= qc_expression.out.pseudobulk_normalised.flatten(),
        gene_locations= pseudobulk_singlecell.out.gene_locations
    )
    combine_eqtls(eqtls= run_matrixeQTL.out.eqtl_results.collect())
    if(params.report){
        final_report(
            eqtl_results_filtered = combine_eqtls.out.mateqtlouts_FDR_filtered,
            eqtl_results = combine_eqtls.out.mateqtlouts,
            report_file = params.quarto_report
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

