nextflow.enable.dsl=2

// default inputs
params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds"
params.single_cell_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/roche_ms_decontx.rds"
params.cov_file="none"
params.help = false
params.workflow = 'matrixeqtl'

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
params.optimize_pcs = true // Whether to optimize PCs or use a fixed number
params.fixed_pcs = 10 // Number of PCs to use when not optimizing
params.report = false

// default parameters 
// Add a new parameter for specifying which covariates to include
params.covariates_to_include = "all" // Default to include all covariates

// Help message similar to scQC-flow style
def helpMessage() {
        log.info """
        ========================================
        sceQTL-flow
        ========================================

        Usage:
            nextflow run main.nf \\
                --workflow matrixeqtl|tensorqtl \\
                --gds_file <path.gds> \\
                --single_cell_file <seurat.rds> \\
                --outdir <output_dir> \\
                --celltype_column <column> \\
                --individual_column <column> \\
                [--cov_file covariates.csv] \\
                [--covariates_to_include <comma list|all>] \\
                [--optimize_pcs true|false] \\
                [--fixed_pcs 10]

        Flags:
            --help              Show this message
            --report            Render the HTML report

        Notes:
            Residuals are calculated automatically when --cov_file is provided.
            Set --optimize_pcs false to force a fixed number of PCs (default 10).
            tensorqTL workflow is a placeholder.
        """
}

if (params.help) {
    helpMessage()
    System.exit(0)
}

include { matrixeqtl } from './workflows/matrixeqtl.nf'
include { tensorqtl } from './workflows/tensorqtl.nf'

workflow {
    def wf = params.workflow?.toLowerCase() ?: 'matrixeqtl'
    switch(wf) {
        case 'tensorqtl':
            log.info "Running tensorQTL workflow (placeholder)."
            tensorqtl()
            break
        case 'matrixeqtl':
            log.info "Running matrixeQTL workflow."
            matrixeqtl()
            break
        default:
            error "Unknown workflow '${wf}'. Choose 'matrixeqtl' or 'tensorqtl'."
    }
}

