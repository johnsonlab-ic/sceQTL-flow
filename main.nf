nextflow.enable.dsl=2

// default inputs
params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds"
params.single_cell_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/roche_ms_decontx.rds"
params.single_cell_file_list="none"
params.cov_file="none"
params.help = false
params.workflow = 'matrixeqtl'

// source functions for easy troubleshooting
params.genotype_source_functions="${baseDir}/R/genotype_functions/genotype_functions.r"
params.pseudobulk_source_functions="${baseDir}/R/expression_functions/pseudobulk_functions.r"
params.atac_source_functions="${baseDir}/R/atac_functions/pseudobulk_functions.r"
params.eqtl_source_functions="${baseDir}/R/MatrixEQTL_functions/matrixeqtl_source.r"
params.quarto_report="${baseDir}/R/rmarkdown_reports/unified_final_report.Rmd"
params.quarto_report_atac="${baseDir}/R/rmarkdown_reports/unified_final_report_atac.Rmd"

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
params.data_type = "RNA" // RNA (default) or ATAC

// PC optimization strategy parameters
params.pc_max_fraction = 0.5
params.pc_max_cap = 100
params.pc_min = 2
params.pc_coarse_step = 10
params.pc_fine_step = 2
params.pc_fine_window = 10
params.pc_elbow_tol = 0.02 // 2% within max
params.pc_early_stop_tol = 0.01 // 1% improvement threshold
params.pc_early_stop_patience = 2

// default parameters 
// Add a new parameter for specifying which covariates to include
params.covariates_to_include = "all" // Default to include all covariates
params.subset_column = "none" // Column to subset by (e.g., "Diagnosis")
params.subset_values = "none" // Values to keep (e.g., "Control" or "Control,AD")

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
                --single_cell_file <seurat.rds> OR --single_cell_file_list <file1.rds,file2.rds,...> \\
                --outdir <output_dir> \\
                --celltype_column <column> \\
                --individual_column <column> \\
                --data_type RNA|ATAC \\
                [--cov_file covariates.csv] \\
                [--covariates_to_include <comma list|all>] \\
                [--optimize_pcs true|false] \\
                [--fixed_pcs 10] \\
                [--pc_coarse_step 10] \\
                [--pc_fine_step 2] \\
                [--pc_fine_window 10] \\
                [--pc_elbow_tol 0.02] \\
                [--pc_early_stop_tol 0.01] \\
                [--pc_early_stop_patience 2]

        Flags:
            --help              Show this message
            --report            Render the HTML report

        Notes:
            Residuals are calculated automatically when --cov_file is provided.
            Set --optimize_pcs false to force a fixed number of PCs (default 10).
            tensorqTL workflow currently runs with fixed PCs only (no PC optimization, no report).
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
            log.info "Running tensorQTL workflow (fixed PCs only; no optimization/report)."
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

