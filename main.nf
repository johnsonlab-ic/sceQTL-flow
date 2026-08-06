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
params.pseudobulk_anndata_script="${baseDir}/R/expression_functions/pseudobulk_anndata.py"
params.pseudobulk_seurat_script="${baseDir}/R/expression_functions/pseudobulk_seurat.R"
params.combine_pseudobulk_script="${baseDir}/R/expression_functions/combine_pseudobulk.R"
params.check_overlap_script="${baseDir}/R/expression_functions/check_overlap.py"
params.eqtl_source_functions="${baseDir}/R/MatrixEQTL_functions/matrixeqtl_source.r"
params.quarto_report="${baseDir}/R/rmarkdown_reports/unified_report.qmd"
params.pc_optimization_report="${baseDir}/R/rmarkdown_reports/pc_optimization.qmd"

// cell-level metadata + grouping (optional external metadata file)
params.cell_metadata_file="none"   // CSV/.gz keyed by metadata_id_col; if "none", labels read from object
params.metadata_id_col="cell_id"
params.overlap_warn_frac=0.5       // warn if genotype<->single-cell overlap below this fraction

// optional genotype sample relabelling (default: expect IDs to already match)
params.sample_map="none"           // CSV mapping genotype IDs -> individual IDs
params.sample_map_from="Sample_ID"
params.sample_map_to="caseid"

params.min_cells=5
params.min_expression=0.1
params.celltype_column="celltype" // comma-separated to pseudobulk on multiple columns/resolutions independently
params.individual_column="individual"
params.counts_assay="RNA"
params.counts_slot="counts"
params.cis_distance=1000000
params.fdr_threshold=0.05
params.save_full_eqtl=false  // also persist full per-celltype cis stats (for coloc/mashr); never concatenated
params.filter_chr = "all" // Optional parameter for filtering by chromosome. use "chr6"
params.optimize_pcs = true // Whether to optimize PCs or use a fixed number
params.fixed_pcs = 10 // Number of PCs to use when not optimizing
params.report = false

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
                --celltype_column <column|column1,column2,...> \\
                --individual_column <column> \\
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
            tensorqTL workflow is a placeholder.
        """
}

include { matrixeqtl } from './workflows/matrixeqtl.nf'
include { tensorqtl } from './workflows/tensorqtl.nf'

workflow {
    if (params.help) {
        helpMessage()
        System.exit(0)
    }

    def wf = params.workflow?.toLowerCase() ?: 'matrixeqtl'
    if (wf == 'tensorqtl') {
        log.info "Running tensorQTL workflow (placeholder)."
        tensorqtl()
    } else if (wf == 'matrixeqtl') {
        log.info "Running matrixeQTL workflow."
        matrixeqtl()
    } else {
        error "Unknown workflow '${wf}'. Choose 'matrixeqtl' or 'tensorqtl'."
    }
}


