# Copilot Instructions for sceQTL-flow

## Project Context

Before answering questions about this project, read `/workspace/.github/PROJECT_OVERVIEW.md` for full context including:
- Pipeline architecture and modules
- Input/output data formats
- Configuration profiles
- Directory structure and key files

## Key Information

- **Project:** sceQTL-flow — Single-cell eQTL analysis pipeline
- **Framework:** Nextflow DSL2 with modular architecture
- **Languages:** Nextflow, R, Python, Shell
- **Containers:** Docker (local), Singularity (HPC)
- **Compute:** PBS Pro (Imperial HPC) and local execution
- **Key Dependencies:** Seurat v5, MatrixEQTL, SeqArray, TensorQTL

## Code Style Preferences

### Nextflow
- Use DSL2 syntax with modular processes
- Separate processes into logical workflow files in `NEXTFLOW/`
- Use process labels for resource allocation (e.g., `label 'process_high_memory'`)
- Support multiple execution profiles (imperial, standard, offline)
- Use `publishDir` with explicit modes (`copy`, `symlink`)
- Channel operators should be explicit and documented

### R Scripts
- Use tidyverse/dplyr for data manipulation
- Use Seurat v5 conventions for single-cell analysis
- Prefer `file.path()` for path construction
- Use explicit variable names (e.g., `celltype_column` not `ct`)
- Add comments for complex logic
- Handle errors gracefully with informative messages
- Load required packages at the start of scripts

### Python Scripts
- Follow PEP 8 style guidelines
- Use type hints where appropriate
- Prefer pandas for data manipulation
- Use pathlib for file paths
- Add docstrings for functions

### Shell Scripts
- Use `#!/bin/bash` shebang
- Use UPPERCASE for variables (e.g., `OUTPUT_DIR`)
- Always quote variables: `"$var"`
- Use `set -e` for error handling in critical scripts
- Add comments for pipeline steps

## File Locations

- **Main workflow:** [main.nf](../main.nf)
- **Pipeline workflows:** [workflows/](../workflows/) (matrixeqtl.nf, tensorqtl.nf)
- **Pipeline modules:** [modules/](../modules/) (genotype/, expression/, eqtl/, residuals/, optimize/, reports/, qc/)
- **R functions & reports:** [R/](../R/) (expression_functions/, genotype_functions/, MatrixEQTL_functions/, rmarkdown_reports/, quarto_reports/)
- **Docker images:** [Images/](../Images/)
- **Config files:** [nextflow.config](../nextflow.config)
- **Example run scripts:** [personal/](../personal/)

## Common Tasks

When asked to:

1. **Add a new pipeline module** — Create `.nf` file in `modules/{category}/`, add process with appropriate labels, import in `workflows/matrixeqtl.nf`
2. **Modify genotype processing** — Update [modules/genotype/genotype.nf](../modules/genotype/genotype.nf) and related R functions in [R/genotype_functions/](../R/genotype_functions/)
3. **Change pseudobulk logic** — Update [R/expression_functions/pseudobulk_functions.r](../R/expression_functions/pseudobulk_functions.r)
4. **Add eQTL testing methods** — Update [modules/eqtl/matrixeqtl.nf](../modules/eqtl/matrixeqtl.nf) or create new module in [modules/eqtl/](../modules/eqtl/)
5. **Modify PC optimization** — Update [modules/optimize/](../modules/optimize/) files
6. **Modify resource allocation** — Update process labels in `nextflow.config` or create profile-specific configs
7. **Update reporting** — Modify [R/rmarkdown_reports/unified_final_report.Rmd](../R/rmarkdown_reports/unified_final_report.Rmd)

### PC Optimization Parameters (when `--optimize_pcs true`)
- `--pc_coarse_step`: Step size for coarse grid (default: 10)
- `--pc_fine_step`: Step size for fine grid (default: 2)
- `--pc_fine_window`: Window around best coarse PC (default: 10)
- `--pc_elbow_tol`: Elbow tolerance as fraction of max (default: 0.02 = 2%)
- `--pc_early_stop_tol`: Early-stop tolerance for gains (default: 0.01 = 1%)
- `--pc_early_stop_patience`: Consecutive low-gain steps to halt coarse grid (default: 2)
- `--gds_file`: Genotype file in GDS format (or VCF with auto-conversion)
- `--single_cell_file`: Seurat object with single-cell RNA-seq data
- `--celltype_column`: Metadata column name for cell type annotation
- `--individual_column`: Metadata column name for individual/donor IDs

### Optional Inputs
- `--cov_file`: External covariate file (CSV format)
- `--covariates_to_include`: Comma-separated covariate names or "all" (default: "all")
- `--subset_column`: Column name to subset samples (e.g., "Diagnosis")
- `--subset_values`: Values to keep for subsetting (comma-separated)
- `--counts_assay`: Seurat assay name (default: "RNA")
- `--counts_slot`: Seurat slot name (default: "counts")
- `--optimize_pcs`: Enable PC optimization workflow (default: true)
- `--fixed_pcs`: Fixed number of PCs when `--optimize_pcs false` (default: 10)
- `--report`: Generate unified HTML report (default: false)

## Pipeline Workflow

1. **Genotype Processing** → Filter and QC genotype data
2. **Pseudobulk Aggregation** → Sum counts by cell type and individual
3. **Sample Subsetting** → Optional filtering by condition/diagnosis (preflight check)
4. **Covariate Preparation** → Extract covariates from metadata and/or file
5. **Residual Calculation** → Regress out covariates if specified
6. **PC Optimization (conditional)** → Two-stage coarse-to-fine grid search for optimal PCs per cell type
   - Coarse stage: Early stopping at 1% improvement plateau
   - Fine stage: Elbow tolerance rule to select final PCs
7. **eQTL Testing** → Run MatrixEQTL for variant-expression associations per cell type
8. **Combine Results** → Merge eQTL results across all cell types
9. **Report Generation** → Create unified HTML report (eQTL + PC optimization results)

## Configuration Profiles

| Profile | Executor | Container Engine | Use Case |
|---------|----------|------------------|----------|
| `standard` | local | docker | Local execution with Docker |
| `imperial` | pbspro | singularity | Imperial HPC with PBS Pro |
| `offline` | local | none | Local execution without containers |

## Important Notes

- GDS format is preferred for genotype data (VCF will be auto-converted)
- Each cell type is processed independently for eQTL testing
- Minimum individuals per cell type should be ≥20 for robust eQTL detection
- Pseudobulk aggregation sums counts per individual per cell type
- **PC Optimization** (when enabled):
  - Two-stage grid search: coarse (large steps) then fine (small steps)
  - Coarse stage uses early stopping to detect plateau
  - Fine stage applies elbow rule to balance accuracy and parsimony
  - Summary CSVs saved to `{outdir}/optimization/` for visualization
  - Each cell type can have different optimal PC counts

## Before Making Changes

- Identify which module(s) will be affected
- Consider impact on upstream/downstream processes
- Test changes on small dataset before full runs
- Update documentation if adding new parameters or features
- Check container dependencies if modifying R/Python code

## Response Guidelines

**Keep code focused:**
- Provide minimal working implementations
- Avoid unnecessary complexity or extra features
- Only add configuration options if explicitly requested
- Present alternative approaches concisely (2-3 options max)

**Explain clearly:**
- Describe changes and their purpose in 2-3 sentences
- Explain how modifications fit into the pipeline workflow
- Note any dependencies or prerequisites
- Mention expected outputs or behavior changes

## Data Format Expectations

### Genotype Data
- Format: GDS (preferred) or VCF
- Variants: SNPs, indels
- Required fields: CHROM, POS, ID, REF, ALT, genotypes

### Single-Cell Data
- Format: Seurat object (RDS)
- Required assays: RNA counts
- Required metadata: cell type, individual ID
- Optional metadata: batch, covariates

### Covariates
- Format: CSV with individual IDs in first column
- Column names must match covariate names
- Numeric or categorical variables

## Testing Recommendations

When modifying the pipeline:
1. Test on small dataset (1-2 cell types, <1000 cells)
2. Verify outputs are generated correctly
3. Check log files for errors or warnings
4. Run full dataset only after validation
5. Compare results with previous runs when applicable

## Container Management

- Docker images defined in [../Images/](../Images/)
- Update container versions in `nextflow.config` when modifying images
- Rebuild containers after changing Dockerfiles
- Push to registry before running on HPC

## HPC-Specific Notes

### PBS Pro (Imperial HPC)
- Queue selection in process directives
- Walltime limits: adjust `time` in process labels
- Storage: Use `$EPHEMERAL` for temporary files, persistent storage for outputs
- Module loading: Specified in profile config
