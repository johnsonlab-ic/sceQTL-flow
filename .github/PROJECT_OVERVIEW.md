# sceQTL-flow — Project Overview

## Project Overview

**Pipeline:** Nextflow workflow for single-cell expression quantitative trait loci (sc-eQTL) analysis  
**Purpose:** Map genetic variants to gene expression changes at single-cell resolution  
**Platform:** Nextflow DSL2 with containerized modules (Docker/Singularity)  
**Analysis Framework:** MatrixEQTL, TensorQTL, R, Python

---

## Pipeline Overview

| Component | Description |
|-----------|-------------|
| **Input Data** | VCF/GDS genotype files + Seurat single-cell objects |
| **Preprocessing** | Pseudobulk aggregation, covariate preparation, genotype filtering |
| **eQTL Testing** | MatrixEQTL and TensorQTL for variant-expression associations |
| **Outputs** | eQTL results, QC reports, visualizations |

---

## Pipeline Modules

### Core Workflow Modules

| Module | File | Description |
|--------|------|-------------|
| **Genotype Processing** | `modules/genotype/genotype.nf` | VCF/GDS conversion, filtering, quality control |
| **Pseudobulk** | `modules/expression/pseudobulk.nf` | Aggregate single-cell counts to pseudobulk by cell type |
| **Expression QC** | `modules/expression/qc_expression.nf` | Normalize and QC pseudobulk matrices |
| **Residuals** | `modules/residuals/get_residuals.nf` | Calculate expression residuals after covariate correction |
| **MatrixEQTL** | `modules/eqtl/matrixeqtl.nf` | Run MatrixEQTL for eQTL mapping |
| **Combine Results** | `modules/eqtl/combine_eqtls.nf` | Merge eQTL results across cell types |
| **Final Report** | `modules/reports/final_report.nf` | Render HTML report (R Markdown) |

### Optimization Modules

| Module | File | Description |
|--------|------|-------------|
| **Count Individuals** | `modules/optimize/count_individuals.nf` | Count individuals per cell type to set coarse PC grid |
| **Optimize PCs (coarse/fine)** | `modules/optimize/optimize_pcs.nf` | Evaluate eGene counts across coarse/fine PC grids |
| **Select PCs (coarse)** | `modules/optimize/select_pcs_coarse.nf` | Early-stop coarse grid and generate fine grid |
| **Select PCs (fine)** | `modules/optimize/select_pcs.nf` | Apply elbow rule to select final PCs |
| **Organize PC Summaries** | `modules/reports/organize_pc_optimization.nf` | Collect coarse/fine summaries for reporting |

## Key Input Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--gds_file` | Genotype file (GDS format) | `genotypes.gds` |
| `--single_cell_file` | Seurat object with single-cell data | `seurat_object.rds` |
| `--single_cell_file_list` | Comma-separated list of Seurat objects | `file1.rds,file2.rds` |
| `--celltype_column` | Metadata column for cell types | `"cell_type"` |
| `--individual_column` | Metadata column for individual IDs | `"Individual_ID"` |
| `--cov_file` | Covariate file (optional) | `covariates.csv` |
| `--covariates_to_include` | Comma-separated covariates or `all` | `age,sex` |
| `--subset_column` | Metadata column to subset on | `"Diagnosis"` |
| `--subset_values` | Values to keep (comma-separated) | `"Control,AD"` |
| `--counts_assay` | Seurat assay to use | `"RNA"` |
| `--counts_slot` | Seurat slot to use | `"counts"` |
| `--cis_distance` | cis-eQTL window size | `1000000` |
| `--fdr_threshold` | FDR cutoff | `0.05` |
| `--filter_chr` | Chromosome filter | `"chr6"` |
| `--optimize_pcs` | Enable PC optimization | `true` |
| `--fixed_pcs` | Fixed PCs when not optimizing | `10` |
| `--pc_coarse_step` | Coarse grid step | `10` |
| `--pc_fine_step` | Fine grid step | `2` |
| `--pc_fine_window` | Fine window around best coarse | `10` |
| `--pc_elbow_tol` | Elbow tolerance | `0.02` |
| `--pc_early_stop_tol` | Early-stop tolerance | `0.01` |
| `--pc_early_stop_patience` | Early-stop patience | `2` |
| `--report` | Render final HTML report | `true` |

## Configuration Profiles

Located in `nextflow.config`:

| Profile | Description |
|---------|-------------|
| `imperial` | PBS Pro execution (Imperial HPC) |
| `standard` | Local execution with Docker |
| `offline` | Local execution without containers |

## Directory Structure

```
├── main.nf                      # Main pipeline workflow
├── nextflow.config              # Global configuration
├── readme.md                    # Pipeline documentation
├── run_hpc.sh                   # Example HPC execution scripts
├── run_dsi.sh                   # Example local execution script
├── modules/                     # Pipeline modules (DSL2)
│   ├── genotype/
│   ├── expression/
│   ├── residuals/
│   ├── eqtl/
│   ├── reports/
│   ├── optimize/
│   └── qc/
├── workflows/                   # High-level workflows
│   ├── matrixeqtl.nf
│   └── tensorqtl.nf
├── R/                           # R scripts and functions
│   ├── expression_functions/
│   ├── genotype_functions/
│   ├── MatrixEQTL_functions/
│   └── rmarkdown_reports/
│       └── unified_final_report.Rmd  # R Markdown report template
├── Images/                      # Docker containers
│   ├── Dockerfile.expression
│   ├── Dockerfile.genotype
│   ├── Dockerfile.python
│   ├── Dockerfile.reports
│   ├── Dockerfile.tensorqtl
│   └── setup_registry.sh
└── .github/                     # GitHub files
  ├── PROJECT_OVERVIEW.md      # This file
  └── copilot-instructions.md  # Copilot instructions
```

---

## Key Output Files

| File | Description |
|------|-------------|
| `pseudobulk_{celltype}.rds` | Pseudobulk expression matrices per cell type |
| `genotype_filtered.gds` | QC-filtered genotype data |
| `covariates_{celltype}.csv` | Covariate matrices per cell type |
| `eqtl_results_{celltype}.csv` | eQTL results per cell type |
| `pipeline_report.html` | Final QC and results report |

---

## Compute Resources

### PBS Pro (Imperial HPC)

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| Genotype processing | 4 | 16 GB | 4 h |
| Pseudobulk | 8 | 32 GB | 6 h |
| MatrixEQTL | 16 | 64 GB | 12 h |

---

## Methods

- **Genotype Processing:** SNPRelate, SeqArray (GDS format)
- **Pseudobulk Aggregation:** Seurat v5 aggregation by cell type and individual
- **eQTL Testing:** MatrixEQTL (cis and trans), TensorQTL (optional)
- **Covariate Correction:** PCA-based covariates, batch effect correction
- **Reporting:** R Markdown reports

---

## Example Usage

### Local execution (with Docker)
```bash
nextflow run main.nf \
  --gds_file genotypes.gds \
  --single_cell_file seurat_object.rds \
  --celltype_column "cell_type" \
  --individual_column "Individual_ID" \
  --outdir ./results
```

### HPC execution (PBS Pro)
```bash
nextflow run -c nextflow.config main.nf \
  -profile imperial \
  --gds_file genotypes.gds \
  --single_cell_file seurat_object.rds \
  --celltype_column "cell_type" \
  --individual_column "Individual_ID" \
  --outdir /path/to/results
```

---

## Current Development Status

This is an active Nextflow pipeline for sc-eQTL analysis. Key features include:
- ✅ Modular Nextflow DSL2 architecture
- ✅ Docker/Singularity containerization
- ✅ HPC scheduler support (PBS Pro)
- ✅ MatrixEQTL integration
- ✅ Automated reporting with R Markdown
- 🔄 TensorQTL support (in development)
- ✅ Coarse-to-fine PC optimization workflows
