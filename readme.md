# 🧬 sceQTL-flow

<img src="https://img.shields.io/badge/Nextflow-v22.10.0+-green.svg" alt="Nextflow Version">
<img src="https://img.shields.io/badge/R-v4.3.0+-blue.svg" alt="R Version">
<img src="https://img.shields.io/badge/Containers-Docker%2FSingularity-orange.svg" alt="Container Support">


**A Nextflow pipeline for single-cell expression Quantitative Trait Loci analysis**

This pipeline integrates genotype data with single-cell RNA sequencing data to identify genetic variants that influence gene expression at the single-cell level.


## 🚀 Quick Start

### Installation

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
```

### Running the Pipeline

```bash
nextflow run johnsonlab-ic/sceQTL-flow \
  --gds_file path_to_genotype_gds_file \
  --single_cell_file path_to_single_cell_file \
  --outdir path_to_output_directory 
```

## 📁 Input Files & Parameters

### Input and Output Paths

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--outdir` | Output directory | `/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/` |
| `--gds_file` | Genotype data file | `/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds` |
| `--single_cell_file` | Seurat object file | `/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/roche_ms_decontx.rds` |

### Input Formats

`--single_cell_file` (one file) or `--single_cell_file_list` (comma-separated) accept either Seurat `.rds` or AnnData `.h5ad` objects — mixed lists of both are fine. Each file is pseudobulked natively (no merge-into-one-object step), then combined.

### Analysis Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--counts_assay` | Assay for counts (Seurat input only) | `RNA` |
| `--counts_slot` | Slot/layer for raw counts | `counts` |
| `--celltype_column` | Column(s) for cell types. Comma-separated to pseudobulk multiple annotation columns/resolutions independently in one run (e.g. `col1,col2`) | `celltype` |
| `--individual_column` | Column for individual IDs | `individual` |
| `--cell_metadata_file` | Optional external CSV/.gz of cell-level labels (celltype + individual), keyed by `--metadata_id_col`. If omitted, labels are read from the object's own metadata | `none` |
| `--metadata_id_col` | Cell-id column in `--cell_metadata_file` (must match the object's cell/barcode names) | `cell_id` |
| `--sample_map` | Optional CSV to relabel genotype sample IDs to individual IDs, when they don't already match | `none` |
| `--sample_map_from` / `--sample_map_to` | Column names in `--sample_map` for the relabeling | `Sample_ID` / `caseid` |
| `--overlap_warn_frac` | Warn if genotype↔single-cell ID overlap falls below this fraction (only runs when `--cell_metadata_file` is set) | `0.5` |
| `--min_cells` | Min total cells for an individual to be kept in pseudobulking | `10` |
| `--min_expression` | Min expression percentage | `0.05` |
| `--cis_distance` | Cis distance for eQTL analysis | `1e6` |
| `--fdr_threshold` | FDR threshold | `0.05` |
| `--save_full_eqtl` | Also persist the full (not just FDR-significant) per-celltype cis association table, as separate un-combined `<celltype>_cis_MatrixEQTLout.rds` files | `false` |
| `--optimize_pcs` | Optimize principal components | `true` |

Celltypes with fewer than 15 individuals remaining after pseudobulking are dropped automatically (logged, not fatal).

### Runtime Options

| Option | Description | Default |
|--------|-------------|---------|
| `-w` | Working directory | `/rds/general/user/$USER/ephemeral/` |
| `-N` | Email for notifications | _none_ |

---

---

## 📋 Overview

The pipeline performs these key steps:

1. **Data Preprocessing** 
   * Reads and processes input files
   * Performs quality control on genotype data

2. **Pseudobulking**
   * Aggregates single-cell data by cell type and individual
   * Normalizes expression data

3. **eQTL Analysis**
   * Identifies genetic variants associated with gene expression
   * Conducts statistical testing with MatrixEQTL

4. **Result Optimization**
   * Optimizes principal components for each cell type
   * Enhances detection power and accuracy

5. **Output Generation**
   * Produces significant eQTLs lists
   * Creates summary statistics and diagnostic plots


## ⚠️ Notes & Warnings

> **System Requirements**: This pipeline is computationally intensive and best run on HPC systems.

This pipeline is optimized for the Imperial College HPC system due to its memory-intensive operations, but it can be adapted to other systems with sufficient resources.

---

## 🐳 Docker Images

The pipeline uses two containers:

| Container | Purpose | Repository |
|-----------|---------|------------|
| **genetics** | Everything by default: genotype processing, PC optimization, matrixQTL, combine/report | `ghcr.io/johnsonlab-ic/genetics:latest` |
| **landmark-sc_image** | Single-cell-heavy steps only: `pseudobulk_anndata`, `pseudobulk_seurat`, `check_overlap` | `ghcr.io/johnsonlab-ic/landmark-sc_image:latest` |

These images are built and published outside this repository. To use custom images, change the `container` entries in `nextflow.config`.

---

## 📚 Repository Information

<img src="https://img.shields.io/badge/GitHub-sceQTL--flow-lightgrey?logo=github" alt="GitHub Repo">

This pipeline is maintained in a public repository:
- [johnsonlab-ic/sceQTL-flow](https://github.com/johnsonlab-ic/sceQTL-flow)

### Contributing

We welcome contributions! Please follow these steps:

1. 🔍 **Open an issue** describing the feature or bug
2. 🍴 **Fork** the repository
3. 🌿 **Create a branch** for your changes
4. 🔄 **Submit a pull request**

For major changes, please discuss them first via issues.

---

## 📊 Output Details

Under `<outdir>/eQTL_outputs/`:

- `mateqtlouts_FDR_filtered.rds` — significant (FDR-passing) cis associations per celltype
- `eqtl_summary.rds` / `.csv` — per-celltype counts (n_individuals, n_tests, n_sig_pairs, n_egenes, ...)
- `<celltype>_cis_MatrixEQTLout.rds` — full unfiltered per-celltype association table, only if `--save_full_eqtl true` (one file per celltype, never combined)
- `eqtl_report.html` — the unified report (celltype QC, PC optimization, cells-per-individual chart, results), if `--report true`

Also written: `<outdir>/QC/overlap_report.txt` (genotype↔single-cell ID overlap, if `--cell_metadata_file` set) and `<outdir>/run_params.txt` (a provenance manifest of the exact parameters and revision used).

---

## 📞 Support

For questions or issues, please:
- Open an issue on the [GitHub repository](https://github.com/johnsonlab-ic/sceQTL-flow/issues)
- Contact the Johnson Lab at Imperial College London
