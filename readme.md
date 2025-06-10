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

### Analysis Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--counts_assay` | Assay for counts | `RNA` |
| `--counts_slot` | Slot for counts | `counts` |
| `--celltype_column` | Column for cell types | `CellType` |
| `--individual_column` | Column for individual IDs | `Individual_ID` |
| `--min_cells` | Min cells for pseudobulking | `10` |
| `--min_expression` | Min expression percentage | `0.05` |
| `--cis_distance` | Cis distance for eQTL analysis | `1e6` |
| `--fdr_threshold` | FDR threshold | `0.05` |
| `--optimize_pcs` | Optimize principal components | `true` |

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

The pipeline uses containerization to ensure reproducibility across environments:

| Container | Purpose | Repository |
|-----------|---------|------------|
| **eqtl-genotype** | Genotype processing & eQTL analysis | `ghcr.io/johnsonlab-ic/eqtl-genotype:latest` |
| **eqtl-expression** | Single-cell data & pseudobulking | `ghcr.io/johnsonlab-ic/eqtl-expression:latest` |
| **eqtl-report** | Report generation & visualization | `ghcr.io/johnsonlab-ic/eqtl-report:latest` |

### Building Custom Containers

```bash
# Build containers locally
cd Images
docker build -t local/eqtl-expression:latest -f Dockerfile.expression .
docker build -t local/eqtl-genotype:latest -f Dockerfile.genotype .
docker build -t local/eqtl-report:latest -f Dockerfile.reports .
```

To use custom images, modify the corresponding entries in `nextflow.config`.

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

The pipeline generates these key outputs:

- **eQTL results**: Lists of significant eQTLs for each cell type
- **Visualization reports**: Interactive HTML reports with plots
- **Optimization data**: PC optimization results if enabled
- **QC metrics**: Quality control information for genotype and expression data

---

## 📞 Support

For questions or issues, please:
- Open an issue on the [GitHub repository](https://github.com/johnsonlab-ic/sceQTL-flow/issues)
- Contact the Johnson Lab at Imperial College London
