


## Welcome

This pipeline is designed for eQTL (expression Quantitative Trait Loci) analysis using single-cell RNA sequencing data. It integrates genotype data with single-cell expression data to identify genetic variants that influence gene expression at the single-cell level.

### Required Input Files

- **Genotype GDS File**: A GDS (Genomic Data Structure) file containing genotype data.
- **Single-Cell Data File**: An RDS file containing single-cell RNA sequencing data in a Seurat object format.

### Pipeline Functionality

The pipeline performs the following steps:

1. **Data Preprocessing**: Reads and preprocesses the input genotype and single-cell data files.
2. **Pseudobulking**: Aggregates single-cell data into pseudobulk samples based on cell types and individuals.
3. **eQTL Analysis**: Conducts eQTL analysis to identify genetic variants associated with gene expression levels.
4. **Result Optimization**: Optionally optimizes the number of principal components for the eQTL analysis.
5. **Output Generation**: Produces results including significant eQTLs, summary statistics, and diagnostic plots.

### Usage

Make sure you have nextflow installed! `curl -s https://get.nextflow.io | bash` 

To run the pipeline:

1. Clone this repository. Run;

    ```sh
    git clone https://github.com/johnsonlab-ic/sc-eQTL-pipepeline/ 
    cd sc-eQTL-pipepeline
    ```

2. Run the following command in the terminal, ensuring you provide the correct paths for the inputs:

    ```sh
    git pull; nextflow run -c nextflow.config main.nf -profile imperial \
    --gds_file path_to_genotype_gds_file \
    --single_cell_file path_to_single_cell_file \
    --outdir path_to_output_directory \
    --N email@adress
    ```

### Input Files / Output Directories

- **outdir**: Output directory for the pipeline results. Default: `/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/`
- **gds_file**: Path to the GDS file containing genotype data. Default: `/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds`
- **single_cell_file**: Path to the single-cell data file. Default: `/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/roche_ms_decontx.rds`

### Run Parameters / Flags

You don't need to add the parameters necessarily, there are default values.

- **--counts_assay**: Assay used for counts. Default: `RNA`
- **--counts_slot**: Slot used for counts. Default: `counts`
- **--celltype_column**: Column name for cell types within the Seurat object. Default: `CellType`
- **--individual_column**: Column name for individual IDs. Default: `Individual_ID`
- **--min_cells**: Minimum number of cells for pseudobulking. Default: `10`
- **--min_expression**: Minimum expression percentage for genes. Genes expressed in less than this percentage of individuals excluded. Default: `0.05`
- **--cis_distance**: Cis distance for eQTL analysis. Default: `1e6`
- **--fdr_threshold**: FDR threshold for eQTL results. Default: `0.05`
- **--optimize_pcs**: Boolean to optimize principal components. This iteratively runs MatrixEQTL to find optimal PCs. Default: `true`

### Other Useful Flags

- **-w**: Working directory. This is where nextflow stages input files. Needs lots of space! Default: `/rds/general/user/$USER/ephemeral/`
- **-N <email.address>**: Email address to send notifications to on run completion. For example: `a.haglund19@imperial.ac.uk`.



### Notes / warnings!!

This is designed to be run on the Imperial College HPC system (for the time being) due to being very memory-intensive but will vary by dataset.