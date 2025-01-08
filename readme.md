Hi, first test!

To run the pipeline:

- Clone this repo
- Run the following command in the terminal

```sh
git pull; nextflow run -c nextflow.config main.nf --gds_file path_to_genotype_gds_file --single_cell_file path_to_single_cell_file --outdir path_to_output_directory -w $ephemeral
```

### Input files / output dirs


- **outdir**: Output directory for the pipeline results. Default: `/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/`
- **gds_file**: Path to the GDS file containing genotype data. Default: `/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds`
- **single_cell_file**: Path to the single-cell data file. Default: `/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/roche_ms_decontx.rds`

### Run Parameters

You don't need to add the parameters necessarily, there are default values.

- **counts_assay**: Assay used for counts. Default: `RNA`
- **counts_slot**: Slot used for counts. Default: `counts`
- **celltype_column**: Column name for cell types within the Seurat object. Default: `CellType`
- **individual_column**: Column name for individual IDs. Default: `Individual_ID`
- **min_cells**: Minimum number of cells for pseudobulking. Default: `10`
- **min_expression**: Minimum expression percentage for genes. Genes expressed in less than this percentage of individuals excluded. Default: `0.05`
- **cis_distance**: Cis distance for eQTL analysis. Default: `1e6`
- **fdr_threshold**: FDR threshold for eQTL results. Default: `0.05`
- **optimize_pcs**: Boolean to optimize principal components. This iteratively runs MatrixEQTL to find optimal PCs. Default: `true`

### Other useful flags

- **-w**: Working directory. This is where nextflow stages input files. Needs lots of space! Default: `/rds/general/user/$USER/ephemeral/`
- **-N <email.address>**: Email address to send notifications to on run completion. For example; "a.haglund19@imperial.ac.uk". 

