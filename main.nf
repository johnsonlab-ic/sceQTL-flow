nextflow.enable.dsl=2

params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds"
params.inputfile="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/eqtl_pipeline_dev/eQTL_PIPELINE/testfile.txt"

params.genotype_source_functions="${baseDir}/R/genotype_functions/genotype_functions.r"
params.pseudobulk_source_functions="${baseDir}/R/expression_functions/pseudobulk_functions.r"
params.eqtl_source_functions="${baseDir}/R/MatrixEQTL_functions/matrixeqtl_source.r"

params.single_cell_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/roche_ms_decontx.rds"



process create_genotype {

    publishDir "${params.outdir}/", mode: "copy"

    input:
    path gds_file 

    output:
    path "genotype_012mat.csv", emit: genotype_mat
    path "snp_chromlocations.csv", emit: snp_chromlocations
    path "MAF_mat.csv", emit: maf_mat


    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    source("${params.genotype_source_functions}")
    generate_genotype_matrix(gds_file="$gds_file")


    """

}


process pseudobulk_singlecell{

   publishDir "${params.outdir}/", mode: "copy"

   input: 
   path single_cell_file

   output:
   path "*pseudobulk.csv", emit: pseudobulk_counts
   path "gene_locations.csv", emit: gene_locations

   script:
    """
    #!/usr/bin/env Rscript
    
    library(Seurat)
    library(BPCells)
    library(dplyr)
    source("${params.pseudobulk_source_functions}")

    seuratobj=readRDS("$single_cell_file")
    celltypelist=Seurat::SplitObject(seuratobj,split.by="CellType")

    aggregated_counts_list=pseudobulk_counts(celltypelist,
    min.cells=10,
    indiv_col="Individual_ID",
    assay="decontXcounts")
    
    gene_locations=get_gene_locations(aggregated_counts_list[[1]])
    data.table::fwrite(gene_locations,"gene_locations.csv")


    """

}


process run_matrixeQTL{
    
    input:
    path genotype_mat
    path snp_locations
    path expression_mat
    path gene_locations 

    output:
    path "*"


    script:
    """
    #!/usr/bin/env Rscript
    
    source("${params.eqtl_source_functions}")
    
    pseudobulk_data=fread("$expression_mat")
    genotype_data=fread("$genotype_mat")
    snp_locations=fread("$snp_locations")
    gene_locations=fread("$gene_locations")



    """

}

process qc_expression{
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path pseudobulk_file
    val min_expression

    output:
    path "*pseudobulk_normalised.csv" , emit: pseudobulk_normalised

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    pseudobulk_data <- fread("$pseudobulk_file")

    min_percentage <- $min_expression
    min_individuals <- min_percentage * ncol(pseudobulk_data)
    pseudobulk_data <- pseudobulk_data[rowSums(pseudobulk_data > 0) >= min_individuals, ]

    pseudobulk_data=log2(edgeR::cpm(pseudobulk_data)+1)
    fwrite(pseudobulk_data,"pseudobulk_normalised.csv")



    """
}

workflow{

  println """
    ========================================
    Welcome to the Nextflow eQTL pipeline
    ========================================
    Output Directory: ${params.outdir}
    GDS File: ${params.gds_file}
    Input Seurat File: ${params.single_cell_file}
    WorkDir: ${workflow.workDir}
    ========================================

    !WARNING - This pipeline is still in development and may not work as expected!

    """
    create_genotype(gds_file=params.gds_file)
    pseudobulk_singlecell(single_cell_file=params.single_cell_file)
    find_top_genes(pseudobulk_file=pseudobulk_singlecell.out.pseudobulk_counts.flatten())


    

}

workflow_expression{

    //aggregate counts
    pseudobulk_singlecell(single_cell_file=params.single_cell_file)

    //QC and normalisation
    find_top_genes(pseudobulk_file=pseudobulk_singlecell.out.pseudobulk_counts.flatten())
    

}

workflow.onComplete {

    println """
    ========================================
    Pipeline Completed!! Hello!
    ========================================

    """
}

