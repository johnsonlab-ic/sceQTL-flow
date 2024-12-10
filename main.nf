nextflow.enable.dsl=2

params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds"
params.inputfile="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/eqtl_pipeline_dev/eQTL_PIPELINE/testfile.txt"

params.genotype_source_functions="${baseDir}/R/genotype_functions/genotype_functions.r"
params.pseudobulk_source_functions="${baseDir}/R/expression_functions/pseudobulk_functions.r"
params.eqtl_source_functions="${baseDir}/R/MatrixEQTL_functions/matrixeqtl_source.r"
params.quarto_report="${baseDir}/R/quarto_reports/run_report.qmd"

params.single_cell_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/roche_ms_decontx.rds"


/// Expression metrics
params.min_cells=10
params.min_expression=0.1


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
    min.cells=as.numeric(${params.min_cells}),
    indiv_col="Individual_ID",
    assay="decontXcounts")

    for (i in 1:length(aggregated_counts_list)) {

        df=aggregated_counts_list[[i]] %>% mutate(geneid=row.names(.)) 
        
        # Write the data frame to a CSV file
        data.table::fwrite(df, paste0(names(aggregated_counts_list[i]), "_pseudobulk.csv"))
    }

    gene_locations=get_gene_locations(aggregated_counts_list[[1]])
    data.table::fwrite(gene_locations,"gene_locations.csv")


    """

}




process qc_expression{
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path pseudobulk_file

    output:
    path "*pseudobulk_normalised.csv" , emit: pseudobulk_normalised

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(tibble)
    pseudobulk_data <- fread("$pseudobulk_file")
    pseudobulk_data <- pseudobulk_data %>% column_to_rownames(var="geneid")

    min_percentage <- as.numeric(${params.min_expression})
    min_individuals <- min_percentage * ncol(pseudobulk_data)
    pseudobulk_data <- pseudobulk_data[rowSums(pseudobulk_data > 0) >= min_individuals, ]

    cell_type_name <- gsub("_pseudobulk.csv", "", "$pseudobulk_file")

    pseudobulk_data=log2(edgeR::cpm(pseudobulk_data)+1)
    pseudobulk_data = pseudobulk_data %>% mutate(geneid=row.names(.))

    # Save the normalized data
    fwrite(pseudobulk_data, paste0(cell_type_name, "_pseudobulk_normalised.csv"))
    """
}


process qc_genotype {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path genotype_mat
    path snp_locations

    output:
    path "*"

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)

    genotype_data <- fread("$genotype_mat")
    snp_locations <- fread("$snp_locations")


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

    library(data.table)
    library(dplyr)
    
    exp_mat=fread("$expression_mat")
    geno_mat=fread("$genotype_mat")
    geno_loc=fread("$snp_locations")
    exp_loc=fread("$gene_locations")

    common_samples <- intersect(colnames(exp_mat), colnames(geno_mat))
    exp_mat <- exp_mat %>% select(all_of(common_samples))
    geno_mat <- geno_mat %>% select(all_of(common_samples))

    common_genes <- intersect(exp_loc %>% pull(geneid), rownames(exp_mat))

    # Subset the data frame using dplyr's filter function
    exp_mat <- exp_mat %>% filter(rownames(exp_mat) %in% common_genes)

    calculate_ciseqtl(exp_mat=exp_mat,
    exp_loc=exp_loc,
    geno_mat=geno_mat,
    geno_loc=geno_loc)

    """

}


process final_report{

    publishDir "${params.outdir}", mode: 'copy'

    input: 
    path pseudobulk_file_list
    path genotype_file
    path report_file

    output: 

    path "report*"


    script:


    """
    #!/bin/bash

    quarto render $report_file --output-dir ./ \
    -P genotype_file:$genotype_file

    """

}

workflow{

  println """
    ========================================
    Welcome to the Nextflow eQTL pipeline
    ========================================
    
    File inputs:

    Output Directory: ${params.outdir}
    GDS File: ${params.gds_file}
    Input Seurat File: ${params.single_cell_file}
    WorkDir: ${workflow.workDir}

    ========================================

    Expression QC metrics:

    Min cells for pseudobulking: ${params.min_cells}
    Min percentage for genes: ${params.min_expression}

    ========================================

    !WARNING - This pipeline is still in development and may not work as expected!

    """
    create_genotype(gds_file= params.gds_file)

    //aggregate counts
    pseudobulk_singlecell(single_cell_file= params.single_cell_file)

    //QC and normalisation
    qc_expression(pseudobulk_file= pseudobulk_singlecell.out.pseudobulk_counts.flatten())

    //run matrix eQTL
    run_matrixeQTL(
        genotype_mat= create_genotype.out.genotype_mat,
        snp_locations= create_genotype.out.snp_chromlocations,
        expression_mat= qc_expression.out.pseudobulk_normalised.flatten(),
        gene_locations= pseudobulk_singlecell.out.gene_locations
    )
    
    // final_report(
    //     pseudobulk_file_list= qc_expression.out.collect(),
    //     genotype_file= create_genotype.out.genotype_mat,
    //     report_file=params.quarto_report
    // )



}


workflow.onComplete {

    println """
    ========================================
    Pipeline Completed!
    ========================================

    """
}

