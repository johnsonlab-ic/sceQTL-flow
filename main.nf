nextflow.enable.dsl=2

params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds"

params.genotype_source_functions="${baseDir}/R/genotype_functions/genotype_functions.r"
params.pseudobulk_source_functions="${baseDir}/R/expression_functions/pseudobulk_functions.r"
params.eqtl_source_functions="${baseDir}/R/MatrixEQTL_functions/matrixeqtl_source.r"
params.quarto_report="${baseDir}/R/quarto_reports/run_report.qmd"

params.single_cell_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/roche_ms_decontx.rds"


/// Expression metrics
params.counts_assay="RNA"
params.counts_slot="counts"
params.celltype_column="CellType"
params.individual_column="Individual_ID"
params.min_cells=10
params.min_expression=0.05

// eQTL parameters
params.cis_distance=1e6
params.fdr_threshold=0.05
params.optimize_pcs=true


process create_genotype {

    label "process_high_memory"

    publishDir "${params.outdir}/genotype_files/", mode: "copy"

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
    generate_genotype_matrix(gds_file="$gds_file", 
    chain_fpath="/usr/local/src/hg19ToHg38.over.chain")


    """

}


process pseudobulk_singlecell{

   label "process_high_memory"

   publishDir "${params.outdir}/expression_matrices/", mode: "copy"

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
    celltypelist=Seurat::SplitObject(seuratobj,split.by="${params.celltype_column}")

    aggregated_counts_list=pseudobulk_counts(celltypelist,
    min.cells=as.numeric(${params.min_cells}),
    indiv_col="${params.individual_column}",
    assay="${params.counts_assay}",
    slot="${params.counts_slot}")

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

    label "process_single"

    publishDir "${params.outdir}/expression_matrices/", mode: 'copy'

    input:
    path pseudobulk_file

    output:
    path "*pseudobulk_normalised.csv" , emit: pseudobulk_normalised

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(dplyr)

    pseudobulk_data <- fread("$pseudobulk_file")
    pseudobulk_data <- pseudobulk_data %>% tibble::column_to_rownames(var="geneid")

    min_percentage <- as.numeric(${params.min_expression})
    min_individuals <- min_percentage * ncol(pseudobulk_data)
    pseudobulk_data <- pseudobulk_data[rowSums(pseudobulk_data > 0) >= min_individuals, ]

    cell_type_name <- gsub("_pseudobulk.csv", "", "$pseudobulk_file")

    pseudobulk_data=log2(edgeR::cpm(pseudobulk_data)+1) %>% as.data.frame()
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

    label "process_high_memory"

    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path genotype_mat
    path snp_locations
    path expression_mat
    path gene_locations 

    output:
    path "*_cis_MatrixEQTLout.rds", emit: eqtl_results
    path "*"


    script:
    """
    #!/usr/bin/env Rscript
    
    source("${params.eqtl_source_functions}")

    library(data.table)
    library(dplyr)
    
    ##read in data
    exp_mat=fread("$expression_mat") %>% tibble::column_to_rownames(var="geneid")
    geno_mat=fread("$genotype_mat") %>% tibble::column_to_rownames(var="snp")
    geno_loc=fread("$snp_locations")
    exp_loc=fread("$gene_locations")
    celltype=gsub("_pseudobulk_normalised.csv","","$expression_mat")


    ##keep same samples
    common_samples <- intersect(colnames(exp_mat), colnames(geno_mat))
    exp_mat <- exp_mat %>% select(all_of(common_samples))
    geno_mat <- geno_mat %>% select(all_of(common_samples))

    #harmonise gene_loc and snp_loc
    common_genes <- intersect(exp_loc %>% pull(geneid), rownames(exp_mat))
    exp_mat <- exp_mat %>% filter(rownames(exp_mat) %in% common_genes)

    
    geno_loc<-geno_loc[,c("annot","chrom","position")] %>% tibble::column_to_rownames(var="annot")
    geno_mat<-geno_mat[rownames(geno_loc),]
    geno_mat<-geno_mat[complete.cases(geno_mat),]
    geno_loc<-geno_loc[rownames(geno_mat),]
    geno_loc=geno_loc %>% mutate(annot=rownames(geno_loc)) %>% select(annot,chrom,position)



    calculate_ciseqtl(exp_mat=exp_mat,
    exp_loc=exp_loc,
    geno_mat=geno_mat,
    geno_loc=geno_loc,
    name=celltype,
    cisDist=${params.cis_distance},
    optimize_pcs=as.logical("${params.optimize_pcs}"))

    """

}

process combine_eqtls{
    label "process_high"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path eqtls

    output:
    path "mateqtlouts.rds"
    path "mateqtlouts_FDR_filtered.rds"


    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(dplyr)

    eqtls=as.character("$eqtls")
    eqtls=unlist(strsplit(eqtls, " "))
    celltypes=gsub("_cis_MatrixEQTLout.rds","",eqtls)

    #first, create a list of all the eqtl results at 5% FDR 
    eqtl_list=lapply(eqtls, function(x) {
        eqtl=as.data.frame(readRDS(x))
        eqtl=eqtl %>% filter(FDR<=${params.fdr_threshold})
        eqtl
    })
    names(eqtl_list)=celltypes
    saveRDS(eqtl_list, "mateqtlouts_FDR_filtered.rds")

    #now, create a list of all associations without cutoff
    eqtl_list=lapply(eqtls, function(x) {
        eqtl=as.data.frame(readRDS(x))
        eqtl
    })
    names(eqtl_list)=celltypes
    saveRDS(eqtl_list, "mateqtlouts.rds")


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

    ==============================================

                RUN PARAMETERS

    Expression/Pseudobulking parameters:

    Min cells for pseudobulking: ${params.min_cells}
    Min percentage for genes: ${params.min_expression}
    Cell-type column: ${params.celltype_column}
    Individual column: ${params.individual_column}
    Assay used: ${params.counts_assay}
    Slot used: ${params.counts_slot}

    eQTL parameters:

    Cis distance: ${params.cis_distance}
    FDR threshold: ${params.fdr_threshold}
    Optimize PCs: ${params.optimize_pcs}

    If "optimize PCs" is set to TRUE, the pipeline will run longer.

    ==============================================

    This is the stable pipeline. 

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
    
    combine_eqtls(eqtls= run_matrixeQTL.out.eqtl_results.collect())


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

