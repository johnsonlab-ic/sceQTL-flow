nextflow.enable.dsl=2

params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds"
params.inputfile="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/eqtl_pipeline_dev/eQTL_PIPELINE/testfile.txt"
params.email="ah3918@ic.ac.uk"

params.genotype_source_functions="${baseDir}/../R/genotype_functions/genotype_functions.r"
params.pseudobulk_source_functions="${baseDir}/../R/expression_functions/pseudobulk_functions.r"
params.eqtl_source_functions="${baseDir}/../R/MatrixEQTL_functions/matrixeqtl_source.r"

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
    
    aggregated_counts_list=lapply(aggregated_counts_list,function(x){
        x=log2(edgeR::cpm(x)+1)
        return(x)
    })

    for(i in 1:length(aggregated_counts_list)){
        data.table::fwrite(aggregated_counts_list[[i]],paste0(names(aggregated_counts_list[i]),"_pseudobulk.csv"))
    }

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

process find_top_genes {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path pseudobulk_file

    output:
    path "*top_genes_list.txt"

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    pseudobulk_data <- fread("$pseudobulk_file")
    gene_sums <- rowSums(pseudobulk_data[, -1, with=FALSE])
    top_genes <- head(order(gene_sums, decreasing=TRUE), 10)
    celltype <- gsub("_pseudobulk.csv", "", basename("$pseudobulk_file"))
    write.table(data.frame(celltype=celltype, top_genes=top_genes), 
    file=paste0(celltype,"top_genes_list.txt"), row.names=FALSE, col.names=TRUE, sep="\t", append=TRUE)
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

    This is the stable pipeline. 

    """
    create_genotype(gds_file=params.gds_file)
    pseudobulk_singlecell(single_cell_file=params.single_cell_file)
    find_top_genes(pseudobulk_file=pseudobulk_singlecell.out.pseudobulk_counts.flatten())


    

}


workflow.onComplete {

    println """
    ========================================
    Pipeline Completed!! Hello!
    ========================================

    """
}

