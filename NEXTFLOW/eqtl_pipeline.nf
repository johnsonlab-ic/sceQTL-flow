nextflow.enable.dsl=2

params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/test_geno.gds"
params.inputfile="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/eqtl_pipeline_dev/eQTL_PIPELINE/testfile.txt"
params.local=true


params.genotype_source_functions="${baseDir}/../genotype_functions/genotype_functions.r"
params.pseudobulk_source_functions="${baseDir}/../expression_functions/pseudobulk_functions.r"

params.single_cell_file="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/TEST_DATA/roche_ms_decontx.rds"


// process create_genotype_qsub {

//     publishDir "${params.outdir}/", mode: "copy"
//     executor="pbspro"
//     clusterOptions = "-lselect=1:ncpus=20:mem=240gb -l walltime=04:00:00"

//     input:
//     path gds_file
//     path genotype_source_functions

//     output:
//     path "genotype_012mat.csv"
//     path "snp_chromlocations.csv"
//     path "MAF_mat.csv"


//     script:
//     """
//     #!/usr/bin/env Rscript
//     library(dplyr)
//     source("$genotype_source_functions")
//     generate_genotype_matrix(gds_file="$gds_file")


//     """


// }

process create_genotype {

    publishDir "${params.outdir}/", mode: "copy"

    input:
    path gds_file 
    path genotype_source_functions

    output:
    path "genotype_012mat.csv", emit: genotype_mat
    path "snp_chromlocations.csv", emit: snp_chromlocations
    path "MAF_mat.csv", emit: maf_mat


    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    source("$genotype_source_functions")
    generate_genotype_matrix(gds_file="$gds_file")


    """

}

process count_snps{

    input:
    path genotype_mat
    

    output:
    stdout 

    script:
    """
    #!/usr/bin/env Rscript
    
    geno_mat=data.table::fread("$genotype_mat")
    print(nrow(geno_mat))


    """


}

process pseudobulk_singlecell{

   publishDir "${params.outdir}/", mode: "copy"

   input: 
   path single_cell_file
   path pseudobulk_source_functions

   output:
   path "*_aggregated_counts.csv"

   script:
    """
    #!/usr/bin/env Rscript
    
    library(Seurat)
    library(BPCells)
    source("$pseudobulk_source_functions")

    seuratobj=readRDS("$single_cell_file")
    celltypelist=Seurat::SplitObject(seuratobj,split.by=celltype_colname)

    aggregated_counts_list=pseudobulk_counts(celltypelist,
    min.cells=10,
    indiv_col="Individual_ID",
    assay="decontXcounts")

    for(i in 1:length(aggregated_counts_list)){
        write.table(aggregated_counts_list[[i]],paste0(names(aggregated_counts_list[i]),"_aggregated_counts.csv"))
    }


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
    Local Execution: ${params.local}
    ========================================
    """

    create_genotype(gds_file=params.gds_file,genotype_source_functions=params.genotype_source_functions)
    
    count_snps(genotype_mat=create_genotype.out.genotype_mat)


    pseudobulk_singlecell(single_cell_file=params.single_cell_file,pseudobulk_source_functions=params.pseudobulk_source_functions)

}