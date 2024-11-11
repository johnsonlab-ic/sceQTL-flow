nextflow.enable.dsl=2

params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/roche/live/ALEX//PROCESSED_DATA/PROCESSED_GENOTYPE/FINAL/final_geno_440samples.gds"
params.inputfile="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/eqtl_pipeline_dev/eQTL_PIPELINE/testfile.txt"
params.local=true
params.genotype_source_functions="${baseDir}/../genotype_functions/genotype_functions.r"
params.single_cell_file=""


process create_genotype_qsub {

    publishDir "${params.outdir}/", mode: "copy"
    executor="pbspro"
    clusterOptions = "-lselect=1:ncpus=20:mem=240gb -l walltime=04:00:00"

    input:
    path gds_file
    path genotype_source_functions

    output:
    path "genotype_mat.csv"
    path "snp_chromlocations.csv"
    path "MAF_mat.csv"


    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    source("$genotype_source_functions")
    generate_genotype_matrix(gds_file="$gds_file")


    """


}

process create_genotype {

    publishDir "${params.outdir}/", mode: "copy"

    input:
    path gds_file 
    path genotype_source_functions

    output:
    path "genotype_mat.csv"
    path "snp_chromlocations.csv"
    path "MAF_mat.csv"


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

// process pseudobulk_singlecell{

//     input: rds_file from file(params.singlecell_file)

//    script:
//     """
//     #!/usr/bin/env Rscript
    
//     library(Seurat)


//     """

// }

workflow{

    if(params.local){
        create_genotype(gds_file=params.gds_file,genotype_source_functions=params.genotype_source_functions)
    }else{
        create_genotype_qsub(gds_file=params.gds_file,genotype_source_functions=params.genotype_source_functions)
    }
    
    count_snps(genotype_mat=create_genotype.genotype_mat)
}