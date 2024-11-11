nextflow.enable.dsl=2

params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/roche/live/ALEX//PROCESSED_DATA/PROCESSED_GENOTYPE/FINAL/final_geno_440samples.gds"
params.inputfile="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/eqtl_pipeline_dev/eQTL_PIPELINE/testfile.txt"
params.local=true
params.genotype_source_functions="${baseDir}/../genotype_functions/genotype_functions.r"



process create_genotype_qsub {

    publishDir "${params.outdir}/", mode: "copy"
    executor="pbspro"
    clusterOptions = "-lselect=1:ncpus=20:mem=240gb -l walltime=04:00:00"

    input:
    path gds_file
    path genotype_source_functions

    output:
    path "*"


    script:
    """
    #!/usr/bin/env Rscript
    
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
    path "*"


    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    source("$genotype_source_functions")
    generate_genotype_matrix(gds_file="$gds_file")


    """

}

workflow{

    if(params.local){
        create_genotype(gds_file=params.gds_file,genotype_source_functions=params.genotype_source_functions)
    }else{
        create_genotype_qsub(gds_file=params.gds_file,genotype_source_functions=params.genotype_source_functions)
    }
    
}