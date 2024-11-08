nextflow.enable.dsl=2

params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/roche/live/ALEX//PROCESSED_DATA/PROCESSED_GENOTYPE/FINAL/final_geno_440samples.gds"
params.inputfile="/rds/general/user/ah3918/projects/puklandmarkproject/live/Users/Alex/pipelines/eqtl_pipeline_dev/eQTL_PIPELINE/testfile.txt"
// process create_genotype{

//     publishDir "${params.outdir}/MatrixEQTL_IO", mode: "copy"

//     // input:
//     // path genofile

//     output:
//     path "*"

//     // container "${baseDir}.container"
//    script:
//     """
//     Rscript -e 'source("..//genotype_functions/genotype_functions.r"); generate_genotype_matrix(gds_file="${params.gds_file}")'
//     """
//     // """
//     // Rscript -e 'library(GenomicRanges)'
//     // """

// }


// process double_file_length {

//     publishDir "${params.outdir}/", mode: "copy"

//     input:
//     path input_file

//     output:
//     path "doubled_${input_file.name}"

//     script:
//     """
//     cat ${input_file} ${input_file} > doubled_${input_file.name}
//     """
// }

process genotype {

    publishDir "${params.outdir}/", mode: "copy"

    input:
    path gds_file

    output:
    path "*"

    script:
    """
    #!/usr/bin/env R

    genofile=seqArray::seqOpen($gds_file)
    sample.id<-SeqArray::seqGetData(genofile,"sample.id")
    writeLines(sample.id,"samples.txt")

    """

}

workflow{
    // create_genotype()
    // double_file_length(input_file=params.inputfile)
    genotype(gds_file=params.gds_file)
}