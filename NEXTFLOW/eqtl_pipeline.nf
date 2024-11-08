

params.outdir="/rds/general/user/ah3918/projects/puklandmarkproject/ephemeral/tmp/"
params.gds_file="/rds/general/user/ah3918/projects/roche/live/ALEX//PROCESSED_DATA/PROCESSED_GENOTYPE/FINAL/final_geno_440samples.gds"


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

process print_dir {

    publishDir "${params.outdir}/", mode: "copy"

    output:
    path "dir_structure.txt"

    script:
    """
    ls -R > dir_structure.txt
    """
}

workflow{
    // create_genotype()
    print_dir()
}