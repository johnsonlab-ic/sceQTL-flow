

params.outdir="./"
params.gds_file="/rds/general/user/ah3918/projects/roche/live/ALEX//PROCESSED_DATA/PROCESSED_GENOTYPE/FINAL/final_geno_440samples.gds"


process create_genotype{

    publishDir "${params.outdir}/MatrixEQTL_IO", mode: "copy"

    // input:
    // path genofile

    output:
    path "*"

    // container "${baseDir}.container"

   script:
    """
    Rscript -e 'source("/mnt/genotype_functions/genotype_functions.r"); get_genotype_matrix(gds_file="${params.genofile}")'
    """
    // """
    // Rscript -e 'library(GenomicRanges)'
    // """

}


workflow{
    create_genotype()
}