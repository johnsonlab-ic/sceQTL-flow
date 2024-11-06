

params.outdir="./"
params.gds_file=""


process create_genotype{

    // publishDir "${params.outdir}/MatrixEQTL_IO", mode: "copy"

    // input:
    // path genofile

    // output:
    // path "genotype_012mat.csv", emit: genomat
    // path "snp_chromlocations.csv", emit: snplocs 
    // path "MAF_mat.csv", emit: mafmat

    // container "${baseDir}.container"

   script:
    // """
    // Rscript -e 'source("${baseDir}/../genotype_functions/genotype_functions.r"); get_genotype_matrix(gds_file="${params.genofile}")'
    // """
    """
    Rscript -e 'library(SeqArray)'
    """

}


workflow{
    create_genotype()

}