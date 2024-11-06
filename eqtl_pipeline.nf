

params.outdir="./"
params.gds_file=""
singularity.enabled = true
params.container = "/rds/general/user/ah3918/home/CONTAINERS/genotype_cont_latest.sif"

process create_genotype{

    publishDir "${params.outdir}/MatrixEQTL_IO", mode: "copy"

    input:
    path genofile

    output:
    path "genotype_012mat.csv", emit: genomat
    path "snp_chromlocations.csv", emit: snplocs 
    path "MAF_mat.csv", emit: mafmat

    container "params.container"

   script:
    """
    Rscript -e 'source("${baseDir}/genotype_functions.r"); get_genotype_matrix(gds_file="${params.genofile}")'
    """


}


workflow{
    create_genotype()

}