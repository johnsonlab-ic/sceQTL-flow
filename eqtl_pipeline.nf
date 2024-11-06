

params.outdir="./"
params.gds_file=""


process create_genotype{

    publishDir "${params.outdir}/MatrixEQTL_IO", mode: "copy"

    input:
    path genofile

    output:
    path "genotype_012mat.csv", emit: genomat
    path "snp_chromlocations.csv", emit: snplocs 
    path "MAF_mat.csv", emit: mafmat



    """
    #!/rds/general/user/ah3918/home/anaconda3/envs/OSIRIS/bin/Rscript

    source("/rds/general/user/ah3918/projects/roche/live/ALEX/SCRIPTS/GITHUB/eQTL_scripts/nextflow_pipelines/singlecell/qtl_pipeline/cellQTL_source.r")
    
    get_genotype_matrix(gds_file="${params.genofile}")



    """


}


workflow{
    create_genotype()

}