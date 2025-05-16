process process_covariates{

    label "process_low_memory"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path covmat
    tuple val(celltype), path(egenes_files)

    output:
    tuple val(celltype), path("${celltype}_covmat.txt"), emit: optimal_pcs


    script:
    """
    #!/usr/bin/env Rscript

    # Combine all files into a single data frame
    cov_file="$cov_file"
    if(file.size(cov_file) > 0){
        covmat=read.table(covmat, header=TRUE, row.names=1)
    }else{
        covmat = NULL
    }
    


    """



}