process run_matrixeQTL {

  
    tag "${expression_mat}"
    label "process_low"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path genotype_mat
    path expression_mat
    path cov_file

    output:
    path "*_cis_MatrixEQTLout.rds", emit: eqtl_results
    path "*"

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(dplyr)
    


    """
}
