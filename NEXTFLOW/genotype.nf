process create_genotype {
    label "process_high_memory"
    publishDir "${params.outdir}/genotype_files/", mode: "copy"
    input:
    path gds_file 
    output:
    path "genotype_012mat.csv", emit: genotype_mat
    path "snp_chromlocations.csv", emit: snp_chromlocations
    path "MAF_mat.csv", emit: maf_mat
    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    source("${params.genotype_source_functions}")
    generate_genotype_matrix(gds_file="$gds_file", 
    chain_fpath="/usr/local/src/hg19ToHg38.over.chain")
    """
}
