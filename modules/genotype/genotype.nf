process create_genotype {
    tag "${gds_file}"
    label "process_high_memory"
    publishDir "${params.outdir}/genotype_files/", mode: "copy"

    input:
    path gds_file
    path source_R

    output:
    path "genotype_012mat.csv", emit: genotype_mat
    path "snp_chromlocations.csv", emit: snp_chromlocations
    path "MAF_mat.csv", emit: maf_mat
    
    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    source("$source_R")
    generate_genotype_matrix(gds_file="$gds_file", 
    chain_fpath="/usr/local/src/hg19ToHg38.over.chain")
    """
}

process qc_genotype {

    input:
    path genotype_mat
    path snp_chromlocations
    val filter_chr

    output:
    path "qc_genotype_mat.csv", emit: qc_genotype_mat
    path "qc_snp_chromlocations.csv", emit: qc_snp_chromlocations

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)

    geno_mat = data.table::fread("$genotype_mat")
    snp_chromlocations = data.table::fread("$snp_chromlocations")
    
    # Optional filtering based on params.filter_chr
    if (!("$filter_chr" == "all")) {
            message("Filtering genotype matrix and snp_chromlocations by chromosome $filter_chr")
            filter_chr_vector = unlist(strsplit("$filter_chr", ","))
            snp_chromlocations = snp_chromlocations %>% filter(chrom %in% filter_chr_vector)
            geno_mat = geno_mat %>% filter(snp %in% snp_chromlocations\$annot)
        }
    data.table::fwrite(geno_mat, "qc_genotype_mat.csv")
    data.table::fwrite(snp_chromlocations, "qc_snp_chromlocations.csv")
    """
}
