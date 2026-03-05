process combine_eqtls {

    label "process_high"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path eqtls
    path maf_mat

    output:
    path "mateqtlouts.rds", emit: mateqtlouts
    path "mateqtlouts_FDR_filtered.rds", emit: mateqtlouts_FDR_filtered
    path "genes_tested.csv", emit: genes_tested
    
    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(dplyr)

    maf_df <- fread("$maf_mat", data.table = FALSE)
    required_maf_cols <- c("ref", "alt", "snp")
    missing_maf_cols <- setdiff(required_maf_cols, colnames(maf_df))
    if (length(missing_maf_cols) > 0) {
        stop(paste0("MAF file is missing required columns: ", paste(missing_maf_cols, collapse = ", ")))
    }
    maf_df <- maf_df %>%
        select(snp, ref, alt) %>%
        distinct(snp, .keep_all = TRUE)

    eqtls <- list.files(pattern = "_cis_MatrixEQTLout\\\\.rds", full.names = TRUE)
    eqtls <- sort(eqtls)
    if (length(eqtls) == 0) {
        stop("No *_cis_MatrixEQTLout.rds files found in work directory.")
    }
    celltypes <- gsub("_cis_MatrixEQTLout.rds", "", basename(eqtls))
    eqtl_list=lapply(eqtls, function(x) {
        eqtl=as.data.frame(readRDS(x))
        eqtl=eqtl %>%
            left_join(maf_df, by = c("SNP" = "snp")) %>%
            rename(effect_allele = ref, other_allele = alt) %>%
            filter(FDR<=${params.fdr_threshold})
        eqtl
    })
    names(eqtl_list)=celltypes
    saveRDS(eqtl_list, "mateqtlouts_FDR_filtered.rds")
    eqtl_list=lapply(eqtls, function(x) {
        eqtl=as.data.frame(readRDS(x))
        eqtl=eqtl %>%
            left_join(maf_df, by = c("SNP" = "snp")) %>%
            rename(effect_allele = ref, other_allele = alt)
        eqtl
    })
    names(eqtl_list)=celltypes
    saveRDS(eqtl_list, "mateqtlouts.rds")

    genes_tested <- tibble(
        celltype = names(eqtl_list),
        genes_tested = vapply(eqtl_list, function(df) dplyr::n_distinct(df\$gene), numeric(1))
    )
    fwrite(genes_tested, "genes_tested.csv")
    """
}
