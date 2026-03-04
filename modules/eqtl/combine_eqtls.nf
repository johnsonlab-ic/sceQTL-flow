process combine_eqtls {

    label "process_high"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path eqtls

    output:
    path "mateqtlouts.rds", emit: mateqtlouts
    path "mateqtlouts_FDR_filtered.rds", emit: mateqtlouts_FDR_filtered
    
    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(dplyr)
    eqtls <- list.files(pattern = "_cis_MatrixEQTLout\\\\.rds", full.names = TRUE)
    eqtls <- sort(eqtls)
    if (length(eqtls) == 0) {
        stop("No *_cis_MatrixEQTLout.rds files found in work directory.")
    }
    celltypes <- gsub("_cis_MatrixEQTLout.rds", "", basename(eqtls))
    eqtl_list=lapply(eqtls, function(x) {
        eqtl=as.data.frame(readRDS(x))
        eqtl=eqtl %>% filter(FDR<=${params.fdr_threshold})
        eqtl
    })
    names(eqtl_list)=celltypes
    saveRDS(eqtl_list, "mateqtlouts_FDR_filtered.rds")
    eqtl_list=lapply(eqtls, function(x) {
        eqtl=as.data.frame(readRDS(x))
        eqtl
    })
    names(eqtl_list)=celltypes
    saveRDS(eqtl_list, "mateqtlouts.rds")
    """
}
