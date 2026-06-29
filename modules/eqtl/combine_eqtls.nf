process combine_eqtls {

    label "process_low"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path sig_files
    path summary_files

    output:
    path "mateqtlouts_FDR_filtered.rds", emit: mateqtlouts_FDR_filtered
    path "eqtl_summary.rds", emit: summary
    path "eqtl_summary.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)

    # Already FDR-filtered upstream; just collate the small per-celltype tables.
    sig_files = unlist(strsplit("$sig_files", " "))
    celltypes = gsub("_cis_eqtl_sig.rds", "", basename(sig_files))
    eqtl_list = lapply(sig_files, function(x) as.data.frame(readRDS(x)))
    names(eqtl_list) = celltypes
    saveRDS(eqtl_list, "mateqtlouts_FDR_filtered.rds")

    summary_files = unlist(strsplit("$summary_files", " "))
    summ = do.call(rbind, lapply(summary_files, function(x) as.data.frame(readRDS(x))))
    saveRDS(summ, "eqtl_summary.rds")
    fwrite(summ, "eqtl_summary.csv")
    """
}
