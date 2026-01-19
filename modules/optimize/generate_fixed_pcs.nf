process generate_fixed_pcs {
    tag "Generate fixed PCs for ${expression_file}"
    label "process_low_memory"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    path expression_file
    val n_pcs

    output:
    tuple val(celltype), path("*_pcs.txt"), emit: exp_pcs

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(dplyr)

    exp_mat <- fread("${expression_file}") %>% tibble::column_to_rownames(var = "geneid")

    # Derive celltype from filename
    celltype <- basename("${expression_file}")
    celltype <- gsub("_residuals.csv", "", celltype)
    celltype <- gsub("_pseudobulk_normalised.csv", "", celltype)
    celltype <- gsub(".csv\$", "", celltype)

    # Compute PCs
    pcs <- prcomp(t(exp_mat), scale. = TRUE)$x[, 1:${n_pcs}]
    pcs <- as.data.frame(pcs)
    colnames(pcs) <- paste0("PC", 1:${n_pcs})
    pcs <- t(pcs)

    out_file <- paste0(celltype, "_", ${n_pcs}, "_pcs.txt")
    write.table(pcs, file = out_file, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
    """
}
