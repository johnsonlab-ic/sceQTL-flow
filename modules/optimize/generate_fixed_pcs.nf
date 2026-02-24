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
    celltype = expression_file.baseName.replaceAll('_residuals$', '').replaceAll('_pseudobulk_normalised$', '')
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(dplyr)

    exp_mat <- fread("${expression_file}") %>% tibble::column_to_rownames(var = "geneid")

    pc_obj <- prcomp(t(exp_mat), scale. = TRUE)
    pcs <- pc_obj\$x[, 1:${n_pcs}]
    pcs <- as.data.frame(pcs)
    colnames(pcs) <- paste0("PC", 1:${n_pcs})
    pcs <- t(pcs)

    write.table(pcs, file = "${celltype}_${n_pcs}_pcs.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
    """
}
