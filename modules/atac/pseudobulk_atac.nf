process pseudobulk_atac {
    tag "${single_cell_file}"
    label "process_high_memory"
    publishDir "${params.outdir}/expression_matrices/", mode: "copy"

    input:
    path single_cell_file
    path source_R

    output:
    path "ct_names.txt", emit: ct_names
    path "*pseudobulk.csv", emit: pseudobulk_counts
    path "gene_locations.csv", emit: gene_locations

    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    library(dplyr)
    library(Matrix)
    library(data.table)

    source("$source_R")

    seuratobj = readRDS("$single_cell_file")

    # Diagnostics: print object structure, assays, default assay, and layers
    cat("[PB] ==== Seurat Object Diagnostics ====\n")
    cat("[PB] File:", "$single_cell_file", "\n")
    cat("[PB] Class:", paste(class(seuratobj), collapse=","), "\n")
    cat("[PB] Assays:", paste(names(seuratobj@assays), collapse=","), "\n")
    cat("[PB] Default assay:", Seurat::DefaultAssay(seuratobj), "\n")
    lyr_info <- tryCatch(SeuratObject::Layers(seuratobj[["${params.counts_assay}"]]), error=function(e) NULL)
    if (is.null(lyr_info)) {
        cat("[PB] Layers in assay ${params.counts_assay}:", "<none or not available>", "\n")
    } else {
        cat("[PB] Layers in assay ${params.counts_assay}:", paste(lyr_info, collapse=","), "\n")
    }
    cat("[PB] Cells:", ncol(seuratobj), "\n")

    counts_preview <- tryCatch(Seurat::GetAssayData(seuratobj, assay="${params.counts_assay}", slot="${params.counts_slot}"),
        error=function(e) {
            cat("[PB] Counts preview retrieval error:", conditionMessage(e), "\n"); NULL
        })
    if (!is.null(counts_preview)) {
        cat("[PB] Counts preview dims:", nrow(counts_preview), "peaks x", ncol(counts_preview), "cells\n")
    }

    celltypelist = Seurat::SplitObject(seuratobj, split.by="${params.celltype_column}")
    aggregated_counts_list = pseudobulk_peaks(celltypelist,
        min.cells=as.numeric(${params.min_cells}),
        indiv_col="${params.individual_column}",
        assay="${params.counts_assay}")

    for (i in 1:length(aggregated_counts_list)) {
        df = aggregated_counts_list[[i]] %>% mutate(geneid=row.names(.))
        celltype = names(aggregated_counts_list[i])
        data.table::fwrite(df, paste0(names(aggregated_counts_list[i]), "_pseudobulk.csv"))
    }

    gene_locations = get_peak_locations(rownames(aggregated_counts_list[[1]]))
    data.table::fwrite(gene_locations, "gene_locations.csv")
    writeLines(names(aggregated_counts_list), "ct_names.txt")
    """
}