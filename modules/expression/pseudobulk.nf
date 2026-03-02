process pseudobulk_singlecell {
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
    library(BPCells)
    library(dplyr)

    source("$source_R")

    seuratobj=readRDS("$single_cell_file")

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
    # Safe counts preview (Seurat v5 layer-aware)
    counts_preview <- tryCatch({
        assay_obj <- seuratobj[["${params.counts_assay}"]]
        all_layers <- tryCatch(SeuratObject::Layers(assay_obj), error=function(e) NULL)
        if (is.null(all_layers)) {
            Seurat::GetAssayData(seuratobj, slot="${params.counts_slot}")
        } else if ("${params.counts_slot}" %in% all_layers) {
            SeuratObject::LayerData(assay_obj, layer="${params.counts_slot}")
        } else {
            slot_layers <- grep(paste0("^", "${params.counts_slot}", "\\\\."), all_layers, value=TRUE)
            if (length(slot_layers) > 0) {
                tmp_obj <- SeuratObject::JoinLayers(seuratobj, assay="${params.counts_assay}", layers=slot_layers, new="${params.counts_slot}")
                SeuratObject::LayerData(tmp_obj[["${params.counts_assay}"]], layer="${params.counts_slot}")
            } else {
                stop(paste0("Could not find layer '", "${params.counts_slot}", "' in assay '", "${params.counts_assay}", "'."))
            }
        }
    }, error=function(e) {
        cat("[PB] Counts preview retrieval error:", conditionMessage(e), "\n"); NULL
    })
    if (!is.null(counts_preview)) {
        cat("[PB] Counts preview dims:", nrow(counts_preview), "genes x", ncol(counts_preview), "cells\n")
    }

    celltypelist=Seurat::SplitObject(seuratobj,split.by="${params.celltype_column}")
    aggregated_counts_list=pseudobulk_counts(celltypelist,
    min.cells=as.numeric(${params.min_cells}),
    indiv_col="${params.individual_column}",
    assay="${params.counts_assay}",
    slot="${params.counts_slot}")
    for (i in 1:length(aggregated_counts_list)) {
        df=aggregated_counts_list[[i]] %>% mutate(geneid=row.names(.)) 
        celltype=names(aggregated_counts_list[i])
        celltype_sanitized=sanitize_celltype_name(celltype)
        data.table::fwrite(df, paste0(celltype_sanitized, "_pseudobulk.csv"))
    }
    gene_locations=get_gene_locations(aggregated_counts_list[[1]])
    data.table::fwrite(gene_locations,"gene_locations.csv")
    sanitized_names <- sapply(names(aggregated_counts_list), sanitize_celltype_name)
    writeLines(sanitized_names, "ct_names.txt")
    """
}
