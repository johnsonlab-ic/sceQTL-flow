process pseudobulk_singlecell {
    label "process_high_memory"
    publishDir "${params.outdir}/expression_matrices/", mode: "copy"
    input: 
    path single_cell_file
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
    source("${params.pseudobulk_source_functions}")
    seuratobj=readRDS("$single_cell_file")
    celltypelist=Seurat::SplitObject(seuratobj,split.by="${params.celltype_column}")
    aggregated_counts_list=pseudobulk_counts(celltypelist,
    min.cells=as.numeric(${params.min_cells}),
    indiv_col="${params.individual_column}",
    assay="${params.counts_assay}",
    slot="${params.counts_slot}")
    for (i in 1:length(aggregated_counts_list)) {
        df=aggregated_counts_list[[i]] %>% mutate(geneid=row.names(.)) 
        celltype=names(aggregated_counts_list[i])
        data.table::fwrite(df, paste0(names(aggregated_counts_list[i]), "_pseudobulk.csv"))
    }
    gene_locations=get_gene_locations(aggregated_counts_list[[1]])
    data.table::fwrite(gene_locations,"gene_locations.csv")
    writeLines(names(aggregated_counts_list), "ct_names.txt")
    """
}
