#!/usr/bin/env Rscript
library(Seurat)
library(BPCells)
library(dplyr)

source("pseudobulk_functions.r")

seuratobj=readRDS("roche_ms_decontx.rds")
celltypelist=Seurat::SplitObject(seuratobj,split.by="CellType")
aggregated_counts_list=pseudobulk_counts(celltypelist,
min.cells=as.numeric(5),
indiv_col="Individual_ID",
assay="decontXcounts",
slot="counts")
for (i in 1:length(aggregated_counts_list)) {
    df=aggregated_counts_list[[i]] %>% mutate(geneid=row.names(.)) 
    celltype=names(aggregated_counts_list[i])
    data.table::fwrite(df, paste0(names(aggregated_counts_list[i]), "_pseudobulk.csv"))
}
gene_locations=get_gene_locations(aggregated_counts_list[[1]])
data.table::fwrite(gene_locations,"gene_locations.csv")
writeLines(names(aggregated_counts_list), "ct_names.txt")
