#!/usr/bin/env Rscript
# Native per-file pseudobulk for one Seurat .rds object.
# Same partial contract as pseudobulk_anndata.py: --celltype-col accepts a
# comma-separated list of columns to pseudobulk on independently. Individual/
# celltype labels come from an external cell-level metadata file
# (--cell-metadata, keyed by --id-col) OR the object's own meta.data. Counts
# are summed cells->individuals within each category via a single sparse
# matmul. Output files are prefixed "<column>-<category>__..." so categories
# from different columns can't collide (and get silently summed together)
# downstream.
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i)) args[i + 1] else default
}

rds          <- getarg("--rds")
celltype_col <- getarg("--celltype-col")
indiv_col    <- getarg("--indiv-col")
counts_assay <- getarg("--counts-assay", "RNA")
counts_slot  <- getarg("--counts-slot", "counts")
cell_meta    <- getarg("--cell-metadata", NULL)
id_col       <- getarg("--id-col", "cell_id")
outdir       <- getarg("--outdir", ".")
tag          <- getarg("--tag", sub("\\.rds$", "", basename(rds), ignore.case = TRUE))
SEP          <- "\x1f"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(rds)
DefaultAssay(obj) <- counts_assay

assay_obj <- obj[[counts_assay]]
layers <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) NULL)
if (!is.null(layers) && length(layers) > 1) {
  message(sprintf("[PB] joining %d layers in assay %s", length(layers), counts_assay))
  obj <- SeuratObject::JoinLayers(obj, assay = counts_assay)
}
counts <- tryCatch(
  SeuratObject::GetAssayData(obj, assay = counts_assay, layer = counts_slot),
  error = function(e) Seurat::GetAssayData(obj, slot = counts_slot)
)  # genes x cells
cell_ids <- colnames(counts)

celltype_cols <- trimws(strsplit(celltype_col, ",")[[1]])
celltype_cols <- celltype_cols[nzchar(celltype_cols)]

# resolve grouping table for cells present in the label source
if (!is.null(cell_meta) && nchar(cell_meta) > 0 && cell_meta != "NO_FILE") {
  meta <- as.data.frame(fread(cell_meta))
  meta <- meta[!duplicated(meta[[id_col]]), ]
  rownames(meta) <- as.character(meta[[id_col]])
  present <- cell_ids[cell_ids %in% rownames(meta)]
  ndrop <- length(cell_ids) - length(present)
  if (ndrop > 0) message(sprintf("[PB] %d/%d cells not in metadata; dropped", ndrop, length(cell_ids)))
  if (length(present) == 0) stop("no cells matched between counts and metadata")
  ind <- as.character(meta[present, indiv_col])
  meta_sub <- meta[present, celltype_cols, drop = FALSE]
  counts <- counts[, present, drop = FALSE]
} else {
  md <- obj@meta.data
  ind <- as.character(md[cell_ids, indiv_col])
  meta_sub <- md[cell_ids, celltype_cols, drop = FALSE]
}

ok_ind <- !is.na(ind)
if (any(!ok_ind)) {
  counts <- counts[, ok_ind, drop = FALSE]
  ind <- ind[ok_ind]
  meta_sub <- meta_sub[ok_ind, , drop = FALSE]
}

san <- function(s) gsub("[^A-Za-z0-9._+-]+", "_", s)

for (col in celltype_cols) {
  ct <- as.character(meta_sub[[col]])
  ok_ct <- !is.na(ct)
  n_invalid <- sum(!ok_ct)
  if (n_invalid > 0) message(sprintf("[PB] [%s] %d cells with missing celltype label; dropped", col, n_invalid))

  ct_c <- ct[ok_ct]
  ind_c <- ind[ok_ct]
  counts_c <- counts[, ok_ct, drop = FALSE]

  grp <- factor(paste(ct_c, ind_c, sep = SEP))
  design_t <- Matrix::fac2sparse(grp)                  # groups x cells
  summed <- as.matrix(counts_c %*% Matrix::t(design_t))  # genes x groups
  storage.mode(summed) <- "integer"
  ncells <- as.integer(Matrix::rowSums(design_t))
  gl <- levels(grp)
  g_ct  <- sub(paste0(SEP, ".*$"), "", gl)
  g_ind <- sub(paste0("^.*", SEP), "", gl)
  genes <- rownames(counts_c)

  for (c in unique(g_ct)) {
    m <- g_ct == c
    mat <- summed[, m, drop = FALSE]
    colnames(mat) <- g_ind[m]
    df <- data.frame(geneid = genes, mat, check.names = FALSE)
    label <- san(paste0(col, "-", c))
    fwrite(df, file.path(outdir, paste0(label, "__", tag, "_partial.csv")))
    fwrite(data.frame(individual = g_ind[m], n_cells = ncells[m]),
           file.path(outdir, paste0(label, "__", tag, "_ncells.csv")))
    message(sprintf("[PB] [%s] %s: %d genes x %d individuals (tag=%s)", col, c, nrow(mat), ncol(mat), tag))
  }
}
