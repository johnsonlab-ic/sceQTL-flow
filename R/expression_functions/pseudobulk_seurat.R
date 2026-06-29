#!/usr/bin/env Rscript
# Native per-file pseudobulk for one Seurat .rds object.
# Same partial contract as pseudobulk_anndata.py: grouping (celltype + individual)
# comes from an external cell-level metadata file (--cell-metadata, keyed by
# --id-col) OR the object's own meta.data. Counts are summed cells->individuals
# within each cell type via a single sparse matmul.
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

# resolve grouping table for cells present in the label source
if (!is.null(cell_meta) && nchar(cell_meta) > 0 && cell_meta != "NO_FILE") {
  meta <- as.data.frame(fread(cell_meta))
  meta <- meta[!duplicated(meta[[id_col]]), ]
  rownames(meta) <- as.character(meta[[id_col]])
  present <- cell_ids[cell_ids %in% rownames(meta)]
  ndrop <- length(cell_ids) - length(present)
  if (ndrop > 0) message(sprintf("[PB] %d/%d cells not in metadata; dropped", ndrop, length(cell_ids)))
  if (length(present) == 0) stop("no cells matched between counts and metadata")
  ct  <- as.character(meta[present, celltype_col])
  ind <- as.character(meta[present, indiv_col])
  counts <- counts[, present, drop = FALSE]
} else {
  md <- obj@meta.data
  ct  <- as.character(md[cell_ids, celltype_col])
  ind <- as.character(md[cell_ids, indiv_col])
}

ok <- !is.na(ct) & !is.na(ind)
if (any(!ok)) counts <- counts[, ok, drop = FALSE]
ct <- ct[ok]; ind <- ind[ok]

grp <- factor(paste(ct, ind, sep = SEP))
design_t <- Matrix::fac2sparse(grp)          # groups x cells
summed <- as.matrix(counts %*% Matrix::t(design_t))  # genes x groups
storage.mode(summed) <- "integer"
ncells <- as.integer(Matrix::rowSums(design_t))
gl <- levels(grp)
g_ct  <- sub(paste0(SEP, ".*$"), "", gl)
g_ind <- sub(paste0("^.*", SEP), "", gl)
genes <- rownames(counts)

san <- function(s) gsub("[^A-Za-z0-9._+-]+", "_", s)
for (c in unique(g_ct)) {
  m <- g_ct == c
  mat <- summed[, m, drop = FALSE]
  colnames(mat) <- g_ind[m]
  df <- data.frame(geneid = genes, mat, check.names = FALSE)
  fwrite(df, file.path(outdir, paste0(san(c), "__", tag, "_partial.csv")))
  fwrite(data.frame(individual = g_ind[m], n_cells = ncells[m]),
         file.path(outdir, paste0(san(c), "__", tag, "_ncells.csv")))
  message(sprintf("[PB] %s: %d genes x %d individuals (tag=%s)", c, nrow(mat), ncol(mat), tag))
}
