#!/usr/bin/env Rscript
# Combine per-file pseudobulk partials into final per-celltype matrices.
# Per cell type: intersect genes across files, concatenate individuals (summing
# any individual that appears in more than one file), drop individuals whose
# total cell count <= --min-cells. Writes <celltype>_pseudobulk.csv (geneid +
# one column per individual) plus a single gene_locations.csv.
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i)) args[i + 1] else default
}
indir     <- getarg("--indir", ".")
outdir    <- getarg("--outdir", ".")
min_cells <- as.numeric(getarg("--min-cells", "5"))
source_fn <- getarg("--source", NULL)
if (!is.null(source_fn)) source(source_fn)  # provides get_gene_locations()
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

partials <- list.files(indir, pattern = "__.*_partial\\.csv$", full.names = TRUE)
if (!length(partials)) stop("no *_partial.csv files found in ", indir)
ct_of <- function(p) sub("__.*$", "", basename(p))
celltypes <- unique(vapply(partials, ct_of, character(1)))

all_genes <- character(0)
for (ct in celltypes) {
  files <- partials[vapply(partials, ct_of, character(1)) == ct]

  ncell_tot <- numeric(0)
  mats <- list()
  genes_common <- NULL
  for (f in files) {
    d <- as.data.frame(fread(f))
    rn <- d[["geneid"]]; d[["geneid"]] <- NULL
    rownames(d) <- rn
    genes_common <- if (is.null(genes_common)) rn else intersect(genes_common, rn)
    mats[[f]] <- d
    nc <- as.data.frame(fread(sub("_partial\\.csv$", "_ncells.csv", f)))
    for (k in seq_len(nrow(nc))) {
      ind <- as.character(nc[["individual"]][k])
      prev <- ncell_tot[ind]
      ncell_tot[ind] <- (if (is.na(prev)) 0 else prev) + as.numeric(nc[["n_cells"]][k])
    }
  }

  combined <- mats[[files[1]]][genes_common, , drop = FALSE]
  if (length(files) > 1) {
    for (f in files[-1]) {
      d <- mats[[f]][genes_common, , drop = FALSE]
      for (cn in colnames(d)) {
        if (cn %in% colnames(combined)) combined[[cn]] <- combined[[cn]] + d[[cn]]
        else combined[[cn]] <- d[[cn]]
      }
    }
  }

  keep <- colnames(combined)[vapply(colnames(combined), function(cn) {
    tot <- ncell_tot[cn]; !is.na(tot) && tot > min_cells
  }, logical(1))]
  dropped <- setdiff(colnames(combined), keep)
  combined <- combined[, keep, drop = FALSE]

  # eQTL mapping is underpowered below ~15 individuals; drop the celltype
  # entirely rather than let it crash downstream residual/PC steps.
  min_individuals <- 15
  if (ncol(combined) < min_individuals) {
    message(sprintf("[COMBINE] %s: dropped entirely - only %d individuals post pseudobulk (< %d minimum)",
                    ct, ncol(combined), min_individuals))
    next
  }

  out <- data.frame(geneid = rownames(combined), combined, check.names = FALSE)
  fwrite(out, file.path(outdir, paste0(ct, "_pseudobulk.csv")))
  all_genes <- union(all_genes, rownames(combined))
  message(sprintf("[COMBINE] %s: %d genes x %d individuals; dropped %d below min_cells=%s",
                  ct, nrow(combined), ncol(combined), length(dropped), min_cells))
}

if (exists("get_gene_locations")) {
  ref <- data.frame(row.names = all_genes)
  gl <- get_gene_locations(ref)
  fwrite(gl, file.path(outdir, "gene_locations.csv"))
  message(sprintf("[COMBINE] gene_locations: %d genes mapped", nrow(gl)))
}
