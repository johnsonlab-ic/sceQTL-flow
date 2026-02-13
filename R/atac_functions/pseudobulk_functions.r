# Import the pipe operator
library(dplyr)
library(Matrix)

# Function to sanitize cell type names for use in file paths
sanitize_celltype_name <- function(celltype_name) {
  sanitized <- gsub("[/:*?\"<>|\\\\]", "_", celltype_name)
  sanitized <- gsub("_+", "_", sanitized)
  sanitized <- gsub("^_+|_+$", "", sanitized)
  return(sanitized)
}

# Function to get peak coordinates from peak names (chr:start-end format)
get_peak_locations <- function(peak_names) {
  parse_peak <- function(peak_name) {
    tokens <- unlist(strsplit(peak_name, "[:\\-]"))
    if (length(tokens) < 3) {
      return(c(NA, NA, NA))
    }
    return(tokens[1:3])
  }

  peak_mat <- do.call(rbind, lapply(peak_names, parse_peak))
  peak_df <- as.data.frame(peak_mat, stringsAsFactors = FALSE)
  colnames(peak_df) <- c("chr", "left", "right")
  peak_df$left <- suppressWarnings(as.numeric(peak_df$left))
  peak_df$right <- suppressWarnings(as.numeric(peak_df$right))
  peak_df$geneid <- peak_names

  peak_df <- peak_df %>% dplyr::select(geneid, chr, left, right)
  return(peak_df)
}

# Function to pseudobulk peaks (ATAC ChromatinAssay)
pseudobulk_peaks <- function(seuratlist, min.cells = 100, indiv_col = "sample_id", assay = "ATAC") {
  agg_peak_list <- lapply(seuratlist, function(x) {
    Seurat::DefaultAssay(x) <- assay
    metadata <- x[[]]

    assay_obj <- x[[assay]]
    counts <- tryCatch(
      Seurat::GetAssayData(x, assay = assay, slot = "counts"),
      error = function(e) {
        assay_obj@counts
      }
    )

    if (is.null(counts)) {
      stop("No counts matrix found in ", assay, " assay for this cell type")
    }

    unique_ids <- unique(metadata[[indiv_col]])
    indiv_table <- metadata %>% dplyr::count(get(indiv_col))
    indiv_table <- indiv_table %>% dplyr::filter(n > min.cells)
    colnames(indiv_table) <- c("unique_ids", "n")
    unique_ids <- indiv_table$unique_ids
    n_indiv_start <- length(unique_ids)

    if (nrow(indiv_table) == 0) {
      emptydf <- data.frame(matrix(ncol = 0, nrow = nrow(counts)))
      rownames(emptydf) <- rownames(counts)
      return(emptydf)
    } else {
      cells <- rownames(metadata[metadata[[indiv_col]] == unique_ids[1], ])
      counts_2 <- Matrix::rowSums(counts[, cells])
      finalmat <- data.frame(counts_2)
      samplenames <- unique_ids[1]

      if (length(unique_ids) > 1) {
        for (i in 2:length(unique_ids)) {
          cells <- rownames(metadata[metadata[[indiv_col]] == unique_ids[i], ])
          sample <- unique_ids[i]
          counts_2 <- Matrix::rowSums(counts[, cells])
          finalmat <- cbind(finalmat, counts_2)
          samplenames <- c(samplenames, sample)
        }
      }
      colnames(finalmat) <- samplenames
      indiv_remain <- n_indiv_start - length(samplenames)

      return(finalmat)
    }
  })
  return(agg_peak_list)
}