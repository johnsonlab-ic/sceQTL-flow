process merge_seurat_objects {
    tag "merge_${seurat_files.size()}_objects"
    label "process_high_memory"
    publishDir "${params.outdir}/merged_objects/", mode: "copy"

    input: 
    path seurat_files
    path source_R

    output:
    path "merged_seurat.rds", emit: merged_seurat
    
    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    library(BPCells)
    library(dplyr)

    # Helper to print memory usage (process RSS and system memory)
    print_mem <- function(stage) {
        pid <- Sys.getpid()
        rss_kb <- suppressWarnings(as.numeric(system2("ps", args = c("-o", "rss=", "-p", pid), stdout = TRUE)))
        if (!is.na(rss_kb)) {
            cat(sprintf("[MEM] %s | RSS: %.2f GB\n", stage, rss_kb / 1024 / 1024))
        } else {
            cat(sprintf("[MEM] %s | RSS: unknown\n", stage))
        }
        free_out <- tryCatch(system2("free", args = c("-m"), stdout = TRUE), error = function(e) NULL)
        if (!is.null(free_out)) {
            cat("[MEM] free -m output:\n")
            cat(paste(free_out, collapse = "\n"), "\n")
        }
        cat("[MEM] gc() summary:\n")
        print(gc())
    }

    # Helper to validate integer-valued counts
    is_integer_counts <- function(mat, tol = 1e-8) {
        if (inherits(mat, "dgCMatrix")) {
            x_vals <- mat@x
            return(all(abs(x_vals - round(x_vals)) < tol))
        } else if (is.matrix(mat)) {
            return(all(abs(mat - round(mat)) < tol))
        } else {
            # Try to coerce to matrix for a strict check
            m <- tryCatch(as.matrix(mat), error = function(e) NULL)
            if (is.null(m)) return(FALSE)
            return(all(abs(m - round(m)) < tol))
        }
    }

    source("$source_R")

    # Get all RDS files
    seurat_files <- list.files(pattern = "\\\\.rds\$", full.names = TRUE)
    print_mem("start: after listing RDS files")
    
    cat("Found", length(seurat_files), "Seurat objects to merge\\n")
    cat("Files:", paste(seurat_files, collapse = ", "), "\\n")
    
    # =====================================================================
    # PASS 1: Quick scan to get gene and celltype intersection across all objects
    # =====================================================================
    cat("\\n=== PASS 1: Scanning for common genes and celltypes ===\\n")
    all_genes <- NULL
    all_celltypes <- NULL
    
    for (i in seq_along(seurat_files)) {
        cat(sprintf("Scanning object %d: %s\\n", i, basename(seurat_files[i])))
        obj <- readRDS(seurat_files[i])
        
        # Validate metadata columns and assay exist
        if (!("${params.celltype_column}" %in% colnames(obj@meta.data))) {
            stop(sprintf("Object %d missing celltype column: ${params.celltype_column}", i))
        }
        if (!("${params.individual_column}" %in% colnames(obj@meta.data))) {
            stop(sprintf("Object %d missing individual column: ${params.individual_column}", i))
        }
        if (!("${params.counts_assay}" %in% names(obj@assays))) {
            stop(sprintf("Object %d does not have assay: ${params.counts_assay}", i))
        }
        
        # Get gene names
        genes <- rownames(obj[["${params.counts_assay}"]]\$counts)
        cat(sprintf("  Object %d: %d cells, %d genes\\n", i, ncol(obj), length(genes)))
        
        # Get celltypes
        celltypes <- unique(obj@meta.data[["${params.celltype_column}"]])
        cat(sprintf("  Object %d celltypes: %s\\n", i, paste(celltypes, collapse = ", ")))
        
        # Intersect genes
        if (is.null(all_genes)) {
            all_genes <- genes
        } else {
            all_genes <- intersect(all_genes, genes)
        }
        cat(sprintf("  Common genes after object %d: %d\\n", i, length(all_genes)))
        
        # Intersect celltypes
        if (is.null(all_celltypes)) {
            all_celltypes <- celltypes
        } else {
            all_celltypes <- intersect(all_celltypes, celltypes)
        }
        cat(sprintf("  Common celltypes after object %d: %s\\n", i, paste(all_celltypes, collapse = ", ")))
        
        rm(obj)
        gc()
    }
    
    if (length(all_genes) == 0) {
        stop("No common genes found across all objects!")
    }
    if (length(all_celltypes) == 0) {
        stop("No common celltypes found across all objects!")
    }
    
    cat(sprintf("\\n=== Final common gene set: %d genes ===\\n", length(all_genes)))
    cat(sprintf("=== Final common celltype set: %s ===\\n\\n", paste(all_celltypes, collapse = ", ")))
    
    # =====================================================================
    # PASS 2: Extract counts matrices (subset to common genes and celltypes), write to disk
    # =====================================================================
    cat("=== PASS 2: Writing BPCells matrices with common genes and celltypes ===\\n")
    metadata_list <- list()
    bpcells_dirs <- character(length(seurat_files))
    
    for (i in seq_along(seurat_files)) {
        cat(sprintf("\\n=== Processing object %d: %s ===\\n", i, basename(seurat_files[i])))
        print_mem(sprintf("before readRDS %d", i))
        obj <- readRDS(seurat_files[i])
        print_mem(sprintf("after readRDS %d", i))
        
        cat(sprintf("  Original: %d cells, %d genes\\n", ncol(obj), nrow(obj)))
        
        # Subset to common celltypes first
        cells_to_keep <- obj@meta.data[["${params.celltype_column}"]] %in% all_celltypes
        obj <- obj[, cells_to_keep]
        cat(sprintf("  After subsetting to common celltypes: %d cells\\n", ncol(obj)))
        
        # Extract counts matrix and subset to common genes
        counts_matrix <- obj[["${params.counts_assay}"]]\$counts
        cat(sprintf("  Counts matrix dimensions: %d genes x %d cells\\n", nrow(counts_matrix), ncol(counts_matrix)))
        
        # Subset to common genes (preserves order from all_genes)
        counts_matrix <- counts_matrix[all_genes, ]
        cat(sprintf("  After subsetting to common genes: %d genes x %d cells\\n", nrow(counts_matrix), ncol(counts_matrix)))
        
        # Validate integer counts
        if (!is_integer_counts(counts_matrix)) {
            print_mem(sprintf("validation failure %d", i))
            stop(paste0(
                "Counts matrix contains non-integer values. Raw integer counts required in 'counts' slot of assay '",
                "${params.counts_assay}", "'. Offending file: ", basename(seurat_files[i]),
                ". Ensure you provided unnormalized counts (not data/scaled values)."
            ))
        }
        cat("  ✓ Counts validated as integer-valued\\n")
        
        # Write counts to BPCells on-disk format
        counts_dir <- paste0("bpcells_obj", i, "_counts")
        cat(sprintf("  Writing BPCells matrix to: %s\\n", counts_dir))
        print_mem(sprintf("before write_matrix_dir %d", i))
        write_matrix_dir(mat = counts_matrix, dir = counts_dir)
        print_mem(sprintf("after write_matrix_dir %d", i))
        bpcells_dirs[i] <- counts_dir
        
        # Extract and standardize metadata - keep only required columns plus any extras
        meta <- obj@meta.data
        meta\$orig.file <- basename(seurat_files[i])
        
        # Ensure consistent column order and presence
        required_cols <- c("${params.celltype_column}", "${params.individual_column}", "orig.file")
        missing_cols <- setdiff(required_cols, colnames(meta))
        if (length(missing_cols) > 0) {
            stop(sprintf("Object %d missing required columns: %s", i, paste(missing_cols, collapse = ", ")))
        }
        
        metadata_list[[i]] <- meta
        cat(sprintf("  Extracted metadata: %d rows, %d columns\\n", nrow(meta), ncol(meta)))
        
        # Explicitly remove object to free memory
        rm(obj, counts_matrix, meta)
        gc()
        print_mem(sprintf("after removing object %d", i))
    }
    
    cat("\\n=== PASS 2 COMPLETE: All matrices written to disk ===\\n\\n")
    
    # =====================================================================
    # PASS 3: Open all BPCells matrices and create Seurat object from named list
    # =====================================================================
    cat("=== PASS 3: Opening BPCells matrices from disk ===\\n")
    data_list <- list()
    
    for (i in seq_along(bpcells_dirs)) {
        cat(sprintf("Opening %s\\n", bpcells_dirs[i]))
        print_mem(sprintf("before open_matrix_dir %d", i))
        data_list[[i]] <- open_matrix_dir(dir = bpcells_dirs[i])
        print_mem(sprintf("after open_matrix_dir %d", i))
    }
    
    # Name the list
    names(data_list) <- paste0("obj", seq_along(data_list))
    cat(sprintf("Opened %d on-disk BPCells matrices: %s\\n", length(data_list), paste(names(data_list), collapse = ", ")))
    
    # Verify all matrices have same gene set
    gene_counts <- sapply(data_list, nrow)
    cat(sprintf("Gene counts per matrix: %s\\n", paste(gene_counts, collapse = ", ")))
    if (length(unique(gene_counts)) > 1) {
        warning(sprintf("Matrices have different numbers of genes! %s", paste(gene_counts, collapse = ", ")))
    }
    
    # Combine metadata - find common columns across all metadata
    cat("\\n=== Combining metadata ===\\n")
    print_mem("before metadata alignment")
    
    # Get intersection of columns across all metadata data.frames
    common_cols <- Reduce(intersect, lapply(metadata_list, colnames))
    cat(sprintf("Common metadata columns: %s\\n", paste(common_cols, collapse = ", ")))
    
    # Subset each metadata to common columns in same order
    metadata_list <- lapply(metadata_list, function(x) x[, common_cols, drop = FALSE])
    
    # Combine metadata with rbind
    metadata <- Reduce(rbind, metadata_list)
    print_mem("after rbind metadata")
    cat(sprintf("Combined metadata: %d rows, %d columns\\n", 
                nrow(metadata), ncol(metadata)))
    
    # Create Seurat object from named list of matrices (as per standard BPCells merge approach)
    cat("\\n=== Creating merged Seurat object ===\\n")
    print_mem("before CreateSeuratObject")
    merged_seurat <- CreateSeuratObject(
        counts = data_list,
        meta.data = metadata,
        assay = "${params.counts_assay}"
    )
    print_mem("after CreateSeuratObject")
    
    cat("\\n=== Merge Summary ===\\n")
    cat(sprintf("Merged object: %d cells, %d features\\n", ncol(merged_seurat), nrow(merged_seurat)))
    cat(sprintf("Cell types: %s\\n", paste(unique(merged_seurat@meta.data[["${params.celltype_column}"]]), collapse = ", ")))
    cat(sprintf("Individuals: %s\\n", paste(unique(merged_seurat@meta.data[["${params.individual_column}"]]), collapse = ", ")))
    
    # Save merged object (BPCells-backed assays will be on-disk)
    cat("Saving merged object...\\n")
    print_mem("before saveRDS")
    saveRDS(merged_seurat, "merged_seurat.rds")
    print_mem("after saveRDS")
    cat("Merge complete!\\n")
    """
}
