process preflight_check {
    tag "Individual Sample Check"
    label "process_single"
    publishDir "${params.outdir}/QC/", mode: 'copy'

    input:
    path genotype_mat
    path cov_file
    path pseudobulk_files
    val has_cov_file

    output:
    path "preflight_check_report.txt", emit: preflight_report

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(dplyr)

    report <- list()
    report_text <- c()
    
    # =============================================
    # 1. Check Genotype File
    # =============================================
    geno <- fread("${genotype_mat}") %>% tibble::column_to_rownames(var="snp")
    n_geno_samples <- ncol(geno)
    geno_samples <- colnames(geno)
    
    report_text <- c(report_text, 
                     "===========================================",
                     "PREFLIGHT CHECK: Individual Sample Inventory",
                     "===========================================",
                     "",
                     sprintf("[GENOTYPE] %d individuals", n_geno_samples),
                     sprintf("  Sample IDs: %s", paste(head(geno_samples, 3), collapse=", "))
                     )
    
    # =============================================
    # 2. Check Covariate File (if provided)
    # =============================================
    has_cov <- as.logical("${has_cov_file}")
    if (has_cov) {
        cov <- fread("${cov_file}")
        n_cov_samples <- nrow(cov)
        cov_samples <- cov[[1]]  # Assume first column is Individual_ID
        
        report_text <- c(report_text,
                         "",
                         sprintf("[COVARIATES] %d individuals", n_cov_samples),
                         sprintf("  Sample IDs: %s", paste(head(cov_samples, 3), collapse=", "))
                         )
        
        # Compare genotype vs covariates
        geno_in_cov <- sum(geno_samples %in% cov_samples)
        cov_in_geno <- sum(cov_samples %in% geno_samples)
        
        report_text <- c(report_text,
                         sprintf("  [MATCH] %d/%d genotype samples in covariates (%.1f%%)", 
                                 geno_in_cov, n_geno_samples, 100*geno_in_cov/n_geno_samples),
                         sprintf("  [MATCH] %d/%d covariate samples in genotype (%.1f%%)", 
                                 cov_in_geno, n_cov_samples, 100*cov_in_geno/n_cov_samples)
                         )
        
        if (geno_in_cov < 0.9 * n_geno_samples) {
            report_text <- c(report_text,
                             "  ⚠️  WARNING: <90% genotype samples matched to covariates!")
        }
    }
    
    # =============================================
    # 3. Check Pseudobulk Files
    # =============================================
    pb_files <- list.files(pattern = "_pseudobulk\\\\.csv\$")
    
    if (length(pb_files) > 0) {
        report_text <- c(report_text, "")
        report_text <- c(report_text, "[PSEUDOBULK MATRICES]")
        
        pb_sample_counts <- data.frame(celltype = character(), n_samples = integer(), stringsAsFactors = FALSE)
        
        for (file in pb_files) {
            pb <- fread(file) %>% tibble::column_to_rownames(var="geneid")
            n_pb_samples <- ncol(pb)
            celltype <- gsub("_pseudobulk\\\\.csv", "", file)
            pb_sample_counts <- rbind(pb_sample_counts, 
                                      data.frame(celltype = celltype, n_samples = n_pb_samples))
            
            pb_samples <- colnames(pb)
            report_text <- c(report_text,
                             sprintf("  %s: %d individuals", celltype, n_pb_samples),
                             sprintf("    Sample IDs: %s", paste(head(pb_samples, 3), collapse=", "))
                             )
            
            # Check match to genotype
            pb_in_geno <- sum(pb_samples %in% geno_samples)
            report_text <- c(report_text,
                             sprintf("    [MATCH] %d/%d pseudobulk samples in genotype (%.1f%%)", 
                                     pb_in_geno, n_pb_samples, 100*pb_in_geno/n_pb_samples)
                             )
            
            if (pb_in_geno < 0.9 * n_pb_samples) {
                report_text <- c(report_text,
                                 sprintf("    ⚠️  WARNING: <%d%% pseudobulk samples in genotype!", 90)
                                 )
            }
        }
        
        # Summary stats
        min_pb <- min(pb_sample_counts\$n_samples)
        max_pb <- max(pb_sample_counts\$n_samples)
        report_text <- c(report_text,
                         "",
                         sprintf("  SUMMARY: %d cell types, %d-%d samples per cell type",
                                 nrow(pb_sample_counts), min_pb, max_pb)
                         )
        
        if (max_pb - min_pb > 0.2 * max_pb) {
            report_text <- c(report_text,
                             sprintf("  ⚠️  WARNING: Large variation in sample counts across cell types (%.0f%% difference)",
                                     100 * (max_pb - min_pb) / max_pb)
                             )
        }
    }
    
    report_text <- c(report_text, "", "===========================================")
    
    # Write report
    writeLines(report_text, "preflight_check_report.txt")
    
    # Also print to console
    cat(paste(report_text, collapse = "\\n"))
    """
}
