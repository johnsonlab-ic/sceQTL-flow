process subset_samples {
    tag "Subsetting by ${subset_column}"
    label "process_single"
    publishDir "${params.outdir}/QC/", mode: 'copy'

    input:
    path cov_file
    val subset_column
    val subset_values

    output:
    path "*_subset_covariates.csv", emit: filtered_cov
    path "subset_report.txt", emit: subset_report

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(dplyr)
    
    cov <- fread("$cov_file")
    n_initial <- nrow(cov)
    
    # Parse values to keep
    values_to_keep <- trimws(unlist(strsplit("${subset_values}", ",")))
    
    # Filter rows
    cov_filtered <- cov %>% filter(.data[["${subset_column}"]] %in% values_to_keep)
    n_final <- nrow(cov_filtered)
    
    # Save filtered covariate file
    fwrite(cov_filtered, "subset_subset_covariates.csv")
    
    # Generate report
    report <- c(
        "===========================================",
        "SAMPLE SUBSETTING REPORT",
        "===========================================",
        paste0("Column: ${subset_column}"),
        paste0("Values kept: ", paste(values_to_keep, collapse=", ")),
        paste0("Initial samples: ", n_initial),
        paste0("Filtered samples: ", n_final),
        paste0("Samples removed: ", n_initial - n_final),
        paste0("Retention: ", round(100 * n_final / n_initial, 1), "%"),
        "==========================================="
    )
    writeLines(report, "subset_report.txt")
    cat(paste(report, collapse="\\n"))
    """
}
