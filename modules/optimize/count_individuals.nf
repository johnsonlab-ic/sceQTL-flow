process count_individuals {
    tag "${residuals_file}"
    label "process_low"
    
    input:
    path residuals_file

    output:
    tuple path(residuals_file), path('pc_values_coarse.txt'), emit: residuals_with_pcs

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)

    # Read the residuals file
    exp_mat = fread("$residuals_file")
    
    # Count the number of columns (excluding the geneid column)
    n_samples = ncol(exp_mat) - 1  # Subtract 1 for geneid column
    
    # Calculate the maximum number of PCs as a fraction of sample count, capped
    max_pcs = min(floor(n_samples * ${params.pc_max_fraction}), ${params.pc_max_cap})
    min_pcs = ${params.pc_min}
    coarse_step = ${params.pc_coarse_step}
    
    # Generate coarse PC values
    if (max_pcs >= min_pcs) {
        pc_values = seq(min_pcs, max_pcs, by = coarse_step)
        pc_values = sort(unique(c(pc_values, max_pcs)))
    } else {
        pc_values = max_pcs
    }
    
    # Write the coarse PC values to a file that Nextflow can read
    write.table(pc_values, "pc_values_coarse.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    """
}
