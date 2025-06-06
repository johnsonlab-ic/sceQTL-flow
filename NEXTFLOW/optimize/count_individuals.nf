process count_individuals {
    tag "${residuals_file}"
    label "process_low"
    
    input:
    path residuals_file

    output:
    tuple path(residuals_file), path('pc_values.txt'), emit: residuals_with_pcs

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)

    # Read the residuals file
    exp_mat = fread("$residuals_file")
    
    # Count the number of columns (excluding the geneid column)
    n_samples = ncol(exp_mat) - 1  # Subtract 1 for geneid column
    
    # Calculate the maximum number of PCs as 50% of sample count, capped at 20
    max_pcs = min(floor(n_samples * 0.5), 100)
    
    # Generate PC values in steps of 2 (adjust step size as needed)
    pc_values = seq(2, max_pcs, by = 2)
    
    # If very few samples, ensure at least one PC value is tested
    if(length(pc_values) == 0) {
        pc_values = min(1, max_pcs)
    }
    
    # Write the PC values to a file that Nextflow can read
    write.table(pc_values, "pc_values.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    """
}

