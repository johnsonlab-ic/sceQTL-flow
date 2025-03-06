process select_pcs {

    label "process_low_memory"
    publishDir "${params.outdir}/eQTL_outputs/", mode: 'copy'

    input:
    tuple val(celltype), path(egenes_files)

    output:
    path "optimal_pcs_${celltype}.txt", emit: optimal_pcs
    path "*png"

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(ggplot2)

    # Combine all files into a single data frame
    file_list <- unlist(strsplit("${egenes_files}", " "))
    results <- rbindlist(lapply(file_list, fread))


    #create plot
    plot <- ggplot(results, aes(x=n_pcs, y=n_egenes)) +
        geom_point() +
        geom_line() +
        ggtitle("${celltype}") +
        xlab("Number of Principal Components") +
        ylab("Number of eGenes") +
        theme_minimal()

    ggsave(filename=paste0("egenes_vs_pcs_${celltype}.png"), plot=plot)

    # Find the optimal n_pcs for the cell type
    optimal_pcs <- results[which.max(n_egenes), ]

    # Write the optimal n_pcs to a file
    fwrite(optimal_pcs, "optimal_pcs_${celltype}.txt")
    """

    
}
