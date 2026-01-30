process select_pcs_coarse {
    tag "Coarse PC selection for ${celltype}"
    label "process_low_memory"
    publishDir "${params.outdir}/optimization/", mode: 'copy'

    input:
    tuple val(celltype), path(egenes_files)

    output:
    tuple val(celltype), path("pc_values_fine.txt"), emit: fine_pc_values
    tuple val(celltype), path("*_coarse_summary.csv"), emit: coarse_summary

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    library(dplyr)

    # Combine all egenes files into a single data frame
    file_list <- unlist(strsplit("${egenes_files}", " "))
    results <- rbindlist(lapply(file_list, fread))

    results\$n_pcs <- as.integer(results\$n_pcs)
    results\$n_assoc <- as.numeric(results\$n_assoc)
    results <- results[order(results\$n_pcs)]

    early_tol <- ${params.pc_early_stop_tol}
    patience <- ${params.pc_early_stop_patience}

    stop_idx <- nrow(results)
    if (nrow(results) > 1) {
        gains <- c(Inf, (results\$n_assoc[-1] - results\$n_assoc[-nrow(results)]) / pmax(results\$n_assoc[-nrow(results)], 1))
        consec <- 0
        for (i in 2:nrow(results)) {
            if (is.finite(gains[i]) && gains[i] < early_tol) {
                consec <- consec + 1
                if (consec >= patience) {
                    stop_idx <- i
                    break
                }
            } else {
                consec <- 0
            }
        }
    }

    max_eval_pcs <- results\$n_pcs[stop_idx]
    results_eval <- results[n_pcs <= max_eval_pcs]
    center <- results_eval[which.max(n_assoc), ]\$n_pcs

    fine_window <- ${params.pc_fine_window}
    fine_step <- ${params.pc_fine_step}

    min_pcs <- min(results\$n_pcs)
    max_pcs <- max(results\$n_pcs)

    fine_min <- max(min_pcs, center - fine_window)
    fine_max <- min(max_pcs, center + fine_window)

    fine_vals <- seq(fine_min, fine_max, by = fine_step)
    fine_vals <- sort(unique(c(fine_vals, center, fine_min, fine_max)))

    if (length(fine_vals) == 0) {
        fine_vals <- center
    }

    write.table(fine_vals, "pc_values_fine.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

    writeLines(paste(
        "Cell type:", "${celltype}",
        "\nEarly-stop tol:", early_tol,
        "\nEarly-stop patience:", patience,
        "\nStop at n_pcs:", max_eval_pcs,
        "\nCoarse best n_pcs:", center
    ), "pc_info_${celltype}_coarse.txt")

    # Output detailed summary CSV for visualization
    results_summary <- results
    results_summary\$gain <- c(NA, (results\$n_assoc[-1] - results\$n_assoc[-nrow(results)]) / pmax(results\$n_assoc[-nrow(results)], 1))
    results_summary\$gain_pct <- results_summary\$gain * 100
    results_summary\$below_threshold <- results_summary\$gain_pct < (early_tol * 100)
    results_summary\$is_ceiling <- results_summary\$n_pcs == max_eval_pcs
    results_summary\$is_best <- results_summary\$n_pcs == center
    write.csv(results_summary, file = paste0("${celltype}_coarse_summary.csv"), row.names = FALSE, quote = TRUE)
    """
}
