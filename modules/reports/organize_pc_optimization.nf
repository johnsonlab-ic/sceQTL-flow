process organize_pc_optimization {
    label "process_single"
    publishDir "${params.outdir}/optimization", mode: 'copy', pattern: "*.csv"

    input:
    path coarse_summaries
    path fine_summaries

    output:
    path "*.csv", emit: summary_csvs, optional: true

    script:
    """
    # Organize coarse and fine summary CSVs into optimization directory
    
    # Handle coarse summaries (could be multiple files)
    for file in ${coarse_summaries}; do
        if [ -f "\$file" ]; then
            cp "\$file" .
            echo "Copied coarse summary: \$file"
        fi
    done
    
    # Handle fine summaries (could be multiple files)
    for file in ${fine_summaries}; do
        if [ -f "\$file" ]; then
            cp "\$file" .
            echo "Copied fine summary: \$file"
        fi
    done
    
    # List what was copied for debugging
    echo "Final contents of optimization directory:"
    ls -la *.csv 2>/dev/null || echo "No CSV files found"
    """
}
