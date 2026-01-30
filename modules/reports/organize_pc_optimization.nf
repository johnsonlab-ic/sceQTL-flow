process organize_pc_optimization {
    label "process_single"
    publishDir "${params.outdir}/optimization", mode: 'copy', pattern: "*.csv"

    input:
    path coarse_summaries
    path fine_summaries

    output:
    path "*.csv", emit: summary_csvs optional true

    script:
    """
    # Organize coarse and fine summary CSVs into optimization directory
    if [ -f "${coarse_summaries}" ]; then
        cp "${coarse_summaries}" . || true
    fi
    if [ -f "${fine_summaries}" ]; then
        cp "${fine_summaries}" . || true
    fi
    
    # List what was copied for debugging
    ls -la *.csv 2>/dev/null || echo "No CSV files to organize"
    """
}
