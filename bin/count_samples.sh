#!/bin/bash

# Count the number of columns (samples) in a CSV file
# Usage: count_samples.sh <csv_file>

# Parse CSV header and count fields (excluding the first gene ID column)
if [[ "$1" == *.csv ]]; then
    # For CSV files
    awk -F, 'NR==1 {print NF-1; exit}' "$1"
elif [[ "$1" == *.txt ]]; then
    # For tab-delimited files
    awk -F '\t' 'NR==1 {print NF-1; exit}' "$1"
else
    # Default case
    head -1 "$1" | tr -dc ',' | awk '{print length($0) + 1}'
fi
