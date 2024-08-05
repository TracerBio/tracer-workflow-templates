#!/bin/bash

# Directory and file paths
cd /workspace/tracer-workflow-templates/data

# Run fastqc 10 times and measure execution time
for i in {1..10}; do
    echo "Running FASTQC iteration $i..."
    
    # Record start time for the current run
    start_fastqc=$(date +%s%3N)
    
    # Execute fastqc
    fastqc s1_1.fq s1_2.fq -o .
    
    # Record end time for the current run
    end_fastqc=$(date +%s%3N)
    
    # Calculate the duration in milliseconds
    fastqc_duration=$((end_fastqc - start_fastqc))
    
    # Print the duration for the current run
    echo "Execution duration of FASTQC iteration $i: $fastqc_duration milliseconds"
done
