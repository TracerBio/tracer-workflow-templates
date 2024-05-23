#!/bin/bash
 
# Specify input files and parameters
FASTQ_FILES="/path/to/fastq_files"
GENOME_DIR="/path/to/genome_dir"
OUTPUT_DIR="/path/to/alignment_output"
 
# Start a new Tracer run
tracer start
 
# Record the execution of the STAR tool
tracer tool STAR v2.7.10a
 
# Generate Genome directory
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir hg19/ \
     --genomeFastaFiles hg19.fa \ 
     --sjdbGTFfile hg19.annotated.gtf \ 
     --sjdbOverhang 99

# Align RNA-Seq reads to the genome using STAR
STAR --runThreadN 8 \
     --genomeDir ${GENOME_DIR} \
     --readFilesIn ${FASTQ_FILES}/sample_1.fastq ${FASTQ_FILES}/sample_2.fastq \
     --outFileNamePrefix ${OUTPUT_DIR}/sample \
     --outSAMtype BAM
 
# Log progress
tracer log "STAR alignment completed for sample"
 
# Record another tool for further downstream analysis (e.g., featureCounts)
tracer tool featureCounts v2.0.1
 
# Example command for counting reads per feature using featureCounts
featureCounts -T 8 \
              -a /path/to/annotation.gtf \
              -o "${OUTPUT_DIR}/counts.txt" \
              "${OUTPUT_DIR}/sample_Aligned.sortedByCoord.out.bam"
 
# Log progress
tracer log "Feature counting completed"
 
# Mark the end of a pipeline run
tracer end
 