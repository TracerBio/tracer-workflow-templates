#!/bin/bash
 
# Specify input files and parameters
FASTQ_FILES="$projectDir/data/human"
GENOME_DIR="$projectDir/data/human/hg19"
OUTPUT_DIR="$projectDir/data/human/STAR"
 
# Start a new Tracer run
tracer log "Started Tracer run"
 
# Record the execution of the STAR tool
tracer tool STAR v2.7.10a

# Generate genome index
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir $projectDir/data/human/hg19 \
     --genomeFastaFiles $projectDir/data/human/human.fa \
     --sjdbGTFfile $projectDir/data/human/hg19.refGene.gtf \
     --sjdbOverhang 99

# Log progress
tracer log "Genome index generated"
 
# Align RNA-Seq reads to the genome using STAR
STAR --runThreadN 8 \
     --genomeDir "${GENOME_DIR}" \
     --readFilesIn "${FASTQ_FILES}/aml_1.fq" "${FASTQ_FILES}/aml_2.fq" \
     --outFileNamePrefix "${OUTPUT_DIR}/aml_" \
     --outSAMtype BAM SortedByCoordinate
 
# Log progress
tracer log "STAR alignment completed for sample"
 
# Record another tool for further downstream analysis (e.g., featureCounts)
tracer tool featureCounts v2.0.1
 
# Example command for counting reads per feature using featureCounts
featureCounts -T 8 \
              -a $projectDir/data/human/hg19.refGene.gtf \
              -o "${OUTPUT_DIR}/counts.txt" \
              "${OUTPUT_DIR}/sample_Aligned.sortedByCoord.out.bam"
 
# Log progress
tracer log "Feature counting completed"
 
# Mark the end of a pipeline run
tracer end
 