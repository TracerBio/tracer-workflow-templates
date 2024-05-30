#!/bin/bash
 
# Specify input files and parameters
FASTQ_FILES="./rnaseq/data/human"
GENOME_DIR="./rnaseq/data/human/hg19"
OUTPUT_DIR="./rnaseq/data/human/STAR"
 
# Start a new Tracer run
tracer log "Started Tracer run"
 
# Record the execution of the STAR tool
tracer tool STAR v2.7.10a

# Generate genome index
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir ./rnaseq/data/human/hg19 \
     --genomeSAindexNbases 10 \
     --genomeFastaFiles ./rnaseq//data/human/human.fa \
     --sjdbGTFfile ./rnaseq/data/human/hg19.refGene.gtf \
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
 
# Mark the end of a pipeline run
tracer end
 