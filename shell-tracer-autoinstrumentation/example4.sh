#!/bin/bash
 
# Specify input files and parameters
FASTQ_FILES="./data/"
GENOME_DIR="./data/hg19"
OUTPUT_DIR="./data/STAR"
 
# Generate genome index
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir ./data/hg19 \
     --genomeSAindexNbases 10 \
     --genomeFastaFiles ./data/human.fa \
     --sjdbGTFfile ./data/hg19.refGene.gtf \
     --sjdbOverhang 99

 
# Align RNA-Seq reads to the genome using STAR
STAR --runThreadN 8 \
     --genomeDir "${GENOME_DIR}" \
     --readFilesIn "${FASTQ_FILES}/aml_1.fq" "${FASTQ_FILES}/aml_2.fq" \
     --outFileNamePrefix "${OUTPUT_DIR}/aml_" \
     --outSAMtype BAM SortedByCoordinate