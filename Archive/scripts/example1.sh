#!/bin/bash
 
# Specify input files and parameters
FASTQ_FILES="/workspace/tracer-workflow-templates/data/"
GENOME_DIR="/workspace/tracer-workflow-templates/data/hg19"
OUTPUT_DIR="/workspace/tracer-workflow-templates/data/STAR"
 
# Generate genome index
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /workspace/tracer-workflow-templates/data/hg19 \
     --genomeSAindexNbases 10 \
     --genomeFastaFiles /workspace/tracer-workflow-templates/data/human.fa \
     --sjdbGTFfile /workspace/tracer-workflow-templates/data/hg19.refGene.gtf \
     --sjdbOverhang 99

 
# Align RNA-Seq reads to the genome using STAR
STAR --runThreadN 8 \
     --genomeDir "${GENOME_DIR}" \
     --readFilesIn "${FASTQ_FILES}/aml_1.fq" "${FASTQ_FILES}/aml_2.fq" \
     --outFileNamePrefix "${OUTPUT_DIR}/aml_" \
     --outSAMtype BAM SortedByCoordinate