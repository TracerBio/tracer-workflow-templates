#!/bin/bash
 
# Specify input files and parameters
FASTQ_FILES="/workspaces/tracer-workflow-templates/data/"
GENOME_DIR="/workspaces/tracer-workflow-templates/data/hg19"
OUTPUT_DIR="/workspaces/tracer-workflow-templates/data/STAR"
 
# Start a new Tracer run
tracer log "Started Tracer run"
 
# Record the execution of the STAR tool
tracer tool STAR v2.7.10a

# Generate genome index
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /workspaces/tracer-workflow-templates/data/hg19 \
     --genomeSAindexNbases 10 \
     --genomeFastaFiles /workspaces/tracer-workflow-templates/data/human.fa \
     --sjdbGTFfile /workspaces/tracer-workflow-templates/data/hg19.refGene.gtf \
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
 