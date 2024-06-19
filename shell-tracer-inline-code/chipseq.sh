#!/bin/bash

# Start a new Tracer run
tracer start;
tracer log "Started Tracer run";
 
# FASTQC - QC layer 1
cd ~/workspaces/tracer-workflow-templates/data;

tracer tool fastqc 0.12.1;
fastqc s1_1.fq s1_2.fq -o .;

# Record the execution of the FASTQC tool
tracer tool star-index 1.7.0;

# Creating a human genome index using STAR
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./human_genome --genomeSAindexNbases 10 --genomeFastaFiles human.fa --sjdbGTFfile hg19.refGene.gtf --sjdbOverhang 99;


# Align RNA-Seq reads to the genome using STAR
tracer tool star-align 1.7.0;
STAR --runThreadN 4 --genomeDir ./human_genome --readFilesIn s1_1.fq s1_2.fq --outFileNamePrefix P1s1 --outSAMtype BAM SortedByCoordinate;

# Sort and Index the output bam file
tracer tool samtools 1.17;
samtools sort *.bam -@ 4 -o P1s1.sorted.bam;
samtools index P1s1.sorted.bam;

# Call Peaks for BAM file
tracer tool macs 3.0.1;
macs3 callpeak -t P1s1.sorted.bam -f BAMPE -p 0.05 --outdir ./Peaks;

# Plot Coverage of Peaks
tracer tool deeptools 3.5.5;
plotCoverage -b P1s1.sorted.bam -o coverage.pdf;
tracer end