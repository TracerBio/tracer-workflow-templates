#!/bin/bash

# FASTQC - QC layer 1
cd /workflows/tracer-workflow-templates/data

# Record the execution of the FASTQC tool
fastqc control1_1.fq control1_2.fq -o .;
fastqc test1_1.fq test1_2.fq -o .;

# Creating a human genome index using bowtie2
bowtie2-build human.fa human_index;

# Align RNA-Seq reads of control to the genome using bowtie2
bowtie2 --local -x human_index -1 control_1.fq -2 control_2.fq -S control.sam;

# Align RNA-Seq reads of test experiment to the genome using bowtie2
bowtie2 --local -x human_index -1 test_1.fq -2 test_2.fq -S test.sam;

# Sort and Index the output sam files
samtools view -@ 4 -b control.sam > control.bam 
samtools sort control.bam -@ 4 -o control.sorted.bam;
samtools index control.sorted.bam;
samtools view -@ 4 -b test.sam > test.bam 
samtools sort test.bam -@ 4 -o test.sorted.bam;
samtools index test.sorted.bam;

# Predict potential transcripts in control and test BAM files using StringTie
stringtie -o control.gtf -G hg19.refGene.gtf control.sorted.bam;
stringtie -o test.gtf -G hg19.refGene.gtf test.sorted.bam;

# Calculate transcript counts from BAM files using featureCounts
featureCounts -p --countReadPairs -B -C -T 4 -a hg19.refGene.gtf -o control_counts control.sorted.bam;
featureCounts -p --countReadPairs -B -C -T 4 -a hg19.refGene.gtf -o test_counts test.sorted.bam;

# Collate logs and perform post-mapping QC with MULTIQC
multiqc .;

# Summary of BAM files using deeptools
multiBamSummary bins --bamfiles control.sorted.bam test.sorted.bam -o rnaseq.npz; 

# PCA analysis for RNASEQ experiments
plotPCA -in rnaseq.npz -o PCA_rnaseq.png;

# RNASEQ data comparison via Fingerprint plots
plotFingerprint -b control.sorted.bam test.sorted.bam --labels Control Test --plotFile fingerprint_rnaseq.png;

# Obtain bigwig files for visualization of data
bamCompare -b1 test.sorted.bam -b2 control.sorted.bam -o differential.bw -of bigwig