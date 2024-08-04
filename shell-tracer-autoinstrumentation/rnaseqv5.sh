#!/bin/bash

# FASTQC - QC layer 1
cd /workspace/tracer-workflow-templates/data

# Creating a human genome index using STAR - STAR-INDEX
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./human_star --genomeSAindexNbases 10 --genomeFastaFiles human.fa --sjdbGTFfile hg19_anno.gtf --sjdbOverhang 99;

# Align RNA-Seq reads of control to the genome using STAR - STAR-CONTROL-MAP
STAR --runThreadN 4 --genomeDir ./human_star --readFilesIn control1_1.fq control1_2.fq --outFileNamePrefix control1_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts;

# Align RNA-Seq reads of test experiment to the genome using STAR - STAR-TEST-MAP
STAR --runThreadN 4 --genomeDir ./human_star --readFilesIn test1_1.fq test1_2.fq --outFileNamePrefix test1_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts;

# Sort and Index the output sam files - SAMTOOLS-CONTROL
mv control1_starAligned.sortedByCoord.out.bam control.sorted.bam;
samtools index control.sorted.bam;

# Sort and Index the output sam files - SAMTOOLS-TEST
mv test1_starAligned.sortedByCoord.out.bam test.sorted.bam;
samtools index test.sorted.bam;

# Predict potential transcripts in control and test BAM files using StringTie - STRINGTIE
stringtie -o control.gtf -G hg19_anno.gtf control.sorted.bam;
stringtie -o test.gtf -G hg19_anno.gtf test.sorted.bam;

# Calculate transcript counts from BAM files using featureCounts - FEATURECOUNTS
featureCounts -p --countReadPairs -B -C -T 4 -a hg19_anno.gtf -o control_counts.txt control.sorted.bam;
featureCounts -p --countReadPairs -B -C -T 4 -a hg19_anno.gtf -o test_counts.txt test.sorted.bam;

# Collate logs and perform post-mapping QC with MULTIQC - MULTIQC
multiqc .;

# Summary of BAM files using deeptools - BAMSUMMARY
multiBamSummary bins --bamfiles control.sorted.bam test.sorted.bam -o rnaseq.npz; 

# PCA analysis for RNASEQ experiments - PCA
plotPCA -in rnaseq.npz -o PCA_rnaseq.png;

# RNASEQ data comparison via Fingerprint plots - FINGERPRINT
plotFingerprint -b control.sorted.bam test.sorted.bam --labels Control Test --plotFile fingerprint_rnaseq.png;

# Obtain bigwig files for visualization of data - BAMCOMPARE
bamCompare -b1 test.sorted.bam -b2 control.sorted.bam -o differential.bw -of bigwig

# Run the R script for edgeR analysis and heatmap creation
Rscript DESeq2.r

echo "Analysis completed. edgeR results and heatmap saved."