#!/bin/bash

# FASTQC - QC layer 1
cd /workspace/tracer-workflow-templates/data

# Creating a human genome index using bowtie2 - BOWTIE2-INDEX
bowtie2-build human.fa human_index;

# Align RNA-Seq reads of control to the genome using bowtie2 - BOWTIE2-CONTROL-MAP
bowtie2 --local -x human_index -1 control1_1.fq -2 control1_2.fq -S control.sam;

# Align RNA-Seq reads of test experiment to the genome using bowtie2 - BOWTIE2-TEST-MAP
bowtie2 --local -x human_index -1 test1_1.fq -2 test1_2.fq -S test.sam;

# Creating a human genome index using STAR - STAR-INDEX
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./human_star --genomeSAindexNbases 10 --genomeFastaFiles human.fa --sjdbGTFfile hg19_anno.gtf --sjdbOverhang 99;

# Align RNA-Seq reads of control to the genome using STAR - STAR-CONTROL-MAP
STAR --runThreadN 4 --genomeDir ./human_star --readFilesIn control1_1.fq control1_2.fq --outFileNamePrefix control1_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts;

# Align RNA-Seq reads of test experiment to the genome using STAR - STAR-TEST-MAP
STAR --runThreadN 4 --genomeDir ./human_star --readFilesIn test1_1.fq test1_2.fq --outFileNamePrefix test1_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts;

# Convert SAM files into transcriptome FASTA files - SAM2FASTA
samtools fasta control.sam > control.fa
samtools fasta test.sam > test.fa

# Creating a transcriptome index using kallisto - KALLISTO-INDEX
kallisto index -t 4 -i control_index control.fa;
kallisto index -t 4 -i test_index test.fa;

# Quantify RNA-Seq reads of control to the genome using kallisto - KALLISTO-CONTROL-QUANT
kallisto quant -t 4 -i control_index -o ./control_quant control1_1.fq control1_2.fq;

# Quantify RNA-Seq reads of test experiment to the genome using kallisto - KALLISTO-TEST-QUANT
kallisto quant -t 4 -i test_index -o ./test_quant test1_1.fq test1_2.fq;

# Creating a human genome index using salmon - SALMON-INDEX
salmon index --threads 4 -t human.fa -i salmon_index;

# Align RNA-Seq reads of control to the genome using salmon - SALMON-CONTROL-QUANT
salmon quant --threads 4 --libType=U -i salmon_index -1 control1_1.fq -2 control1_2.fq -o ./salmon_control;

# Align RNA-Seq reads of test experiment to the genome using salmon - SALMON-TEST-QUANT
salmon quant --threads 4 --libType=U -i salmon_index -1 test1_1.fq -2 test1_2.fq -o ./salmon_test;

# Creating a human genome index using hisat2 - HISAT2-INDEX
hisat2-build human.fa hisat2_index;

# Align RNA-Seq reads of control to the genome using hisat2 - HISAT2-CONTROL-ALIGN
hisat2 --fast -x hisat2_index -1 control1_1.fq -2 control1_2.fq -S control_hs2.sam;

# Align RNA-Seq reads of test experiment to the genome using hisat2 - HISAT2-TEST-ALIGN
hisat2 --fast -x hisat2_index -1 test1_1.fq -2 test1_2.fq -S test_hs2.sam;

# Creating a human genome index using bwa - BWA-INDEX
bwa index -p bwa_index human.fa;

# Align RNA-Seq reads of control to the genome using bwa - BWA-CONTROL-QUANT
bwa mem bwa_index control1_1.fq control1_2.fq > control_reads.sam;

# Align RNA-Seq reads of test experiment to the genome using bwa - BWA-TEST-QUANT
bwa mem bwa_index test1_1.fq test1_2.fq > test_reads.sam;

# Sort and Index the output sam files - SAMTOOLS-CONTROL
samtools view -@ 4 -b control.sam > control.bam 
samtools sort control.bam -@ 4 -o control.sorted.bam;
samtools index control.sorted.bam;

# Sort and Index the output sam files - SAMTOOLS-TEST
samtools view -@ 4 -b test.sam > test.bam 
samtools sort test.bam -@ 4 -o test.sorted.bam;
samtools index test.sorted.bam;

# Predict potential transcripts in control and test BAM files using StringTie - STRINGTIE
stringtie -o control.gtf -G hg19_anno.gtf control.sorted.bam;
stringtie -o test.gtf -G hg19_anno.gtf test.sorted.bam;

# Calculate transcript counts from BAM files using featureCounts - FEATURECOUNTS
featureCounts -p --countReadPairs -B -C -T 4 -a hg19_anno.gtf -o control_counts control.sorted.bam;
featureCounts -p --countReadPairs -B -C -T 4 -a hg19_anno.gtf -o test_counts test.sorted.bam;

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