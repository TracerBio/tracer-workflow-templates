#!/bin/bash

# FASTQC - QC layer 1
cd /workspace/tracer-workflow-templates/data

# Step 2: Trimming adapters and low-quality bases using Trimmomatic
for fq in raw_data/*_R1.fastq.gz; do
    sample=$(basename $fq _R1.fastq.gz)
    trimmomatic PE -phred33 \
        raw_data/${sample}_R1.fastq.gz raw_data/${sample}_R2.fastq.gz \
        trimmed_data/${sample}_R1_paired.fastq.gz trimmed_data/${sample}_R1_unpaired.fastq.gz \
        trimmed_data/${sample}_R2_paired.fastq.gz trimmed_data/${sample}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
done

# Step 3: Alignment to the reference genome using BWA
for fq in trimmed_data/*_R1_paired.fastq.gz; do
    sample=$(basename $fq _R1_paired.fastq.gz)
    bwa mem -t 8 reference_genome.fa \
        trimmed_data/${sample}_R1_paired.fastq.gz trimmed_data/${sample}_R2_paired.fastq.gz > aligned_data/${sample}.sam
done

# Step 4: Convert SAM to BAM, sort, and index using SAMtools
for sam in aligned_data/*.sam; do
    sample=$(basename $sam .sam)
    samtools view -Sb $sam | samtools sort -o bam_files/${sample}.sorted.bam
    samtools index bam_files/${sample}.sorted.bam
done

# Step 5: Mark duplicates using Picard
for bam in bam_files/*.sorted.bam; do
    sample=$(basename $bam .sorted.bam)
    picard MarkDuplicates I=$bam O=bam_files/${sample}.dedup.bam M=logs/${sample}.metrics.txt
    samtools index bam_files/${sample}.dedup.bam
done

# Step 6: Base Quality Score Recalibration (BQSR) using GATK
for bam in bam_files/*.dedup.bam; do
    sample=$(basename $bam .dedup.bam)
    gatk BaseRecalibrator \
        -I $bam \
        -R reference_genome.fa \
        --known-sites known_variants.vcf \
        -O logs/${sample}.recal.table
    
    gatk ApplyBQSR \
        -R reference_genome.fa \
        -I $bam \
        --bqsr-recal-file logs/${sample}.recal.table \
        -O bam_files/${sample}.recal.bam
    
    samtools index bam_files/${sample}.recal.bam
done

# Step 7: Variant Calling using GATK HaplotypeCaller
for bam in bam_files/*.recal.bam; do
    sample=$(basename $bam .recal.bam)
    gatk HaplotypeCaller \
        -R reference_genome.fa \
        -I $bam \
        -O vcf_files/${sample}.g.vcf.gz \
        -ERC GVCF
done

# Step 8: Joint Genotyping using GATK GenotypeGVCFs
gatk GenotypeGVCFs \
    -R reference_genome.fa \
    -V vcf_files/cohort.g.vcf.gz \
    -O vcf_files/cohort.vcf.gz

# Step 9: Variant Annotation using ANNOVAR or SnpEff
for vcf in vcf_files/*.vcf.gz; do
    sample=$(basename $vcf .vcf.gz)
    snpEff ann \
        -v GRCh38.86 \
        $vcf > vcf_files/${sample}.annotated.vcf
done

# Step 10: CNV Analysis using CNVkit
for bam in bam_files/*.recal.bam; do
    sample=$(basename $bam .recal.bam)
    cnvkit.py batch $bam \
        -n \
        -f reference_genome.fa \
        -t exome_targets.bed \
        -m wgs \
        --output-dir cnv_results/ \
        --scatter --diagram
done

# Final message
echo "Pipeline completed. Results are in the respective directories."