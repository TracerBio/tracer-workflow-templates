# Snakefile for RNA-Seq analysis pipeline

# Define the data directory
DATA_DIR = "/workspace/tracer-workflow-templates/data"

# Define the samples to be processed
SAMPLES = ["control", "test"]

# Final target rule that depends on the final outputs you want
rule all:
    input:
        # Final outputs for each sample
        expand(f"{DATA_DIR}/{{sample}}.sorted.bam", sample=SAMPLES),
        expand(f"{DATA_DIR}/{{sample}}.sorted.bam.bai", sample=SAMPLES),
        f"{DATA_DIR}/multiqc_report.html",
        f"{DATA_DIR}/rnaseq.npz",
        f"{DATA_DIR}/PCA_rnaseq.png",
        f"{DATA_DIR}/fingerprint_rnaseq.png",
        f"{DATA_DIR}/differential.bw",
        expand(f"{DATA_DIR}/hisat2_index.{{i}}.ht2", i=range(1, 9)),
        expand(f"{DATA_DIR}/bwa_index.{{suffix}}", suffix=["bwt", "pac", "ann", "amb", "sa"]),
        expand(f"{DATA_DIR}/salmon_{{sample}}", sample=SAMPLES),
        
# Creating a human genome index using STAR - STAR-INDEX
rule star_index:
    input:
        fa=f"{DATA_DIR}/human.fa",
        gtf=f"{DATA_DIR}/hg19_anno.gtf"
    output:
        directory(f"{DATA_DIR}/human_star")
    params:
        sjdbOverhang=99
    shell:
        """
        cd {DATA_DIR}
        STAR --runThreadN 4 --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeSAindexNbases 10 \
             --genomeFastaFiles {input.fa} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.sjdbOverhang}
        """

# Align RNA-Seq reads to the genome using STAR - STAR-MAP
rule star_map:
    input:
        index=f"{DATA_DIR}/human_star",
        fq1=f"{DATA_DIR}/{{sample}}1_1.fq",
        fq2=f"{DATA_DIR}/{{sample}}1_2.fq"
    output:
        bam=f"{DATA_DIR}/{{sample}}_starAligned.sortedByCoord.out.bam",
        counts=f"{DATA_DIR}/{{sample}}_starReadsPerGene.out.tab"
    params:
        prefix=f"{DATA_DIR}/{{sample}}_star"
    shell:
        """
        cd {DATA_DIR}
        STAR --runThreadN 4 --genomeDir {input.index} \
             --readFilesIn {input.fq1} {input.fq2} \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts
        """

# Convert SAM files into transcriptome FASTA files - SAM2FASTA
rule samtools_fasta:
    input:
        sam=f"{DATA_DIR}/{{sample}}.sam"
    output:
        fa=f"{DATA_DIR}/{{sample}}.fa"
    shell:
        """
        cd {DATA_DIR}
        samtools fasta {input.sam} > {output.fa}
        """

# Creating a human genome index using salmon - SALMON-INDEX
rule salmon_index:
    input:
        fa=f"{DATA_DIR}/human.fa"
    output:
        directory(f"{DATA_DIR}/salmon_index")
    shell:
        """
        cd {DATA_DIR}
        salmon index --threads 4 -t {input.fa} -i {output}
        """

# Quantify RNA-Seq reads using salmon - SALMON-QUANT
rule salmon_quant:
    input:
        index=f"{DATA_DIR}/salmon_index",
        fq1=f"{DATA_DIR}/{{sample}}1_1.fq",
        fq2=f"{DATA_DIR}/{{sample}}1_2.fq"
    output:
        directory(f"{DATA_DIR}/salmon_{{sample}}")
    shell:
        """
        cd {DATA_DIR}
        salmon quant --threads 4 --libType=U \
                     -i {input.index} \
                     -1 {input.fq1} -2 {input.fq2} \
                     -o {output}
        """

# Creating a human genome index using hisat2 - HISAT2-INDEX
rule hisat2_index:
    input:
        fa=f"{DATA_DIR}/human.fa"
    output:
        expand(f"{DATA_DIR}/hisat2_index.{{i}}.ht2", i=range(1, 9))
    params:
        prefix=f"{DATA_DIR}/hisat2_index"
    shell:
        """
        cd {DATA_DIR}
        hisat2-build {input.fa} {params.prefix}
        """

# Align RNA-Seq reads to the genome using hisat2 - HISAT2-ALIGN
rule hisat2_align:
    input:
        index=f"{DATA_DIR}/hisat2_index",
        fq1=f"{DATA_DIR}/{{sample}}1_1.fq",
        fq2=f"{DATA_DIR}/{{sample}}1_2.fq"
    output:
        sam=f"{DATA_DIR}/{{sample}}_hs2.sam"
    shell:
        """
        cd {DATA_DIR}
        hisat2 --fast -x {input.index} -1 {input.fq1} -2 {input.fq2} -S {output.sam}
        """

# Creating a human genome index using bwa - BWA-INDEX
rule bwa_index:
    input:
        fa=f"{DATA_DIR}/human.fa"
    output:
        expand(f"{DATA_DIR}/bwa_index.{{suffix}}", suffix=["bwt", "pac", "ann", "amb", "sa"])
    params:
        prefix=f"{DATA_DIR}/bwa_index"
    shell:
        """
        cd {DATA_DIR}
        bwa index -p {params.prefix} {input.fa}
        """

# Align RNA-Seq reads to the genome using bwa - BWA-QUANT
rule bwa_quant:
    input:
        index=f"{DATA_DIR}/bwa_index",
        fq1=f"{DATA_DIR}/{{sample}}1_1.fq",
        fq2=f"{DATA_DIR}/{{sample}}1_2.fq"
    output:
        sam=f"{DATA_DIR}/{{sample}}_reads.sam"
    shell:
        """
        cd {DATA_DIR}
        bwa mem {input.index} {input.fq1} {input.fq2} > {output.sam}
        """

# Sort and index the output SAM files - SAMTOOLS-SORT-INDEX
rule samtools_sort_index:
    input:
        sam=f"{DATA_DIR}/{{sample}}.sam"
    output:
        bam=f"{DATA_DIR}/{{sample}}.sorted.bam",
        bai=f"{DATA_DIR}/{{sample}}.sorted.bam.bai"
    shell:
        """
        cd {DATA_DIR}
        samtools view -@ 4 -b {input.sam} > {wildcards.sample}.bam
        samtools sort {wildcards.sample}.bam -@ 4 -o {output.bam}
        samtools index {output.bam}
        rm {wildcards.sample}.bam
        """

# Calculate transcript counts from BAM files using featureCounts - FEATURECOUNTS
rule featurecounts:
    input:
        bam=f"{DATA_DIR}/{{sample}}.sorted.bam",
        bai=f"{DATA_DIR}/{{sample}}.sorted.bam.bai",
        gtf=f"{DATA_DIR}/hg19_anno.gtf"
    output:
        counts=f"{DATA_DIR}/{{sample}}_counts.txt"
    shell:
        """
        cd {DATA_DIR}
        featureCounts -p --countReadPairs -B -C -T 4 -a {input.gtf} -o {output.counts} {input.bam}
        """

# Collate logs and perform post-mapping QC with MultiQC - MULTIQC
rule multiqc:
    input:
        expand(f"{DATA_DIR}/{{sample}}_starAligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand(f"{DATA_DIR}/{{sample}}_counts.txt", sample=SAMPLES),
        expand(f"{DATA_DIR}/{{sample}}.sorted.bam", sample=SAMPLES),
        expand(f"{DATA_DIR}/{{sample}}.sorted.bam.bai", sample=SAMPLES)
    output:
        html=f"{DATA_DIR}/multiqc_report.html"
    shell:
        """
        cd {DATA_DIR}
        multiqc {DATA_DIR} -o {DATA_DIR} --force
        """

# Summary of BAM files using deepTools - BAMSUMMARY
rule multiBamSummary:
    input:
        bam1=f"{DATA_DIR}/control.sorted.bam",
        bam2=f"{DATA_DIR}/test.sorted.bam"
    output:
        npz=f"{DATA_DIR}/rnaseq.npz"
    shell:
        """
        cd {DATA_DIR}
        multiBamSummary bins --bamfiles {input.bam1} {input.bam2} -o {output.npz}
        """

# PCA analysis for RNA-Seq experiments - PCA
rule plotPCA:
    input:
        npz=f"{DATA_DIR}/rnaseq.npz"
    output:
        png=f"{DATA_DIR}/PCA_rnaseq.png"
    shell:
        """
        cd {DATA_DIR}
        plotPCA -in {input.npz} -o {output.png}
        """

# RNA-Seq data comparison via Fingerprint plots - FINGERPRINT
rule plotFingerprint:
    input:
        bam1=f"{DATA_DIR}/control.sorted.bam",
        bam2=f"{DATA_DIR}/test.sorted.bam"
    output:
        png=f"{DATA_DIR}/fingerprint_rnaseq.png"
    shell:
        """
        cd {DATA_DIR}
        plotFingerprint -b {input.bam1} {input.bam2} --labels Control Test --plotFile {output.png}
        """

# Obtain bigwig files for visualization of data - BAMCOMPARE
rule bamCompare:
    input:
        bam1=f"{DATA_DIR}/test.sorted.bam",
        bam2=f"{DATA_DIR}/control.sorted.bam"
    output:
        bw=f"{DATA_DIR}/differential.bw"
    shell:
        """
        cd {DATA_DIR}
        bamCompare -b1 {input.bam1} -b2 {input.bam2} -o {output.bw} -of bigwig
        """