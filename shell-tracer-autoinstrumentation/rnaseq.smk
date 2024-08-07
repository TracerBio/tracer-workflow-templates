# Snakefile for RNA-Seq analysis pipeline

# Use the data directory from the configuration file
DATA_DIR = "/workspace/tracer-workflow-templates/data"

# Define the samples that will be processed
SAMPLES = ["control"]

# Final target rule that depends on the final outputs you want
rule all:
    input:
        f"{DATA_DIR}/control.sorted.bam",
        f"{DATA_DIR}/control.sorted.bam.bai",
        f"{DATA_DIR}/control.fa",
        f"{DATA_DIR}/control_counts.txt",
        f"{DATA_DIR}/rnaseq.npz",
        f"{DATA_DIR}/PCA_rnaseq.png",
        f"{DATA_DIR}/fingerprint_rnaseq.png",
        f"{DATA_DIR}/differential.bw",
        f"{DATA_DIR}/multiqc_report.html",
        expand(f"{DATA_DIR}/human_index.{{i}}.bt2", i=range(1, 4)),
        expand(f"{DATA_DIR}/human_index.rev.{{i}}.bt2", i=range(1, 2)),
        expand(f"{DATA_DIR}/hisat2_index.{{i}}.ht2", i=range(1, 9)),
        expand(f"{DATA_DIR}/bwa_index.{{suffix}}", suffix=["bwt", "pac", "ann", "amb", "sa"]),
        f"{DATA_DIR}/control_starAligned.sortedByCoord.out.bam",
        f"{DATA_DIR}/control_starAligned.sortedByCoord.out.bam",
        f"{DATA_DIR}/control_starReadsPerGene.out.tab",
        directory(f"{DATA_DIR}/salmon_index"),
        directory(f"{DATA_DIR}/control_salmon_quant"),
        
# Creating a human genome index using bowtie2 - BOWTIE2-INDEX
rule bowtie2_index:
    input:
        fa=f"{DATA_DIR}/human.fa"
    output:
        expand(f"{DATA_DIR}/human_index.{{i}}.bt2", i=range(1, 4)),
        expand(f"{DATA_DIR}/human_index.rev.{{i}}.bt2", i=range(1, 2)),
    params:
        prefix=f"{DATA_DIR}/human_index"
    shell:
        """
        bowtie2-build {input.fa} {params.prefix}
        """

# Align RNA-Seq reads to the genome using bowtie2 - BOWTIE2-MAP
rule bowtie2_map:
    input:
        index=f"{DATA_DIR}/human_index",
        fq1=f"{DATA_DIR}/control1_1.fq",
        fq2=f"{DATA_DIR}/control1_2.fq"
    output:
        f"{DATA_DIR}/{{sample}}.sam"
    shell:
        """
        bowtie2 --local -x {input.index} -1 {input.fq1} -2 {input.fq2} -S {output}
        """

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
        salmon index --threads 4 -t {input.fa} -i {output}
        """

# Quantify RNA-Seq reads using salmon - SALMON-QUANT
rule salmon_quant:
    input:
        index=f"{DATA_DIR}/salmon_index",
        fq1=f"{DATA_DIR}/{{sample}}1_1.fq",
        fq2=f"{DATA_DIR}/{{sample}}1_2.fq"
    output:
        directory(f"{DATA_DIR}/{{sample}}_salmon_quant")
    shell:
        """
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
        hisat2 --fast -x {input.index} -1 {input.fq1} -2 {input.fq2} -S {output}
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
        featureCounts -p --countReadPairs -B -C -T 4 -a {input.gtf} -o {output.counts} {input.bam}
        """

# Collate logs and perform post-mapping QC with MultiQC - MULTIQC
rule multiqc:
    input:
        # Collect all relevant input files for MultiQC
        expand(f"{DATA_DIR}/{{sample}}.sam", sample=SAMPLES),
        expand(f"{DATA_DIR}/{{sample}}_starAligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand(f"{DATA_DIR}/{{sample}}_counts.txt", sample=SAMPLES),
        expand(f"{DATA_DIR}/{{sample}}.sorted.bam", sample=SAMPLES),
        expand(f"{DATA_DIR}/{{sample}}.sorted.bam.bai", sample=SAMPLES)
    output:
        html=f"{DATA_DIR}/multiqc_report.html"
    shell:
        """
        multiqc {DATA_DIR} -o {DATA_DIR} --force
        """

# Summary of BAM files using deepTools - multiBamSummary
rule multiBamSummary:
    input:
        bam=f"{DATA_DIR}/control.sorted.bam"
    output:
        npz=f"{DATA_DIR}/rnaseq.npz"
    shell:
        """
        multiBamSummary bins --bamfiles {input.bam} -o {output.npz}
        """

# PCA analysis for RNA-Seq experiments - PCA
rule plotPCA:
    input:
        npz=f"{DATA_DIR}/rnaseq.npz"
    output:
        png=f"{DATA_DIR}/PCA_rnaseq.png"
    shell:
        """
        plotPCA --corData {input.npz} --plotFile {output.png}
        """

# RNA-Seq data comparison via Fingerprint plots - FINGERPRINT
rule plotFingerprint:
    input:
        bam=f"{DATA_DIR}/control.sorted.bam"
    output:
        png=f"{DATA_DIR}/fingerprint_rnaseq.png"
    shell:
        """
        plotFingerprint -b {input.bam} --labels Control --plotFile {output.png}
        """

# Obtain bigwig files for visualization of data - BAMCOMPARE
rule bamCompare:
    input:
        bam1=f"{DATA_DIR}/control.sorted.bam",
        bam2=f"{DATA_DIR}/control.sorted.bam"  # Ideally you want a second sample, like test.sorted.bam
    output:
        bw=f"{DATA_DIR}/differential.bw"
    shell:
        """
        bamCompare -b1 {input.bam1} -b2 {input.bam2} -o {output.bw} -of bigwig
        """