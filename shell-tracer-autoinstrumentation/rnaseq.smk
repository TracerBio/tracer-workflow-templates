# Snakefile for RNA-Seq analysis pipeline

configfile: "config.yaml"

rule all:
    input:
        # Specify final output files for the entire pipeline
        "results/multiqc/multiqc_report.html",
        "results/deeptools/PCA_rnaseq.png",
        "results/deeptools/fingerprint_rnaseq.png",
        "results/deeptools/differential.bw"

# FASTQC - QC layer 1
rule fastqc:
    input:
        fq1=config["fq1"],
        fq2=config["fq2"]
    output:
        "results/fastqc/{sample}_1_fastqc.html",
        "results/fastqc/{sample}_2_fastqc.html"
    shell:
        """
        fastqc {input.fq1} -o results/fastqc/
        fastqc {input.fq2} -o results/fastqc/
        """

# Creating a human genome index using bowtie2 - BOWTIE2-INDEX
rule bowtie2_index:
    input:
        fa="data/human.fa"
    output:
        "results/bowtie2/human_index.1.bt2"
    params:
        prefix="results/bowtie2/human_index"
    shell:
        """
        bowtie2-build {input.fa} {params.prefix}
        """

# Align RNA-Seq reads to the genome using bowtie2 - BOWTIE2-MAP
rule bowtie2_map:
    input:
        index="results/bowtie2/human_index",
        fq1="data/{sample}_1.fq",
        fq2="data/{sample}_2.fq"
    output:
        "results/bowtie2/{sample}.sam"
    shell:
        """
        bowtie2 --local -x {input.index} -1 {input.fq1} -2 {input.fq2} -S {output}
        """

# Creating a human genome index using STAR - STAR-INDEX
rule star_index:
    input:
        fa="data/human.fa",
        gtf="data/hg19_anno.gtf"
    output:
        directory("results/star/human_star")
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
        index=directory("results/star/human_star"),
        fq1="data/{sample}_1.fq",
        fq2="data/{sample}_2.fq"
    output:
        bam="results/star/{sample}_starAligned.sortedByCoord.out.bam",
        counts="results/star/{sample}_starReadsPerGene.out.tab"
    params:
        prefix="results/star/{sample}_star"
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
        sam="results/bowtie2/{sample}.sam"
    output:
        fa="results/bowtie2/{sample}.fa"
    shell:
        """
        samtools fasta {input.sam} > {output.fa}
        """

# Creating a transcriptome index using kallisto - KALLISTO-INDEX
rule kallisto_index:
    input:
        fa="results/bowtie2/{sample}.fa"
    output:
        index="results/kallisto/{sample}_index"
    shell:
        """
        kallisto index -t 4 -i {output.index} {input.fa}
        """

# Quantify RNA-Seq reads using kallisto - KALLISTO-QUANT
rule kallisto_quant:
    input:
        index="results/kallisto/{sample}_index",
        fq1="data/{sample}_1.fq",
        fq2="data/{sample}_2.fq"
    output:
        directory("results/kallisto/{sample}_quant")
    shell:
        """
        kallisto quant -t 4 -i {input.index} -o {output} {input.fq1} {input.fq2}
        """

# Creating a human genome index using salmon - SALMON-INDEX
rule salmon_index:
    input:
        fa="data/human.fa"
    output:
        directory("results/salmon/salmon_index")
    shell:
        """
        salmon index --threads 4 -t {input.fa} -i {output}
        """

# Quantify RNA-Seq reads using salmon - SALMON-QUANT
rule salmon_quant:
    input:
        index=directory("results/salmon/salmon_index"),
        fq1="data/{sample}_1.fq",
        fq2="data/{sample}_2.fq"
    output:
        directory("results/salmon/{sample}_quant")
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
        fa="data/human.fa"
    output:
        "results/hisat2/hisat2_index.1.ht2"
    params:
        prefix="results/hisat2/hisat2_index"
    shell:
        """
        hisat2-build {input.fa} {params.prefix}
        """

# Align RNA-Seq reads to the genome using hisat2 - HISAT2-ALIGN
rule hisat2_align:
    input:
        index="results/hisat2/hisat2_index",
        fq1="data/{sample}_1.fq",
        fq2="data/{sample}_2.fq"
    output:
        sam="results/hisat2/{sample}_hs2.sam"
    shell:
        """
        hisat2 --fast -x {input.index} -1 {input.fq1} -2 {input.fq2} -S {output.sam}
        """

# Creating a human genome index using bwa - BWA-INDEX
rule bwa_index:
    input:
        fa="data/human.fa"
    output:
        "results/bwa/bwa_index.bwt"
    params:
        prefix="results/bwa/bwa_index"
    shell:
        """
        bwa index -p {params.prefix} {input.fa}
        """

# Align RNA-Seq reads to the genome using bwa - BWA-QUANT
rule bwa_quant:
    input:
        index="results/bwa/bwa_index",
        fq1="data/{sample}_1.fq",
        fq2="data/{sample}_2.fq"
    output:
        sam="results/bwa/{sample}_reads.sam"
    shell:
        """
        bwa mem {input.index} {input.fq1} {input.fq2} > {output.sam}
        """

# Sort and index the output SAM files - SAMTOOLS-SORT-INDEX
rule samtools_sort_index:
    input:
        sam="results/bowtie2/{sample}.sam"
    output:
        bam="results/samtools/{sample}.sorted.bam",
        bai="results/samtools/{sample}.sorted.bam.bai"
    shell:
        """
        samtools view -@ 4 -b {input.sam} > {wildcards.sample}.bam
        samtools sort {wildcards.sample}.bam -@ 4 -o {output.bam}
        samtools index {output.bam}
        rm {wildcards.sample

# Predict potential transcripts in BAM files using StringTie - STRINGTIE
rule stringtie:
    input:
        bam="results/samtools/{sample}.sorted.bam",
        bai="results/samtools/{sample}.sorted.bam.bai",
        gtf="data/hg19_anno.gtf"
    output:
        gtf="results/stringtie/{sample}.gtf"
    shell:
        """
        stringtie {input.bam} -o {output.gtf} -G {input.gtf}
        """

# Calculate transcript counts from BAM files using featureCounts - FEATURECOUNTS
rule featurecounts:
    input:
        bam="results/samtools/{sample}.sorted.bam",
        bai="results/samtools/{sample}.sorted.bam.bai",
        gtf="data/hg19_anno.gtf"
    output:
        counts="results/featurecounts/{sample}_counts.txt"
    shell:
        """
        featureCounts -p --countReadPairs -B -C -T 4 -a {input.gtf} -o {output.counts} {input.bam}
        """

# Collate logs and perform post-mapping QC with MULTIQC - MULTIQC
rule multiqc:
    input:
        expand("results/fastqc/{sample}_1_fastqc.html", sample=config["samples"]),
        expand("results/fastqc/{sample}_2_fastqc.html", sample=config["samples"]),
        expand("results/featurecounts/{sample}_counts.txt", sample=config["samples"])
    output:
        html="results/multiqc/multiqc_report.html"
    shell:
        """
        multiqc results/fastqc results/featurecounts -o results/multiqc/
        """

# Summary of BAM files using deeptools - BAMSUMMARY
rule multiBamSummary:
    input:
        bam1="results/samtools/control.sorted.bam",
        bam2="results/samtools/test.sorted.bam"
    output:
        npz="results/deeptools/rnaseq.npz"
    shell:
        """
        multiBamSummary bins --bamfiles {input.bam1} {input.bam2} -o {output.npz}
        """

# PCA analysis for RNASEQ experiments - PCA
rule plotPCA:
    input:
        npz="results/deeptools/rnaseq.npz"
    output:
        png="results/deeptools/PCA_rnaseq.png"
    shell:
        """
        plotPCA -in {input.npz} -o {output.png}
        """

# RNASEQ data comparison via Fingerprint plots - FINGERPRINT
rule plotFingerprint:
    input:
        bam1="results/samtools/control.sorted.bam",
        bam2="results/samtools/test.sorted.bam"
    output:
        png="results/deeptools/fingerprint_rnaseq.png"
    shell:
        """
        plotFingerprint -b {input.bam1} {input.bam2} --labels Control Test --plotFile {output.png}
        """

# Obtain bigwig files for visualization of data - BAMCOMPARE
rule bamCompare:
    input:
        bam1="results/samtools/test.sorted.bam",
        bam2="results/samtools/control.sorted.bam"
    output:
        bw="results/deeptools/differential.bw"
    shell:
        """
        bamCompare -b1 {input.bam1} -b2 {input.bam2} -o {output.bw} -of bigwig
        """