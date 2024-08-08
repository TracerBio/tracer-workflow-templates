# Use the data directory from the configuration file
DATA_DIR = "/workspace/tracer-workflow-templates/data"

# Final target rule that depends on the final outputs you want
rule all:
    input:
        expand(f"{DATA_DIR}/human_index.{{i}}.bt2", i=range(1, 4)),
        expand(f"{DATA_DIR}/human_index.rev.{{i}}.bt2", i=range(1, 2)),
        expand(f"{DATA_DIR}/hisat2_index.{{i}}.ht2", i=range(1, 9)),
        expand(f"{DATA_DIR}/bwa_index.{{suffix}}", suffix=["bwt", "pac", "ann", "amb", "sa"]),
        f"{DATA_DIR}/control.sorted.bam",
        f"{DATA_DIR}/test.sorted.bam",
        f"{DATA_DIR}/control.sorted.bam.bai",
        f"{DATA_DIR}/test.sorted.bam.bai",
        f"{DATA_DIR}/control_counts.txt",
        f"{DATA_DIR}/test_counts.txt",
        f"{DATA_DIR}/rnaseq.npz",
        f"{DATA_DIR}/PCA_rnaseq.png",
        f"{DATA_DIR}/fingerprint_rnaseq.png",
        f"{DATA_DIR}/differential.bw",
        f"{DATA_DIR}/edger_results.csv",
        f"{DATA_DIR}/diff_heatmap.png"

# Rule for creating a STAR genome index
rule star_index:
    input:
        fasta=f"{DATA_DIR}/human.fa",
        gtf=f"{DATA_DIR}/hg19_anno.gtf"
    output:
        index_dir=directory(f"{DATA_DIR}/human_star")
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate \
            --genomeDir {output.index_dir} --genomeSAindexNbases 10 \
            --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 99
        """

# Rule for aligning control RNA-Seq reads using STAR
rule star_align_control:
    input:
        index_dir=f"{DATA_DIR}/human_star",
        reads1=f"{DATA_DIR}/control1_1.fq",
        reads2=f"{DATA_DIR}/control1_2.fq"
    output:
        bam=f"{DATA_DIR}/control.sorted.bam"
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index_dir} \
            --readFilesIn {input.reads1} {input.reads2} \
            --outFileNamePrefix control_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts;
        mv control_starAligned.sortedByCoord.out.bam {output.bam}
        """

# Rule for aligning test RNA-Seq reads using STAR
rule star_align_test:
    input:
        index_dir=f"{DATA_DIR}/human_star",
        reads1=f"{DATA_DIR}/test1_1.fq",
        reads2=f"{DATA_DIR}/test1_2.fq"
    output:
        bam=f"{DATA_DIR}/test.sorted.bam"
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index_dir} \
            --readFilesIn {input.reads1} {input.reads2} \
            --outFileNamePrefix test_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts;
        mv test_starAligned.sortedByCoord.out.bam {output.bam}
        """

# Rule for indexing control BAM files with SAMtools
rule samtools_index_control:
    input:
        bam=f"{DATA_DIR}/control.sorted.bam"
    output:
        bai=f"{DATA_DIR}/control.sorted.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

# Rule for indexing test BAM files with SAMtools
rule samtools_index_test:
    input:
        bam=f"{DATA_DIR}/test.sorted.bam"
    output:
        bai=f"{DATA_DIR}/test.sorted.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

# Rule for creating a bowtie2 genome index
rule bowtie2_index:
    input:
        fasta=f"{DATA_DIR}/human.fa"
    output:
        expand(f"{DATA_DIR}/human_index.{{i}}.bt2", i=range(1, 4)),
        expand(f"{DATA_DIR}/human_index.rev.{{i}}.bt2", i=range(1, 2)),
    params:
        prefix=f"{DATA_DIR}/human_index"
    shell:
        """
        bowtie2-build {input.fasta} {params.prefix}
        """

# Rule for aligning control RNA-Seq reads using bowtie2
rule bowtie2_align_control:
    input:
        index=f"{DATA_DIR}/human_index",
        reads1=f"{DATA_DIR}/control1_1.fq",
        reads2=f"{DATA_DIR}/control1_2.fq"
    output:
        sam=f"{DATA_DIR}/control_bowtie2.sam"
    threads: 4
    shell:
        """
        bowtie2 --local -x {input.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam}
        """

# Rule for aligning test RNA-Seq reads using bowtie2
rule bowtie2_align_test:
    input:
        index=f"{DATA_DIR}/human_index",
        reads1=f"{DATA_DIR}/test1_1.fq",
        reads2=f"{DATA_DIR}/test1_2.fq"
    output:
        sam=f"{DATA_DIR}/test_bowtie2.sam"
    threads: 4
    shell:
        """
        bowtie2 --local -x {input.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam}
        """
# Rule for creating a hisat2 genome index
rule hisat2_index:
    input:
        fasta=f"{DATA_DIR}/human.fa"
    output:
        expand(f"{DATA_DIR}/hisat2_index.{{i}}.ht2", i=range(1, 9))
    params:
        prefix=f"{DATA_DIR}/hisat2_index"
    shell:
        """
        hisat2-build {input.fasta} {params.prefix}
        """

# Rule for aligning control RNA-Seq reads using hisat2
rule hisat2_align_control:
    input:
        index=f"{DATA_DIR}/hisat2_index",
        reads1=f"{DATA_DIR}/control1_1.fq",
        reads2=f"{DATA_DIR}/control1_2.fq"
    output:
        sam=f"{DATA_DIR}/control_hisat2.sam"
    threads: 4
    shell:
        """
        hisat2 --fast -x {input.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam}
        """

# Rule for aligning test RNA-Seq reads using hisat2
rule hisat2_align_test:
    input:
        index=f"{DATA_DIR}/hisat2_index",
        reads1=f"{DATA_DIR}/test1_1.fq",
        reads2=f"{DATA_DIR}/test1_2.fq"
    output:
        sam=f"{DATA_DIR}/test_hisat2.sam"
    threads: 4
    shell:
        """
        hisat2 --fast -x {input.index} -1 {input.reads1} -2 {input.reads2} -S {output.sam}
        """

# Rule for creating a bwa genome index
rule bwa_index:
    input:
        fasta=f"{DATA_DIR}/human.fa"
    output:
        # Correctly defining BWA index files
        expand(f"{DATA_DIR}/bwa_index.{{suffix}}", suffix=["bwt", "pac", "ann", "amb", "sa"])
    params:
        prefix=f"{DATA_DIR}/bwa_index"
    shell:
        """
        bwa index -p {params.prefix} {input.fasta}
        """

# Rule for aligning control RNA-Seq reads using bwa
rule bwa_align_control:
    input:
        index=f"{DATA_DIR}/bwa_index",
        reads1=f"{DATA_DIR}/control1_1.fq",
        reads2=f"{DATA_DIR}/control1_2.fq"
    output:
        sam=f"{DATA_DIR}/control_bwa.sam"
    shell:
        """
        bwa mem {input.index} {input.reads1} {input.reads2} > {output.sam}
        """

# Rule for aligning test RNA-Seq reads using bwa
rule bwa_align_test:
    input:
        index=f"{DATA_DIR}/bwa_index",
        reads1=f"{DATA_DIR}/test1_1.fq",
        reads2=f"{DATA_DIR}/test1_2.fq"
    output:
        sam=f"{DATA_DIR}/test_bwa.sam"
    shell:
        """
        bwa mem {input.index} {input.reads1} {input.reads2} > {output.sam}
        """

# Rule for calculating transcript counts using featureCounts for control
rule featurecounts_control:
    input:
        bam=f"{DATA_DIR}/control.sorted.bam",
        gtf=f"{DATA_DIR}/hg19_anno.gtf"
    output:
        counts=f"{DATA_DIR}/control_counts.txt"
    threads: 4
    shell:
        """
        featureCounts -p --countReadPairs -B -C -T {threads} -a {input.gtf} -o {output.counts} {input.bam}
        """

# Rule for calculating transcript counts using featureCounts for test
rule featurecounts_test:
    input:
        bam=f"{DATA_DIR}/test.sorted.bam",
        gtf=f"{DATA_DIR}/hg19_anno.gtf"
    output:
        counts=f"{DATA_DIR}/test_counts.txt"
    threads: 4
    shell:
        """
        featureCounts -p --countReadPairs -B -C -T {threads} -a {input.gtf} -o {output.counts} {input.bam}
        """

# Rule for collating logs and performing post-mapping QC with MultiQC
rule multiqc:
    input:
        f"{DATA_DIR}/control_counts.txt",
        f"{DATA_DIR}/test_counts.txt",
        f"{DATA_DIR}/control.sorted.bam",
        f"{DATA_DIR}/test.sorted.bam"
    output:
        report=directory(f"{DATA_DIR}/multiqc_report")  # Specify as directory
    shell:
        """
        multiqc {DATA_DIR} -o {output.report}
        """

# Rule for generating summary of BAM files using deeptools
rule deeptools_summary:
    input:
        bam_control=f"{DATA_DIR}/control.sorted.bam",
        bam_test=f"{DATA_DIR}/test.sorted.bam",
        bai_control=f"{DATA_DIR}/control.sorted.bam.bai",
        bai_test=f"{DATA_DIR}/test.sorted.bam.bai"
    output:
        npz=f"{DATA_DIR}/rnaseq.npz"
    log:
        log=f"{DATA_DIR}/logs/deeptools_summary.log",
        err=f"{DATA_DIR}/logs/deeptools_summary.err"
    shell:
        """
        mkdir -p {DATA_DIR}/logs
        multiBamSummary bins --bamfiles {input.bam_control} {input.bam_test} -o {output.npz} > {log.log} 2> {log.err}
        """

# Rule for PCA analysis using deeptools
rule pca_plot:
    input:
        npz=f"{DATA_DIR}/rnaseq.npz"
    output:
        pca_png=f"{DATA_DIR}/PCA_rnaseq.png"
    shell:
        """
        plotPCA -in {input.npz} -o {output.pca_png}
        """

# Rule for RNA-Seq data comparison via Fingerprint plots using deeptools
rule fingerprint_plot:
    input:
        bam_control=f"{DATA_DIR}/control.sorted.bam",
        bam_test=f"{DATA_DIR}/test.sorted.bam",
        bai_control=f"{DATA_DIR}/control.sorted.bam.bai",
        bai_test=f"{DATA_DIR}/test.sorted.bam.bai"
    output:
        fingerprint_png=f"{DATA_DIR}/fingerprint_rnaseq.png"
    shell:
        """
        plotFingerprint -b {input.bam_control} {input.bam_test} --labels Control Test --plotFile {output.fingerprint_png}
        """

# Rule for obtaining bigwig files for visualization using deeptools
rule bamcompare:
    input:
        bam_control=f"{DATA_DIR}/control.sorted.bam",
        bam_test=f"{DATA_DIR}/test.sorted.bam",
        bai_control=f"{DATA_DIR}/control.sorted.bam.bai",
        bai_test=f"{DATA_DIR}/test.sorted.bam.bai"
    output:
        bw=f"{DATA_DIR}/differential.bw"
    shell:
        """
        bamCompare -b1 {input.bam_test} -b2 {input.bam_control} -o {output.bw} -of bigwig
        """

# Rule for running the R script for edgeR analysis and heatmap creation
rule run_deseq2_r:
    input:
        counts_control=f"{DATA_DIR}/control_counts.txt",
        counts_test=f"{DATA_DIR}/test_counts.txt",
        diff_txt=f"{DATA_DIR}/diff.txt"  # Assuming diff.txt is in the data directory
    output:
        csv=f"{DATA_DIR}/edger_results.csv",
        heatmap_png=f"{DATA_DIR}/diff_heatmap.png"
    shell:
        """
        Rscript ../data/DESeq2.r
        """