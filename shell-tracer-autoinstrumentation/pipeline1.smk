# Rule for creating a STAR genome index
rule star_index:
    input:
        fasta="human.fa",
        gtf="hg19_anno.gtf"
    output:
        index_dir=directory("human_star")
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate \
            --genomeDir {output.index_dir} --genomeSAindexNbases 10 \
            --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 99
        """

# Rule for aligning control RNA-Seq reads using STAR
rule star_align:
    input:
        index_dir="human_star",
        reads1="control_1.fq",
        reads2="control_2.fq"
    output:
        bam="control.sorted.bam"
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index_dir} \
            --readFilesIn {input.reads1} {input.reads2} \
            --outFileNamePrefix {wildcards.sample}_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts;
        mv {wildcards.sample}_starAligned.sortedByCoord.out.bam {output.bam}
        """

# Rule for indexing BAM files with SAMtools
rule samtools_index:
    input:
        bam="control.sorted.bam"
    output:
        bai="control.sorted.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

# Rule for aligning test RNA-Seq reads using STAR
rule star_align:
    input:
        index_dir="human_star",
        reads1="test_1.fq",
        reads2="test_2.fq"
    output:
        bam="test.sorted.bam"
    threads: 4
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index_dir} \
            --readFilesIn {input.reads1} {input.reads2} \
            --outFileNamePrefix {wildcards.sample}_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts;
        mv {wildcards.sample}_starAligned.sortedByCoord.out.bam {output.bam}
        """

# Rule for indexing BAM files with SAMtools
rule samtools_index:
    input:
        bam="test.sorted.bam"
    output:
        bai="test.sorted.bam.bai"
    shell:
        """
        samtools index {input.bam}
        """

# Rule for predicting transcripts using StringTie
rule stringtie:
    input:
        bam="control.sorted.bam",
        gtf="hg19_anno.gtf"
    output:
        gtf="control.gtf"
    shell:
        """
        stringtie -o {output.gtf} -G {input.gtf} {input.bam}
        """

# Rule for calculating transcript counts using featureCounts
rule featurecounts:
    input:
        bam="control.sorted.bam",
        gtf="hg19_anno.gtf"
    output:
        counts="control_counts.txt"
    threads: 4
    shell:
        """
        featureCounts -p --countReadPairs -B -C -T {threads} -a {input.gtf} -o {output.counts} {input.bam}
        """

# Rule for predicting transcripts using StringTie
rule stringtie:
    input:
        bam="test.sorted.bam",
        gtf="hg19_anno.gtf"
    output:
        gtf="test.gtf"
    shell:
        """
        stringtie -o {output.gtf} -G {input.gtf} {input.bam}
        """

# Rule for calculating transcript counts using featureCounts
rule featurecounts:
    input:
        bam="test.sorted.bam",
        gtf="hg19_anno.gtf"
    output:
        counts="test_counts.txt"
    threads: 4
    shell:
        """
        featureCounts -p --countReadPairs -B -C -T {threads} -a {input.gtf} -o {output.counts} {input.bam}
        """

# Rule for collating logs and performing post-mapping QC with MultiQC
rule multiqc:
    input:
        expand("{sample}_counts.txt", sample=["control", "test"]),
        expand("{sample}.gtf", sample=["control", "test"]),
        expand("{sample}.sorted.bam", sample=["control", "test"])
    output:
        report="multiqc_report.html"
    shell:
        """
        multiqc .
        """

# Rule for generating summary of BAM files using deeptools
rule deeptools_summary:
    input:
        bam_control="control.sorted.bam",
        bam_test="test.sorted.bam"
    output:
        npz="rnaseq.npz"
    shell:
        """
        multiBamSummary bins --bamfiles {input.bam_control} {input.bam_test} -o {output.npz}
        """

# Rule for PCA analysis using deeptools
rule pca_plot:
    input:
        npz="rnaseq.npz"
    output:
        pca_png="PCA_rnaseq.png"
    shell:
        """
        plotPCA -in {input.npz} -o {output.pca_png}
        """

# Rule for RNA-Seq data comparison via Fingerprint plots using deeptools
rule fingerprint_plot:
    input:
        bam_control="control.sorted.bam",
        bam_test="test.sorted.bam"
    output:
        fingerprint_png="fingerprint_rnaseq.png"
    shell:
        """
        plotFingerprint -b {input.bam_control} {input.bam_test} --labels Control Test --plotFile {output.fingerprint_png}
        """

# Rule for obtaining bigwig files for visualization using deeptools
rule bamcompare:
    input:
        bam_control="control.sorted.bam",
        bam_test="test.sorted.bam"
    output:
        bw="differential.bw"
    shell:
        """
        bamCompare -b1 {input.bam_test} -b2 {input.bam_control} -o {output.bw} -of bigwig
        """

# Rule for running the R script for edgeR analysis and heatmap creation
rule run_deseq2_r:
    input:
        counts_control="control_counts.txt",
        counts_test="test_counts.txt",
        diff_txt="diff.txt"
    output:
        csv="edger_results.csv",
        heatmap_png="diff_heatmap.png"
    shell:
        """
        Rscript ./DESeq2.r
        """