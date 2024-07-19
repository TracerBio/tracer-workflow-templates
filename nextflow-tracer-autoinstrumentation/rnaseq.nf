nextflow.enable.dsl=2

process Bowtie2Index {
    input:
    path genome_fa

    output:
    path "human_index.*.bt2"
    path "human_index.rev.*.bt2"

    script:
    """
    bowtie2-build $genome_fa human_index
    """
}

process Bowtie2MapCONTROL {
    tag "$sample_id1"
    input:
    path(control1)
    path(control2) 
    path "human_index.*.bt2"
    path "human_index.rev.*.bt2"

    output:
    path ("control.sam")

    script:
    """
    bowtie2 --local -x human_index -1 $control1 -2 $control2 -S control.sam
    """
}

process Bowtie2MapTEST {
    tag "$sample_id2"
    input:
    path(test1)
    path(test2)
    path "human_index.*.bt2"
    path "human_index.rev.*.bt2"

    output:
    path ("test.sam")

    script:
    """
    bowtie2 --local -x human_index -1 $test1 -2 $test2 -S test.sam
    """
}

process STARIndex {
    input:
    path genome_fa
    path genome_gtf

    output:
    path ("human_star")

    script:
    """
    STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./human_star --genomeSAindexNbases 10 --genomeFastaFiles $genome_fa --sjdbGTFfile $genome_gtf --sjdbOverhang 99
    """
}

process STARMapCONTROL {
    tag "$sample_id1"
    input:
    path(control1)
    path(control2)
    path("human_star")

    output:
    path "control_starAligned.sortedByCoord.out.bam"
    path "control_starReadsPerGene.out.tab"

    script:
    """
    STAR --runThreadN 4 --genomeDir ./human_star --readFilesIn $control1 $control2 --outFileNamePrefix control_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
    """
}

process STARMapTEST {
    tag "$sample_id2"
    input:
    path(test1)
    path(test2)
    path("human_star")

    output:
    path "test_starAligned.sortedByCoord.out.bam"
    path "test_starReadsPerGene.out.tab"

    script:
    """
    STAR --runThreadN 4 --genomeDir ./human_star --readFilesIn $test1 $test2 --outFileNamePrefix test_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
    """
}

process SamToFastaCONTROL {
    tag "$sample_id1"
    input:
    path("control.sam")

    output:
    path "control.fa"

    script:
    """
    samtools fasta control.sam > control.fa
    """
}

process SamToFastaTEST {
    tag "$sample_id2"
    input:
    path ("test.sam")

    output:
    path "test.fa"

    script:
    """
    samtools fasta test.sam > test.fa
    """
}

process KallistoIndexCONTROL {
    tag "$sample_id1"
    input:
    path("control.fa")

    output:
    path ("control_index.idx")

    script:
    """
    kallisto index -t 4 -i control_index.idx control.fa
    """
}

process KallistoIndexTEST {
    tag "$sample_id2"
    input:
    path ("test.fa")

    output:
    path ("test_index.idx")

    script:
    """
    kallisto index -t 4 -i test_index.idx test.fa
    """
}

process KallistoQuantCONTROL {
    tag "$sample_id1"
    input:
    path(control1)
    path(control2)
    path ("control_index.idx")

    output:
    path "control_quant"

    script:
    """
    kallisto quant -t 4 -i control_index.idx -o control_quant $control1 $control2
    """
}

process KallistoQuantTEST {
    tag "$sample_id2"
    input:
    path(test1)
    path(test2)
    path ("test_index.idx")
    
    output:
    path "test_quant"

    script:
    """
    kallisto quant -t 4 -i test_index.idx -o test_quant $test1 $test2
    """
}

process SalmonIndex {
    input:
    path genome_fa

    output:
    path "salmon_index/*"

    script:
    """
    salmon index --threads 4 -t $genome_fa -i salmon_index
    """
}

process SalmonQuantCONTROL {
    tag "$sample_id1"
    input:
    path(control1)
    path(control2)
    path("salmon_index/*")

    output:
    path "salmon_control"

    script:
    """
    salmon quant --threads 4 --libType=U -i salmon_index -1 $control1 -2 $control2 -o salmon_control
    """
}

process SalmonQuantTEST {
    tag "$sample_id2"
    input:
    path(test1)
    path(test2)
    path("salmon_index/*")

    output:
    path "salmon_test"

    script:
    """
    salmon quant --threads 4 --libType=U -i salmon_index -1 $test1 -2 $test2 -o salmon_test
    """
}

process Hisat2Index {
    input:
    path genome_fa

    output:
    path "hisat2_index.*.ht2"

    script:
    """
    hisat2-build $genome_fa hisat2_index
    """
}

process Hisat2AlignCONTROL {
    tag "$sample_id1"
    input:
    path(control1)
    path(control2)
    path("hisat2_index.*.ht2")

    output:
    path "control_hs2.sam"

    script:
    """
    hisat2 --fast -x hisat2_index -1 $control1 -2 $control2 -S control_hs2.sam
    """
}

process Hisat2AlignTEST {
    tag "$sample_id2"
    input:
    path(test1)
    path(test2)
    path("hisat2_index.*.ht2")

    output:
    path "test_hs2.sam"

    script:
    """
    hisat2 --fast -x hisat2_index -1 $test1 -2 $test2 -S test_hs2.sam
    """
}

process BWAIndex {
    input:
    path genome_fa

    output:
    path "bwa_index.amb"
    path "bwa_index.ann"
    path "bwa_index.bwt"
    path "bwa_index.pac"
    path "bwa_index.sa"

    script:
    """
    bwa index -p bwa_index $genome_fa
    """
}

process BWAAlignCONTROL {
    tag "$sample_id1"
    input:
    path(control1)
    path(control2)
    path "bwa_index.amb"
    path "bwa_index.ann"
    path "bwa_index.bwt"
    path "bwa_index.pac"
    path "bwa_index.sa"

    output:
    path "control_reads.sam"

    script:
    """
    bwa mem bwa_index $control1 $control2 > control_reads.sam
    """
}

process BWAAlignTEST {
    tag "$sample_id2"
    input:
    path(test1)
    path(test2)
    path "bwa_index.amb"
    path "bwa_index.ann"
    path "bwa_index.bwt"
    path "bwa_index.pac"
    path "bwa_index.sa"

    output:
    path "test_reads.sam"

    script:
    """
    bwa mem bwa_index $test1 $test2 > test_reads.sam
    """
}

process SamToBamCONTROL {
    tag "$sample_id1"
    input:
    path("control.sam")

    output:
    path "control.bam"

    script:
    """
    samtools view -@ 4 -b control.sam > control.bam
    """
}

process SamToBamTEST {
    tag "$sample_id2"
    input:
    path("test.sam")

    output:
    path "test.bam"

    script:
    """
    samtools view -@ 4 -b test.sam > test.bam
    """
}

process BamSortedCONTROL {
    tag "$sample_id1"
    input:
    path("control.bam")

    output:
    path "control.sorted.bam"

    script:
    """
    samtools sort -@ 4 control.bam -o control.sorted.bam
    """
}

process BamSortedTEST {
    tag "$sample_id2"
    input:
    path("test.bam")

    output:
    path "test.sorted.bam"

    script:
    """
    samtools sort -@ 4 test.bam -o test.sorted.bam
    """
}

process IndexCONTROL {
    tag "$sample_id1"
    input:
    path("control.sorted.bam")

    output:
    path "control.sorted.bam.bai"

    script:
    """
    samtools index control.sorted.bam
    """
}

process IndexTEST {
    tag "$sample_id2"
    input:
    path("test.sorted.bam")

    output:
    path "test.sorted.bam.bai"

    script:
    """
    samtools index test.sorted.bam
    """
}

process FeatureCountsCONTROL {
    tag "$sample_id1"
    input:
    path("control.sorted.bam")
    path(genome_gtf)

    output:
    path "control_counts"

    script:
    """
    featureCounts -p --countReadPairs -B -C -T 4 -a $genome_gtf -o control_counts control.sorted.bam
    """
}

process FeatureCountsTEST {
    tag "$sample_id2"
    input:
    path("test.sorted.bam")
    path(genome_gtf)

    output:
    path "test_counts"

    script:
    """
    featureCounts -p --countReadPairs -B -C -T 4 -a $genome_gtf -o test_counts test.sorted.bam
    """
}

process MultiQC {
    input:
    path("*")

    output:
    path "."

    script:
    """
    multiqc .
    """
}

process BamSummary {
    input:
    path("control.sorted.bam")
    path("test.sorted.bam")
    path("control.sorted.bam.bai")
    path("test.sorted.bam.bai")

    output:
    path "rnaseq.npz"

    script:
    """
    multiBamSummary bins --bamfiles control.sorted.bam test.sorted.bam -o rnaseq.npz
    """
}

process PCA {
    input:
    path "rnaseq.npz"

    output:
    path "PCA_rnaseq.png"

    script:
    """
    plotPCA -in rnaseq.npz -o PCA_rnaseq.png
    """
}

process Fingerprint {
    input:
    path("control.sorted.bam")
    path("test.sorted.bam")
    path("control.sorted.bam.bai")
    path("test.sorted.bam.bai")

    output:
    path "fingerprint_rnaseq.png"

    script:
    """
    plotFingerprint -b control.sorted.bam test.sorted.bam --labels Control Test --plotFile fingerprint_rnaseq.png
    """
}

process BamCompare {
    input:
    path("control.sorted.bam")
    path("test.sorted.bam")
    path("control.sorted.bam.bai")
    path("test.sorted.bam.bai")

    output:
    path "differential.bw"

    script:
    """
    bamCompare -b1 control.sorted.bam -b2 test.sorted.bam -o differential.bw -of bigwig
    """
}

workflow {
    params.control1 = "$projectDir/../data/control1_1.fq"
    params.control2 = "$projectDir/../data/control1_2.fq"
    params.test1 = "$projectDir/../data/test1_1.fq"
    params.test2 = "$projectDir/../data/test1_2.fq"
    params.genome_fa = "$projectDir/../data/human.fa"
    params.genome_gtf = "$projectDir/../data/hg19_anno.gtf"
    params.gtf = "$projectDir/../data/hg19.gtf"
    params.sample_id1 = "control"
    params.sample_id2 = "test"
    params.work_dir = "$projectDir/../data/"
    
    Channel
        .fromPath(params.control1)
        .set { control1_ch }

    Channel
        .fromPath(params.control2)
        .set { control2_ch }

    Channel
        .fromPath(params.test1)
        .set { test1_ch }

    Channel
        .fromPath(params.test2)
        .set { test2_ch }
    
    Channel
        .fromPath(params.genome_fa)
        .set { fasta_ch }

    Channel
        .fromPath(params.genome_gtf)
        .set { gtf_ch }

    Channel
        .fromPath(params.gtf)
        .set { gtf2_ch }


    // Build Bowtie2 index
    bowtie2_index = Bowtie2Index(fasta_ch)

    // Map control and test samples using Bowtie2
    bowtie2_map_control = Bowtie2MapCONTROL(control1_ch, control2_ch, bowtie2_index)
    bowtie2_map_test = Bowtie2MapTEST(test1_ch, test2_ch, bowtie2_index)

    // Build STAR index
    star_index = STARIndex(fasta_ch, gtf_ch)

    // Map control and test samples using STAR
    star_map_control = STARMapCONTROL(control1_ch, control2_ch, star_index)
    star_map_test = STARMapTEST(test1_ch, test2_ch, star_index)

    // Convert SAM to FASTA for control and test samples
    sam_to_fasta_control = SamToFastaCONTROL(bowtie2_map_control)
    sam_to_fasta_test = SamToFastaTEST(bowtie2_map_test)

    // Build Kallisto index for control and test samples
    kallisto_index_control = KallistoIndexCONTROL(sam_to_fasta_control)
    kallisto_index_test = KallistoIndexTEST(sam_to_fasta_test)

    // Quantify control and test samples using Kallisto
    kallisto_quant_control = KallistoQuantCONTROL(control1_ch, control2_ch, kallisto_index_control)
    kallisto_quant_test = KallistoQuantTEST(test1_ch, test2_ch, kallisto_index_test)

    // Build Salmon index
    salmon_index = SalmonIndex(fasta_ch)

    // Quantify control and test samples using Salmon
    salmon_quant_control = SalmonQuantCONTROL(control1_ch, control2_ch, salmon_index)
    salmon_quant_test = SalmonQuantTEST(test1_ch, test2_ch, salmon_index)

    // Build Hisat2 index
    hisat2 = Hisat2Index(fasta_ch)

    // Align control and test samples using Hisat2
    hisat2_align_control = Hisat2AlignCONTROL(control1_ch, control2_ch, hisat2)
    hisat2_align_test = Hisat2AlignTEST(test1_ch, test2_ch, hisat2)

    // Build BWA index
    bwa_index = BWAIndex(fasta_ch)

    // Align control and test samples using BWA
    bwa_align_control = BWAAlignCONTROL(control1_ch, control2_ch, bwa_index)
    bwa_align_test = BWAAlignTEST(test1_ch, test2_ch, bwa_index)

    // Sort and index SAM files using Samtools for control and test samples
    samtools_view_control = SamToBamCONTROL(bowtie2_map_control)
    samtools_view_test = SamToBamTEST(bowtie2_map_test)
    samtools_sort_control = BamSortedCONTROL(samtools_view_control)
    samtools_sort_test = BamSortedTEST(samtools_view_test)
    samtools_index_control = IndexCONTROL(samtools_sort_control)
    samtools_index_test = IndexTEST(samtools_sort_test)

    // Count features using FeatureCounts for control and test samples
    featurecounts_control = FeatureCountsCONTROL(samtools_sort_control, gtf_ch)
    featurecounts_test = FeatureCountsTEST(samtools_sort_test, gtf_ch)

    // Run MultiQC to summarize the results
    multiqc = MultiQC(bowtie2_map_control)

    // Generate BAM summary
    bam_summary = BamSummary(samtools_sort_control, samtools_sort_test, samtools_index_control, samtools_index_test)

    // Generate PCA plot
    pca = PCA(bam_summary)

    // Generate fingerprint plot
    fingerprint = Fingerprint(samtools_sort_control, samtools_sort_test, samtools_index_control, samtools_index_test)

    // Compare BAM files
    bam_compare = BamCompare(samtools_sort_control, samtools_sort_test, samtools_index_control, samtools_index_test)
}