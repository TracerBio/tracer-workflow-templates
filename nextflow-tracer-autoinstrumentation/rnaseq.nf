nextflow.enable.dsl=2

params.control1 = "$projectDir/../data/control1_1.fq"
params.control2 = "$projectDir/../data/control1_2.fq"
params.test1 = "$projectDir/../data/test1_1.fq"
params.test2 = "$projectDir/../data/test1_2.fq"
params.test2 = "$projectDir/../data/test1_2.fq"
params.output_dir = "$projectDir/../data/s1_1.fq"
params.genome_fa = "$projectDir/../data/human.fa"
params.genome_gtf = "$projectDir/../data/hg19_anno.gtf"
params.sample_id1 = "control"
params.sample_id2 = "test"

process FastQC-CONTROL {
    tag "$sample_id1"
    input:
    path(control1)
    path(control2)

    output:
    path '$projectDir/../data/'

    script:
    """
    fastqc $control1 $control2 -o .
    """
}

process FastQC-TEST {
    tag "$sample_id2"
    input:
    path(test1)
    path(test2)

    output:
    path '$projectDir/../data/'

    script:
    """
    fastqc $test1 $test2 -o .
    """
}

process Bowtie2Index {
    input:
    path genome_fa

    output:
    path "$projectDir/../data/human_index"

    script:
    """
    bowtie2-build $genome_fa human_index
    """
}

process Bowtie2Map-CONTROL {
    tag "$sample_id1"
    input:
    tuple val(sample_id1), path(control1), path(control2), path 'human_index'

    output:
    path '$projectDir/../data/control.sam'

    script:
    """
    bowtie2 --local -x human_index -1 $control1 -2 $control2 -S control.sam
    """
}

process Bowtie2Map-TEST {
    tag "$sample_id2"
    input:
    tuple val(sample_id2), path(test1), path(test2), path 'human_index'

    output:
    path '$projectDir/../data/test.sam'

    script:
    """
    bowtie2 --local -x human_index -1 $test1 -2 $test2 -S test.sam
    """
}

process STARIndex {
    input:
    path genome_fa, path genome_gtf

    output:
    path '$projectDir/../data/human_star'

    script:
    """
    STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./human_star --genomeSAindexNbases 10 --genomeFastaFiles $genome_fa --sjdbGTFfile $genome_gtf --sjdbOverhang 99
    """
}

process STARMap-CONTROL {
    tag "$sample_id1"
    input:
    tuple val(sample_id1), path(control1), path(control2), path 'human_star'

    output:
    path "$projectDir/../data/control_star_Aligned.sortedByCoord.out.bam"
    path "$projectDir/../data/control_star_ReadsPerGene.out.tab"

    script:
    """
    STAR --runThreadN 4 --genomeDir ./human_star --readFilesIn $control1 $control2 --outFileNamePrefix control_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
    """
}

process STARMap-TEST {
    tag "$sample_id2"
    input:
    tuple val(sample_id2), path(test1), path(test2), path 'human_star'

    output:
    path "$projectDir/../data/test_star_Aligned.sortedByCoord.out.bam"
    path "$projectDir/../data/test_star_ReadsPerGene.out.tab"

    script:
    """
    STAR --runThreadN 4 --genomeDir ./human_star --readFilesIn $test1 $test2 --outFileNamePrefix test_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts
    """
}

process SamToFasta-CONTROL {
    tag "$sample_id1"
    input:
    path 'control.sam'

    output:
    path "$projectDir/../data/control.fa"

    script:
    """
    samtools fasta control.sam > control.fa
    """
}

process SamToFasta-TEST {
    tag "$sample_id2"
    input:
    path 'test.sam'

    output:
    path "$projectDir/../data/test.fa"

    script:
    """
    samtools fasta test.sam > test.fa
    """
}

process KallistoIndex-CONTROL {
    tag "$sample_id1"
    input:
    tuple val(sample_id1), path 'control.fa'

    output:
    path "$projectDir/../data/control_index"

    script:
    """
    kallisto index -t 4 -i control_index control.fa
    """
}

process KallistoIndex-TEST {
    tag "$sample_id2"
    input:
    tuple val(sample_id2), path 'test.fa'

    output:
    path "$projectDir/../data/test_index"

    script:
    """
    kallisto index -t 4 -i test_index test.fa
    """
}

process KallistoQuant-CONTROL {
    tag "$sample_id1"
    input:
    tuple val(sample_id1), path(control1), path(control2), path 'control_index'

    output:
    path "$projectDir/../data/control_quant"

    script:
    """
    kallisto quant -t 4 -i control_index -o ./control_quant $control1 $control2
    """
}

process KallistoQuant-TEST {
    tag "$sample_id2"
    input:
    tuple val(sample_id2), path(test1), path(test2), path 'test_index'

    output:
    path "$projectDir/../data/test_quant"

    script:
    """
    kallisto quant -t 4 -i test_index -o ./test_quant $test1 $test2
    """
}

process SalmonIndex {
    input:
    path genome_fa

    output:
    path "$projectDir/../data/salmon_index"

    script:
    """
    salmon index --threads 4 -t $genome_fa -i salmon_index
    """
}

process SalmonQuant-CONTROL {
    tag "$sample_id1"
    input:
    tuple val(sample_id1), path(control1), path(control2), path 'salmon_index'

    output:
    path "$projectDir/../data/salmon_control"

    script:
    """
    salmon quant --threads 4 --libType=U -i salmon_index -1 $control1 -2 $control2 -o ./salmon_control
    """
}

process SalmonQuant-TEST {
    tag "$sample_id2"
    input:
    tuple val(sample_id2), path(test1), path(test2), path 'salmon_index'

    output:
    path "$projectDir/../data/salmon_test"

    script:
    """
    salmon quant --threads 4 --libType=U -i salmon_index -1 $test1 -2 $test2 -o ./salmon_test
    """
}

process Hisat2Index {
    input:
    path genome_fa

    output:
    path "$projectDir/../data/hisat2_index"

    script:
    """
    hisat2-build $genome_fa hisat2_index
    """
}

process Hisat2Align-CONTROL {
    tag "$sample_id1"
    input:
    tuple val(sample_id1), path(control1), path(control2), path 'hisat2_index'

    output:
    path "$projectDir/../data/control_hs2.sam"

    script:
    """
    hisat2 --fast -x hisat2_index -1 $control1 -2 $control2 -S control_hs2.sam
    """
}

process Hisat2Align-TEST {
    tag "$sample_id2"
    input:
    tuple val(sample_id2), path(test1), path(test2), path 'hisat2_index'

    output:
    path "$projectDir/../data/test_hs2.sam"

    script:
    """
    hisat2 --fast -x hisat2_index -1 $test1 -2 $test2 -S test_hs2.sam
    """
}

process BWAIndex {
    input:
    path genome_fa

    output:
    path "$projectDir/../data/bwa_index"

    script:
    """
    bwa index -p bwa_index $genome_fa
    """
}

process BWAAlign-CONTROL {
    tag "$sample_id1"
    input:
    tuple val(sample_id1), path(control1), path(control2), path 'bwa_index'

    output:
    path "$projectDir/../data/control_reads.sam"

    script:
    """
    bwa mem bwa_index $control1 $control2 > control_reads.sam
    """
}

process BWAAlign-TEST {
    tag "$sample_id2"
    input:
    tuple val(sample_id2), path(test1), path(test2), path 'bwa_index'

    output:
    path "$projectDir/../data/test_reads.sam"

    script:
    """
    bwa mem bwa_index $test1 $test2 > test_reads.sam
    """
}

process SamtoolsSortIndex-CONTROL {
    tag "$sample_id1"
    input:
    path("control.sam")

    output:
    path("$projectDir/../data/control.bam")
    path("$projectDir/../data/control.sorted.bam")

    script:
    """
    samtools view -@ 4 -b control.sam > control.bam
    samtools sort -@ 4 control.bam -o control.sorted.bam
    samtools index control.sorted.bam
    """
}

process SamtoolsSortIndex-TEST {
    tag "$sample_id2"
    input:
    path("test.sam")

    output:
    path("$projectDir/../data/test.bam")
    path("$projectDir/../data/test.sorted.bam")

    script:
    """
    samtools view -@ 4 -b test.sam > test.bam
    samtools sort -@ 4 test.bam -o test.sorted.bam
    samtools index test.sorted.bam
    """
}

process StringTie-CONTROL {
    tag "$sample_id1"
    input:
    path("$projectDir/../data/control.sorted.bam"), path(genome_gtf)

    output:
    path("$projectDir/../data/control.gtf")

    script:
    """
    stringtie -o control.gtf -G $genome_gtf control.sorted.bam
    """
}

process StringTie-TEST {
    tag "$sample_id2"
    input:
    path("$projectDir/../data/test.sorted.bam"), path(genome_gtf)

    output:
    path("$projectDir/../data/test.gtf")

    script:
    """
    stringtie -o test.gtf -G $genome_gtf test.sorted.bam
    """
}

process FeatureCounts-CONTROL {
    tag "$sample_id1"
    input:
    path("control.sorted.bam"), path(genome_gtf)

    output:
    path("$projectDir/../data/control_counts")

    script:
    """
    featureCounts -p --countReadPairs -B -C -T 4 -a $genome_gtf -o control_counts control.sorted.bam
    """
}

process FeatureCounts-TEST {
    tag "$sample_id2"
    input:
    path("test.sorted.bam"), path(genome_gtf)

    output:
    path("$projectDir/../data/test_counts")

    script:
    """
    featureCounts -p --countReadPairs -B -C -T 4 -a $genome_gtf -o test_counts test.sorted.bam
    """
}

process MultiQC {
    input:
    path("$projectDir/../data/*")

    output:
    path("$projectDir/../data/")

    script:
    """
    multiqc .
    """
}

process BamSummary {
    input:
    path("control.sorted.bam"), path("test.sorted.bam")

    output:
    path "$projectDir/../data/rnaseq.npz"

    script:
    """
    multiBamSummary bins --bamfiles control.sorted.bam test.sorted.bam -o rnaseq.npz
    """
}

process PCA {
    input:
    path "$projectDir/../data/rnaseq.npz"

    output:
    path "$projectDir/../data/PCA_rnaseq.png"

    script:
    """
    plotPCA -in rnaseq.npz -o PCA_rnaseq.png
    """
}

process Fingerprint {
    input:
    path("control.sorted.bam"), path("test.sorted.bam")

    output:
    path "$projectDir/../data/fingerprint_rnaseq.png"

    script:
    """
    plotFingerprint -b ${bam_files.join(' ')} --labels Control Test --plotFile fingerprint_rnaseq.png
    """
}

process BamCompare {
    input:
    tuple path(test_bam), path(control_bam)

    output:
    path "$projectDir/../data/differential.bw"

    script:
    """
    bamCompare -b1 control.sorted.bam -b2 test.sorted.bam -o differential.bw -of bigwig
    """
}

workflow {
    // Channel to read input files
    read_pairs = Channel.fromFilePairs("${params.reads_dir}/*_{1,2}.fq", flat: false)

    // FASTQC process
    fastqc_results = read_pairs.map { sample_id, reads -> tuple(sample_id, reads) }
                              .set { read_pairs_ch }
                              .map { sample_id, reads -> tuple(sample_id, reads) }
                              .into { fastqc_ch; fastqc_ch2 }

    // Bowtie2 index
    bowtie2_index = Bowtie2Index(genome_fa: file(params.genome_fa))

    // Bowtie2 mapping
    bowtie2_map = read_pairs_ch.map { sample_id, reads -> tuple(sample_id, reads, bowtie2_index) }
                               .set { bowtie2_map_ch }
    Bowtie2Map(bowtie2_map_ch)

    // STAR index
    star_index = STARIndex(genome_fa: file(params.genome_fa), genome_gtf: file(params.genome_gtf))

    // STAR mapping
    star_map = read_pairs_ch.map { sample_id, reads -> tuple(sample_id, reads, star_index) }
                            .set { star_map_ch }
    STARMap(star_map_ch)

    // Sam to Fasta
    sam_files = Channel.of(["control.sam", "test.sam"])
    SamToFasta(sam_files)

    // Kallisto index
    fasta_files = Channel.of(["control.fa", "test.fa"])
    kallisto_index = fasta_files.map { fasta -> tuple(fasta, fasta) }
    KallistoIndex(kallisto_index)

    // Kallisto quantification
    kallisto_quant = read_pairs_ch.map { sample_id, reads -> tuple(sample_id, reads, kallisto_index) }
                                  .set { kallisto_quant_ch }
    KallistoQuant(kallisto_quant_ch)

    // Salmon index
    salmon_index = SalmonIndex(genome_fa: file(params.genome_fa))

    // Salmon quantification
    salmon_quant = read_pairs_ch.map { sample_id, reads -> tuple(sample_id, reads, salmon_index) }
                                .set { salmon_quant_ch }
    SalmonQuant(salmon_quant_ch)

    // Hisat2 index
    hisat2_index = Hisat2Index(genome_fa: file(params.genome_fa))

    // Hisat2 alignment
    hisat2_align = read_pairs_ch.map { sample_id, reads -> tuple(sample_id, reads, hisat2_index) }
                                .set { hisat2_align_ch }
    Hisat2Align(hisat2_align_ch)

    // BWA index
    bwa_index = BWAIndex(genome_fa: file(params.genome_fa))

    // BWA alignment
    bwa_align = read_pairs_ch.map { sample_id, reads -> tuple(sample_id, reads, bwa_index) }
                             .set { bwa_align_ch }
    BWAAlign(bwa_align_ch)

    // Samtools sort and index
    sam_files_from_alignments = Channel.of(["control.sam", "test.sam"])
    SamtoolsSortIndex(sam_files_from_alignments)

    // StringTie
    bam_files = Channel.of(["control.sorted.bam", "test.sorted.bam"])
    StringTie(bam_files)

    // FeatureCounts
    feature_counts = Channel.of(["control.sorted.bam", "test.sorted.bam"])
    FeatureCounts(feature_counts)

    // MultiQC
    multiqc_results = Channel.of("${params.output_dir}/*")
    MultiQC(multiqc_results)

    // BamSummary
    bam_files_for_summary = Channel.of(["control.sorted.bam", "test.sorted.bam"])
    BamSummary(bam_files_for_summary)

    // PCA
    PCA(bam_files_for_summary)

    // Fingerprint
    Fingerprint(bam_files_for_summary)

    // BamCompare
    BamCompare(bam_files_for_summary.collect())
}
