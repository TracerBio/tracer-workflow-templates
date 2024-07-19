nextflow.enable.dsl=2

process STAR_INDEX {
    input:
    path 'human.fa'
    path 'hg19.refGene.gtf'
    
    output:
    path 'human/*'
    
    script:
    """
    STAR --runThreadN 4 --runMode genomeGenerate --genomeDir human --genomeSAindexNbases 10 --genomeFastaFiles human.fa --sjdbGTFfile hg19.refGene.gtf --sjdbOverhang 99
    """
}

process STAR_ALIGN {
    tag "$sample_id"
    input:
    path 'human/*'
    path sample1
    path sample2
    
    output:
    path 'P1s1*.sorted.bam'
    path 'P1s1*.sorted.bam.bai'

    script:
    """
    STAR --runThreadN 4 --genomeDir human --readFilesIn $sample1 $sample2 --outFileNamePrefix P1s1 --outSAMtype BAM SortedByCoordinate
    samtools sort P1s1*.bam -@ 4 -o P1s1.sorted.bam
    samtools index P1s1.sorted.bam
    """
}

process MACS {
    tag "$sample_id"
    input:
    path 'P1s1.sorted.bam'
    path 'P1s1.sorted.bam.bai'
    
    output:
    path './P1s1_peaks/'
    
    script:
    """
    macs3 callpeak -t P1s1.sorted.bam -f BAMPE -p 0.05 --outdir ./P1s1_peaks
    """
}

process DEEPTOOLS {
    tag "$sample_id"
    input:
    path 'P1s1.sorted.bam'
    path 'P1s1.sorted.bam.bai'
    
    output:
    path 'coverage.pdf'
    
    script:
    """
    plotCoverage -b P1s1.sorted.bam -o coverage.pdf
    """
}

workflow {
    params.sample_id = "P1s1"
    params.sample1 = "$projectDir/../data/s1_1.fq"
    params.sample2 = "$projectDir/../data/s1_2.fq"
    params.genome_fasta = "$projectDir/../data/human.fa"
    params.annotation_gtf = "$projectDir/../data/hg19.refGene.gtf"
    
    Channel
        .fromPath(params.sample1)
        .set { sample1_ch }

    Channel
        .fromPath(params.sample2)
        .set { sample2_ch }

    Channel
        .fromPath(params.genome_fasta)
        .set { genome_fasta_ch }

    Channel
        .fromPath(params.annotation_gtf)
        .set { annotation_gtf_ch }

    star_index_output = STAR_INDEX(genome_fasta_ch, annotation_gtf_ch)
    
    star_align_output = STAR_ALIGN(star_index_output, sample1_ch, sample2_ch)
    
    macs_output = MACS(star_align_output)
    
    deeptools_output = DEEPTOOLS(star_align_output)
}
