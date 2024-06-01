nextflow.enable.dsl=2

process STAR_GENOMEGENERATE {
    tag "$fasta"
    label 'process_high'

    input:
    tuple val(meta) 
    path(fasta)
    tuple val(meta2) 
    path(gtf)

    output:
    tuple val(meta), path("star")  , emit: index
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args_list = args.tokenize()
    def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
        """
        mkdir star
        tracer start
        tracer tool index 2.7.10
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir ./examples/data/human/hg19 \\
            --genomeFastaFiles $fasta \\
            --sjdbGTFfile $gtf \\
            --runThreadN $task.cpus \\
            $memory \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            star: \$(STAR --version | sed -e "s/STAR_//g")
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
        END_VERSIONS
        """

    stub:
    """
    mkdir star
    touch star/Genome
    touch star/Log.out
    touch star/SA
    touch star/SAindex
    touch star/chrLength.txt
    touch star/chrName.txt
    touch star/chrNameLength.txt
    touch star/chrStart.txt
    touch star/exonGeTrInfo.tab
    touch star/exonInfo.tab
    touch star/geneInfo.tab
    touch star/genomeParameters.txt
    touch star/sjdbInfo.txt
    touch star/sjdbList.fromGTF.out.tab
    touch star/sjdbList.out.tab
    touch star/transcriptInfo.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}