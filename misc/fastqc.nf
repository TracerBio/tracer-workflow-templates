#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Define the QC process
 */

process QC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastqc=0.12.1"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"), emit: zip
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[reads, "${prefix}.${reads.extension}"]] : reads.withIndex().collect { entry, index -> [entry, "${prefix}_${index + 1}.${entry.extension}"] }
    def rename_to = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect { old_name, new_name -> new_name }.join(' ')
    """
    printf "%s %s\\n" $rename_to | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done

    tracer log "Initiate FASTQC..."
    tracer tool fastqc 0.12.1
    fastqc \\
        $args \\
        --threads $task.cpus \\
        $renamed_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    touch ${prefix}.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}

/*
 * Define the workflow
 */

workflow {
    params.reads = file(params.reads ?: 'data/*.fastq.gz')
    params.meta = params.meta ?: { id -> [id: id] }
    
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { meta, reads -> tuple(params.meta(meta), reads) }
        .set { read_files }

    read_files
        .flatMap { meta, reads -> reads.collect { [meta, it] } }
        .set { individual_reads }

    process_QC(individual_reads)

    read_files
        .flatMap { meta, reads -> [meta, reads] }
        .set { qc_results }

    qc_results.html.view { "QC report: ${it}" }
    qc_results.zip.view { "QC archive: ${it}" }
}

/*
 * Define helper function
 */

def process_QC(input_channel) {
    input_channel
        .into { qc_input }
    
    QC(qc_input)
}