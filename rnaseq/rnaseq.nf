nextflow.enable.dsl=2

params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome = "$projectDir/data/ggal/ref1.fa"

process INDEX {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    tracer start
    tracer tool index v1
    salmon index --threads ${task.cpus} -t ${transcriptome} -i salmon_index
    """
}

/*
 * Define the QUANTIFICATION process
 */
process QUANTIFICATION {
    input:
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    tracer tool quant v2
    salmon quant --threads ${task.cpus} --libType=U -i ${salmon_index} -1 ${reads[0]} -2 ${reads[1]} -o ${sample_id}
    tracer end
    """
}

workflow {
    // Create a channel from the reads
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // Run the INDEX process
    index_ch = INDEX(params.transcriptome)

    // Run the QUANTIFICATION process
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
}