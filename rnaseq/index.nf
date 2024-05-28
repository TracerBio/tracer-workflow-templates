nextflow.enable.dsl=2

/*
 * Define the INDEX process
 */
process INDEX {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    tracer start
    tracer tool index v1
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}