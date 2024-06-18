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
    tracer tool salmon-index 1.1.0
    salmon index --threads 4 -t $transcriptome -i salmon_index
    """
}