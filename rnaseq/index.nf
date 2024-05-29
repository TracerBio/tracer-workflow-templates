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
    tracer tool index v1
    salmon index --threads 4 -t $transcriptome -i salmon_index
    """
}