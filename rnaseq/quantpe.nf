nextflow.enable.dsl=2

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
    salmon quant --threads 4 --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    tracer end
    """
}