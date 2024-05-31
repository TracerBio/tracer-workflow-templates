nextflow.enable.dsl=2

params.reads = "/workspace/utility-pod/examples/data/ggal/gut_{1,2}.fq"
params.meta = "nextcore"

// Include the INDEX process
include { QC } from '../misc/fastqc.nf'

/*
 * Define the workflow
 */

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { id, reads -> tuple([id: id], reads) }
        .set { read_files }

    read_files
        .flatMap { meta, reads -> reads.collect { [meta, it] } }
        .set { individual_reads }

    qc_ch = QC(individual_reads)
}