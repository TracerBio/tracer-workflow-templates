nextflow.enable.dsl=2

params.reads = "$projectDir/data/human/aml_{1,2}.fq"
params.transcriptome = "$projectDir/data/human/human.fa"

// Include the INDEX process
include { INDEX } from './misc/index.nf'

// Include the QUANTIFICATION process
include { QUANTIFICATION } from './misc/quantpe.nf'

workflow {
    // Create a channel from the reads
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // Run the INDEX process
    index_ch = INDEX(params.transcriptome)

    // Run the QUANTIFICATION process
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
}