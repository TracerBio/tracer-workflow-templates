nextflow.enable.dsl=2

params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome = "$projectDir/data/ggal/ref1.fa"

// Include the INDEX process
include { INDEX } from './index.nf'

// Include the QUANTIFICATION process
include { QUANTIFICATION } from './quantpe.nf'

workflow {
    // Create a channel from the reads
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // Run the INDEX process
    index_ch = INDEX(params.transcriptome)

    // Run the QUANTIFICATION process
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
}