nextflow.enable.dsl=2

params.fasta = "$projectDir/data/human/human.fa"
params.gtf = "$projectDir/data/human/hg19.refGene.gtf"
params.reads = "$projectDir/data/human/RARA.fq" 
params.star_ignore_sjdbgtf = false
params.seq_platform = "Illumina"
params.seq_center = "CRUK"

// Meta information
params.meta1 = [id: "Human"]
params.meta2 = [description: "Neuroblastoma"]

// Include the STAR_GENOMEGENERATE process
include { STAR_GENOMEGENERATE } from '../misc/genomeGenerate.nf'

// Include the STAR_ALIGN process
include { STAR_ALIGN } from '../misc/star-align.nf'

workflow {
    // Generate genome index
    index_ch = STAR_GENOMEGENERATE(tuple(params.meta1, file(params.fasta), params.meta2, file(params.gtf)))

    // Align reads using the generated index
    STAR_ALIGN(
        tuple(params.meta1, file(params.reads)),
        index_ch.index,
        tuple(params.meta2, file(params.gtf)),
        params.star_ignore_sjdbgtf,
        params.seq_platform,
        params.seq_center
    )

    // View the output of alignment
    STAR_ALIGN.out.view { "Alignment output: $it" }
}
