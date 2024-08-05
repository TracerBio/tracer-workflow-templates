nextflow.enable.dsl=2

params.index = "$projectDir/../data/hg19"
params.fasta = "$projectDir/../data/human.fa"
params.gtf = "$projectDir/../data/hg19.refGene.gtf"
params.reads = "$projectDir/../data/RARA.fq"
params.bam = "$projectDir/../data/responses/*.bam" 
params.star_ignore_sjdbgtf = false
params.seq_platform = "Illumina"
params.seq_center = "CRUK"

// Meta information
params.meta1 = [id: "Human"]
params.meta2 = [description: "Neuroblastoma"]
params.meta3 = [antibody: "RARA"]

// Include the STAR_GENOMEGENERATE process
include { STAR_GENOMEGENERATE } from './genomeGenerate.nf'

// Include the STAR_ALIGN process
include { STAR_ALIGN } from './star-align.nf'

workflow {
    // Generate genome index
    index_ch = STAR_GENOMEGENERATE(params.meta1, params.fasta, params.meta2, params.gtf)

    // Align reads using the generated index
    align_ch = STAR_ALIGN(
        tuple(params.meta1, file(params.reads)),
        index_ch.index,
        tuple(params.meta2, file(params.gtf)),
        params.star_ignore_sjdbgtf,
        params.seq_platform,
        params.seq_center
    )
}