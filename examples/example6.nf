nextflow.enable.dsl=2

params.fasta = "$projectDir/data/human/human.fa"
params.gtf = "$projectDir/data/human/hg19.refGene.gtf"
params.meta1 = "Human"
params.meta2 = "Neuroblastoma"

// Include the STAR_GENOMEGENERATE process
include { STAR_GENOMEGENERATE } from '../misc/genomeGenerate.nf'
$projectDir/data/human/human.fa
workflow {
    
    index_ch = STAR_GENOMEGENERATE(params.meta1, params.fasta, params.meta2, params.gtf)

    index_ch.view { "Generated index files: $it" }
    index_ch.versions.view { "Tool versions: $it" }
}