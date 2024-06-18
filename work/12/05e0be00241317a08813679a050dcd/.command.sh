#!/bin/bash -ue
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir human --genomeSAindexNbases 10 --genomeFastaFiles human.fa --sjdbGTFfile hg19.refGene.gtf --sjdbOverhang 99
