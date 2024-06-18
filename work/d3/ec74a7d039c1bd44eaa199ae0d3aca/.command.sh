#!/bin/bash -ue
mkdir star
tracer start
tracer tool index 2.7.10
STAR \
    --runMode genomeGenerate \
    --genomeDir index/ \
    --genomeSAindexNbases 10 \
    --genomeFastaFiles human.fa \
    --sjdbGTFfile hg19.refGene.gtf \
    --runThreadN 1 \
     \


cat <<-END_VERSIONS > versions.yml
"STAR_GENOMEGENERATE":
    star: $(STAR --version | sed -e "s/STAR_//g")
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
    gawk: $(echo $(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*$//')
END_VERSIONS
