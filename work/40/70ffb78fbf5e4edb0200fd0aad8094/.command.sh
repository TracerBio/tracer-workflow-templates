#!/bin/bash -ue
mkdir star
cd star  
STAR \
    --genomeDir index \
    --readFilesIn RARA.fq  \
    --runThreadN 1 \
    --outFileNamePrefix Human. \
    --outSAMtype BAM Unsorted \
    --outSAMattrRGline ID:Human 'CN:CRUK' 'SM:Human' \


samtools sort Human.Aligned.out.bam -@ 4 -o Human.sorted.bam
samtools index Human.sorted.bam

cat <<-END_VERSIONS > versions.yml
"STAR_ALIGN":
    star: $(STAR --version | sed -e "s/STAR_//g")
END_VERSIONS
