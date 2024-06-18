#!/bin/bash -ue
STAR \
    --genomeDir index \
    --readFilesIn RARA.fq  \
    --runThreadN 1 \
    --outFileNamePrefix Human. \
    --outSAMtype BAM Unsorted \
    --outSAMattrRGline ID:Human 'CN:CRUK' 'SM:Human' \

mkdir -p /workspace/utility-pod/examples/data/results
cp -r ./* /workspace/utility-pod/examples/data/results
samtools sort /workspace/utility-pod/examples/data/results/Human.Aligned.out.bam -@ 4 -o /workspace/utility-pod/examples/data/results/Human.sorted.bam
samtools index /workspace/utility-pod/examples/data/results/Human.sorted.bam

cat <<-END_VERSIONS > versions.yml
"STAR_ALIGN":
    star: $(STAR --version | sed -e "s/STAR_//g")
END_VERSIONS
