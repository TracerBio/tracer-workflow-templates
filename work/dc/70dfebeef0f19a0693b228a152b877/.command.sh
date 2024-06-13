#!/bin/bash -ue
tracer tool star-align 1.7.0
STAR --runThreadN 4 --genomeDir human --readFilesIn s1_1.fq s1_2.fq --outFileNamePrefix P1s1 --outSAMtype BAM SortedByCoordinate
tracer tool samtools 1.17
samtools sort P1s1*.bam -@ 4 -o P1s1.sorted.bam
samtools index P1s1.sorted.bam
