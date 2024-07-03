#!/bin/bash
 
# FASTQC - QC layer 1
cd /workspace/tracer-workflow-templates/data;

echo "Running FASTQC...";
start_fastqc=$(date +%s%3N);
fastqc s1_1.fq s1_2.fq -o .;
end_fastqc=$(date +%s%3N);
fastqc_duration=$((end_fastqc - start_fastqc));

# Creating a human genome index using STAR
echo "Creating a human genome index using STAR...";
start_star_genome=$(date +%s%3N);
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./mouse_genome --genomeSAindexNbases 10 --genomeFastaFiles mouse.fa --sjdbGTFfile mm10.gtf --sjdbOverhang 99;
end_star_genome=$(date +%s%3N);
star_genome_duration=$((end_star_genome - start_star_genome));

# Align RNA-Seq reads to the genome using STAR
echo "Aligning RNA-Seq reads to the genome using STAR...";
start_star_align=$(date +%s%3N);
STAR --runThreadN 4 --genomeDir ./mouse_genome --readFilesIn s1_1.fq s1_2.fq --outFileNamePrefix P1s1_mouse --outSAMtype BAM SortedByCoordinate;
end_star_align=$(date +%s%3N);
star_align_duration=$((end_star_align - start_star_align));

# Sort and Index the output bam file
echo "Sorting and indexing the BAM file...";
start_samtools=$(date +%s%3N);
samtools sort *.bam -@ 4 -o P1s1_mouse.sorted.bam;
samtools index P1s1_mouse.sorted.bam;
end_samtools=$(date +%s%3N);
samtools_duration=$((end_samtools - start_samtools));

# Call Peaks for BAM file
echo "Calling peaks for BAM file...";
start_macs=$(date +%s%3N);
macs3 callpeak -t P1s1_mouse.sorted.bam -f BAMPE -p 0.05 --outdir ./Peaks;
end_macs=$(date +%s%3N);
macs_duration=$((end_macs - start_macs));

# Plot Coverage of Peaks
echo "Plotting coverage of peaks...";
start_plotCoverage=$(date +%s%3N);
plotCoverage -b P1s1_mouse.sorted.bam -o coverage.pdf;
end_plotCoverage=$(date +%s%3N); 
plotCoverage_duration=$((end_plotCoverage - start_plotCoverage));

# Record end time of the entire script
end_time=$(date +%s%3N)
total_duration=$((end_time - start_time))

# Print the durations in milliseconds
{
	echo "Execution duration of FASTQC: $fastqc_duration milliseconds"
	echo "Execution duration of STAR genome generation: $star_genome_duration milliseconds"
	echo "Execution duration of STAR alignment: $star_align_duration milliseconds"
	echo "Execution duration of samtools: $samtools_duration milliseconds"
	echo "Execution duration of MACS peak calling: $macs_duration milliseconds"
	echo "Execution duration of plotCoverage: $plotCoverage_duration milliseconds"
	echo "Total execution duration of the script: $total_duration milliseconds"
} > output.txt