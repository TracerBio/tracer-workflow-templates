from airflow import DAG
from airflow.operators.bash import BashOperator
from airflow.utils.dates import days_ago
from airflow.operators.python import PythonOperator
from airflow.hooks.base import BaseHook
import os

# Define the data directory
DATA_DIR = "/workspace/tracer-workflow-templates/data"

# Define the default arguments for the DAG
default_args = {
    'owner': 'airflow',
    'depends_on_past': False,
    'email_on_failure': False,
    'email_on_retry': False,
    'retries': 3,
}

# Initialize the DAG
dag = DAG(
    'rna_seq_pipeline',
    default_args=default_args,
    description='RNA-Seq Analysis Pipeline',
    schedule_interval=None,  # Set as needed
    start_date=days_ago(1),
    catchup=False,
)

# Task: STAR genome index creation
star_index = BashOperator(
    task_id='star_index',
    bash_command=f"""
    STAR --runThreadN 4 --runMode genomeGenerate \
        --genomeDir {DATA_DIR}/human_star --genomeSAindexNbases 10 \
        --genomeFastaFiles {DATA_DIR}/human.fa --sjdbGTFfile {DATA_DIR}/hg19_anno.gtf \
        --sjdbOverhang 99
    """,
    dag=dag,
)

# Task: STAR alignment for control sample
star_align_control = BashOperator(
    task_id='star_align_control',
    bash_command=f"""
    STAR --runThreadN 4 --genomeDir {DATA_DIR}/human_star \
        --readFilesIn {DATA_DIR}/control1_1.fq {DATA_DIR}/control1_2.fq \
        --outFileNamePrefix {DATA_DIR}/control_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts;
    mv {DATA_DIR}/control_starAligned.sortedByCoord.out.bam {DATA_DIR}/control.sorted.bam
    """,
    dag=dag,
)

# Task: STAR alignment for test sample
star_align_test = BashOperator(
    task_id='star_align_test',
    bash_command=f"""
    STAR --runThreadN 4 --genomeDir {DATA_DIR}/human_star \
        --readFilesIn {DATA_DIR}/test1_1.fq {DATA_DIR}/test1_2.fq \
        --outFileNamePrefix {DATA_DIR}/test_star --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts;
    mv {DATA_DIR}/test_starAligned.sortedByCoord.out.bam {DATA_DIR}/test.sorted.bam
    """,
    dag=dag,
)

# Task: Index control BAM files
samtools_index_control = BashOperator(
    task_id='samtools_index_control',
    bash_command=f"samtools index {DATA_DIR}/control.sorted.bam",
    dag=dag,
)

# Task: Index test BAM files
samtools_index_test = BashOperator(
    task_id='samtools_index_test',
    bash_command=f"samtools index {DATA_DIR}/test.sorted.bam",
    dag=dag,
)

# Task: FeatureCounts for control
featurecounts_control = BashOperator(
    task_id='featurecounts_control',
    bash_command=f"""
    featureCounts -p --countReadPairs -B -C -T 4 -a {DATA_DIR}/hg19_anno.gtf \
        -o {DATA_DIR}/control_counts.txt {DATA_DIR}/control.sorted.bam
    """,
    dag=dag,
)

# Task: FeatureCounts for test
featurecounts_test = BashOperator(
    task_id='featurecounts_test',
    bash_command=f"""
    featureCounts -p --countReadPairs -B -C -T 4 -a {DATA_DIR}/hg19_anno.gtf \
        -o {DATA_DIR}/test_counts.txt {DATA_DIR}/test.sorted.bam
    """,
    dag=dag,
)

# Task: MultiQC
multiqc = BashOperator(
    task_id='multiqc',
    bash_command=f"multiqc {DATA_DIR} -o {DATA_DIR}/multiqc_report",
    dag=dag,
)

# Task: Deeptools summary
deeptools_summary = BashOperator(
    task_id='deeptools_summary',
    bash_command=f"""
    multiBamSummary bins --bamfiles {DATA_DIR}/control.sorted.bam {DATA_DIR}/test.sorted.bam -o {DATA_DIR}/rnaseq.npz
    """,
    dag=dag,
)

# Task: PCA plot
pca_plot = BashOperator(
    task_id='pca_plot',
    bash_command=f"plotPCA -in {DATA_DIR}/rnaseq.npz -o {DATA_DIR}/PCA_rnaseq.png",
    dag=dag,
)

# Task: Fingerprint plot
fingerprint_plot = BashOperator(
    task_id='fingerprint_plot',
    bash_command=f"""
    plotFingerprint -b {DATA_DIR}/control.sorted.bam {DATA_DIR}/test.sorted.bam --labels Control Test --plotFile {DATA_DIR}/fingerprint_rnaseq.png
    """,
    dag=dag,
)

# Task: BamCompare
bamcompare = BashOperator(
    task_id='bamcompare',
    bash_command=f"""
    bamCompare -b1 {DATA_DIR}/test.sorted.bam -b2 {DATA_DIR}/control.sorted.bam -o {DATA_DIR}/differential.bw -of bigwig
    """,
    dag=dag,
)

# Task: Run DESeq2 R script
run_deseq2_r = BashOperator(
    task_id='run_deseq2_r',
    bash_command=f"Rscript {DATA_DIR}/DESeq2.r",
    dag=dag,
)

# Set task dependencies
star_index >> [star_align_control, star_align_test]
star_align_control >> samtools_index_control
star_align_test >> samtools_index_test
samtools_index_control >> featurecounts_control
samtools_index_test >> featurecounts_test
[featurecounts_control, featurecounts_test] >> multiqc
multiqc >> deeptools_summary
deeptools_summary >> pca_plot >> fingerprint_plot >> bamcompare >> run_deseq2_r