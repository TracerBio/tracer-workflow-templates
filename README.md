# Tracer Workflow Templates

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/TracerBio/tracer-workflow-templates)

This repository provides varied examples of performing RNASEQ & CHIPSEQ analyses using bash scripts and pipeline-based systems (Nextflow, Snakemake, Airflow),  integrated with our tracking & observability tool - Tracer.


# Example 1: bash CHIPSEQ training example with indexing, mapping & quantification

1. "cd shell-tracer-autoinstrumentation"
2. "sh chipseq.sh"
3. Tools included: STAR, samtools, MACS3, plotCoverage [deeptools]

# Example 2: bash RNASEQ training example with indexing, mapping & quantification

1. "cd shell-tracer-autoinstrumentation"
2. "sh rnaseq.sh"
3. Tools included: bowtie2, STAR, samtools, kallisto, salmon, hisat2, bwa, stringTie, featureCounts, MultiQC, multiBamSummary, plotPCA, plotFingerprint, bamCompare

# Example 3: bash + RScript RNASEQ training example with indexing, mapping, quantification and data visualization

1. "cd shell-tracer-autoinstrumentation"
2. "sh rnaseqv3-0.sh"
3. Tools included: bowtie2, STAR, samtools, kallisto, salmon, hisat2, bwa, stringTie, featureCounts, MultiQC, multiBamSummary, plotPCA, plotFingerprint, bamCompare
4. R libraries included: edgeR, ggplot2, reshape2

# Example 4: Snakemake RNASEQ training examples - indexing, mapping, quantification and data visualization  

1. Run "sh ./config/airflow.sh" to configure Airflow scheduler and webserver at port 8080
2. Enter Password: tracertothestars
3. Open webserver portal with credentials:
    - Username: tracer
    - Password: tracertothestars
4. Trigger run for rnaseq_dag (RNASEQ workflow) 
4. Tools included: bowtie2, STAR, samtools, kallisto, salmon, hisat2, bwa, stringTie, featureCounts, MultiQC, multiBamSummary, plotPCA, plotFingerprint, bamCompare
5. R libraries included: edgeR, ggplot2, reshape2

# Example 5: Airflow RNASEQ training example - indexing, mapping, quantification and data visualization 

1. Generates a human genome index and maps the CHIPSEQ data (.fq) to the index
2. Details of dataset available at "./examples/data/human/datasets.txt"
2. "nextflow run ./nextflow-tracer-inline-code/example5.nf"

# Example 6: Nextflow complete CHIPSEQ workflow - from indexing to calling peaks

1. "cd nextflow-tracer-autoinstrumentation"
2. "nextflow run chipseq.nf"
3. Tools included: STAR, samtools, MACS3, plotCoverage [deeptools]

# Example 7: Nextflow complete RNASEQ workflow - from indexing to gene counts

1. "cd nextflow-tracer-autoinstrumentation"
2. "nextflow run rnaseq.nf"
3. Tools included: bowtie2, STAR, samtools, kallisto, salmon, hisat2, bwa, stringTie, featureCounts, MultiQC, multiBamSummary, plotPCA, plotFingerprint, bamCompare

    ## Files

    - `chipseq.sh`, `chipseq.nf`: The main CHIPSEQ bash and Nextflow workflow script files.
    - `rnaseq.sh`, `rnaseq.nf`: The main RNASEQ bash and Nextflow workflow script files.
    - `rnaseq2-0.sh`: This script combines bash and R to perform RNASEQ analysis - from indexing to data visualization.
    - `rnaseq3-0.sh`: This script combines bash and RScript to perform RNASEQ analysis - from indexing to data visualization.
    - `Snakefile`: Default snakemake RNASEQ pipeline file.
    - `rnaseq.smk`: RNASEQ snakemake pipeline file.
    - `rnaseq_dag.py`: RNASEQ airflow pipeline file.   

    ## Workspace Configuration

    - Conda 
    - r-base
    - Mapping tools: salmon, bowtie2, STAR, hisat2, bwa, kallisto    
    - deeptools, samtools
    - MACS3, stringTie, featureCounts
    - Nextflow
    - Snakemake
    - Airflow
    - JAVA

    ## Archive

    Contains other example shell & Nextflow scripts for CHIPSEQ, RNASEQ.    

    ## Usage

    1. To clone the repository:

    git clone https://github.com/TracerBio/tracer-workflow-templates

    cd tracer-workflow-templates
    

    ## Tracer Integration

    Tracer is run as a daemon in your bioinformatics environment.

    Use `tracer start` to initiate a new trackable run & record every step of your bioinformatics workflow!


    ## License

    This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

    ## Acknowledgments

    - Thanks to the Nextflow community for their excellent documentation and support.

    ## GitPod

    [![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/TracerBio/tracer-workflow-templates)
---


# Tracer-daemon instructions
**how to install tracer**
To install Tracer on your initialized Gitpod workspace:

--> Create your account on Tracer 
--> Copy & Paste one-line installation script at https://app.tracer.bio/app/setup

**how to run tracer daemon** 
```bash
tracer-deamon 
# should return similar to: 1305 ?        00:00:00 tracer-daemon
```
**How to check if Tracer Daemon Is Running**
```bash
ps -e | grep tracer

```

**latest stable version June 21st**
```bash
git checkout cbd49684df61419a08e41f1a3ebe99507a2de0ed
```

This README provides a comprehensive guide for users to understand and execute a nextflow RNA-Seq workflow with Tracer enabled tracking. Adjust the repository URL, license details, and acknowledgments as per your specific project details.
