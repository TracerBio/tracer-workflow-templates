# Tracer Workflow Templates
**how to install tracer**
To install Tracer on your initialized Gitpod workspace:

--> Create your account on Tracer 
--> Copy & Paste one-line installation script at https://app.tracer.bio/app/setup

**how to run tracer daemon** 
```bash
tracer-deamon 
```
**How to check if Tracer Daemon Is Running**
```bash
ps -e | grep tracer

```

# 1.1 Nextflow Tracer Workflow Templates

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/tracer-pod/utility-pod)

This repository provides an example of using Nextflow with our tracking & observability tool - Tracer.


# Example 1: Nextflow RNASEQ training example with indexing & quantification

1. INDEX and QUANTIFICATION processes are run and tracked by Tracer
2. "nextflow run ./nextflow-tracer-inline-code/example1.nf"

# Example 2: Nextflow RNASEQ training example with indexing & quantification - Bacterial GEO dataset

1. INDEX and QUANTIFICATION processes are run and tracked by Tracer using GEO E.coli dataset
2. Details of dataset available at "./examples/data/bacteria/datasets.txt"
2. "nextflow run ./nextflow-tracer-inline-code/example2.nf"

# Example 3: Nextflow RNASEQ training example with indexing & quantification - Human GEO dataset

1. INDEX and QUANTIFICATION processes are run and tracked by Tracer using GEO human AML dataset
2. Details of dataset available at "./examples/data/human/datasets.txt"
2. "nextflow run ./nextflow-tracer-inline-code/example3.nf"

# Example 4: Nextflow Quality control training example using FASTQC - Human GEO dataset 

1. Generate data metrics and QC characteristics - GC%, reads per bp etc. for raw data (.fq files)
2. Details of dataset available at "./examples/data/human/datasets.txt"
2. "nextflow run ./nextflow-tracer-inline-code/example4.nf"

# Example 5: Nextflow CHIPSEQ training example with indexing & mapping - Human GEO dataset 

1. Generates a human genome index and maps the CHIPSEQ data (.fq) to the index
2. Details of dataset available at "./examples/data/human/datasets.txt"
2. "nextflow run ./nextflow-tracer-inline-code/example5.nf"

# Example 6: Nextflow complete CHIPSEQ workflow - from index to peaks - Human GEO dataset

1. Generate a genome index and performs mapping of sample dataset using STAR
2. Calls peaks on mapped data using samtools and MACS3
3. Performs mapping coverage analysis using deeptools
4. "nextflow run ./nextflow-tracer-inline-code/chipseq.nf"

# 1.2 Shell Script Tracer Examples

# Example 1: Shell script tracking for STAR mapper with RNASEQ data

1. Generate a genome index and performs mapping of sample dataset using STAR
2. Details of dataset available at "./examples/data/human/datasets.txt"
2. "sh ./shell-tracer-inline-code/example1.sh"

# Example 2: Complete CHIPSEQ processing example

1. Generate a genome index and performs mapping of sample dataset using STAR
2. Calls peaks on mapped data using samtools and MACS3
3. Performs mapping coverage analysis using deeptools
4. "sh ./shell-tracer-inline-code/chipseq.sh"


    ## Files

    - `example1.nf`, `example2.nf`, `example3.nf`, `example1.sh`, `example5.nf`, `example6.nf`, `chipseq.sh`, `chipseq.nf`: The main workflow script files.
    - `index.nf`: Contains the INDEX process definition [./examples/misc/]
    - `quantpe.nf`: Contains the QUANTIFICATION process definition. [./examples/misc/]
    - `genomeGenerate.nf`: Contains the GENOME_GENERATE process definition. [./examples/misc/]
    - `star-align.nf`: Contains the STAR_ALIGN process definition. [./examples/misc/]
    - `fastqc.nf`: Contains the FASTQC process definition. [./examples/misc/]

    ## Workspace

    - Conda and salmon are pre-configured in the workspace using docker.
    - deeptools, samtools are configured into the docker container.
    - Nextflow is pre-configured in the workspace using docker.
    - JAVA is pre-configured in the workspace using docker.
    - STAR is pre-configured in the workspace using docker.
    - MACS3 is pre-configured in the workspace using docker.  
    - Tracer: The tracking tool used in this workflow. Setup your Tracer account and record your unique API key at https://app.tracer.bio. Enter this one-line installation script prior to running other example scripts.  

    ## Usage

    1. To clone the repository:

    git clone https://github.com/tracer-pod/utility-pod.git

    cd utility-pod
    

    ## Tracer Integration

    Tracer is integrated into the workflow to track the execution and versions of the processes. The workflow initiates Tracer at the beginning and terminates it at the end. 

    ## License

    This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

    ## Acknowledgments

    - Thanks to the Nextflow community for their excellent documentation and support.

    ## GitPod

    [![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/tracer-pod/utility-pod)
---

This README provides a comprehensive guide for users to understand and execute a nextflow RNA-Seq workflow with Tracer enabled tracking. Adjust the repository URL, license details, and acknowledgments as per your specific project details.
