# Nextflow Tracer Use-Case Examples

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/tracer-pod/utility-pod)

This repository provides an example of using Nextflow with our tracking & observability tool - Tracer.

To install Tracer on your initialized Gitpod workspace:

--> Create your account on Tracer 
--> Copy & Paste one-line installation script at https://app.tracer.bio/app/setup

# Example 1: Nextflow RNASEQ training example with indexing & quantification

1. INDEX and QUANTIFICATION processes are run and tracked by Tracer
2. "nextflow run ./rnaseq/rnaseqv1.nf"

# Example 2: Nextflow RNASEQ training example with indexing & quantification - Bacterial GEO dataset

1. INDEX and QUANTIFICATION processes are run and tracked by Tracer using GEO E.coli dataset
2. Details of dataset available at "./rnaseq/data/bacteria/datasets.txt"
2. "nextflow run ./rnaseq/rnaseqv2.nf"

# Example 3: Nextflow RNASEQ training example with indexing & quantification - Human GEO dataset

1. INDEX and QUANTIFICATION processes are run and tracked by Tracer using GEO human AML dataset
2. Details of dataset available at "./rnaseq/data/human/datasets.txt"
2. "nextflow run ./rnaseq/rnaseqv3.nf"

# Example 4: Shell script tracking for STAR mapper with RNASEQ data

1. Generate a genome index and performs mapping of sample dataset using STAR
2. Details of dataset available at "./rnaseq/data/human/datasets.txt"
2. "sh ./rnaseq/STAR.sh"

    ## Files

    - `rnaseqv1.nf`, `rnaseqv2.nf`, `rnaseqv3.nf` : The main workflow script files.
    - `index.nf`: Contains the INDEX process definition.
    - `quantpe.nf`: Contains the QUANTIFICATION process definition.

    ## Workspace

    - Conda and salmon are pre-configured in the workspace using docker.
    - Nextflow will be pre-configured in the workspace using docker.
    - JAVA will be pre-configured in the workspace using docker.
    - STAR will be pre-configured in the workspace using docker.  
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
