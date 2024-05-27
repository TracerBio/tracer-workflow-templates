# Nextflow Tracer Use-Case Example

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/tracerhub/tracer-usecase)

This repository provides an example of using Nextflow with a tracking tool called Tracer. The workflow consists of three main files: `main2.nf`, `index.nf`, and `quantpe.nf`. The key workflow file is `main2.nf`.

## Overview

The workflow performs the following steps:
1. Initializes the Tracer tool.
2. Executes the INDEX process to prepare the transcriptome.
3. Executes the QUANTIFICATION process to quantify the paired-end reads.
4. Terminates the Tracer tool.

The user can run the entire workflow by typing `nextflow run main2.nf`.

## Files

- `main2.nf`: The main workflow script.
- `index.nf`: Contains the INDEX process definition.
- `quantpe.nf`: Contains the QUANTIFICATION process definition.

## Requirements

- Conda and Salmon: Anaconda will be auto-installed along with salmon when you initialise your GitPod workspace.
- Nextflow: Nextflow will be auto-installed on your GitPod workspace) 
- Tracer: The tracking tool used in this workflow. Setup your Tracer account and record your unique API key at https://app.tracer.bio. Enter this API key when prompted for Tracer installation. Exit terminal & Open a new terminal to initiate Tracer.  

## Usage

1. Clone the repository:

    bash
    git clone https://github.com/yourusername/your-repo-name.git
    cd your-repo-name
    

2. Input data:
    - Input files `gut_1.fq` and `gut_2.fq` are present in the `data/ggal/` directory.
    - Transcriptome reference is placed in the `data/ggal/` directory and is named `ref1.fa`.

3. Workflow will be auto-run upon initiating a new terminal post Tracer installation using the following command:

    nextflow run main2.nf

## Workflow Details

### main2.nf

This file contains the main workflow definition:

'groovy'
nextflow.enable.dsl=2

params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome = "$projectDir/data/ggal/ref1.fa"

// Include the INDEX process
include { INDEX } from './index.nf'

// Include the QUANTIFICATION process
include { QUANTIFICATION } from './quantpe.nf'

workflow {
    // Create a channel from the reads
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // Run the INDEX process
    index_ch = INDEX(params.transcriptome)

    // Run the QUANTIFICATION process
    quant_ch = QUANTIFICATION(index_ch, read_pairs_ch)
}


### index.nf

Defines the INDEX process that builds an index from the transcriptome reference:

'groovy'
nextflow.enable.dsl=2

/*
 * Define the INDEX process
 */
process INDEX {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    tracer start
    tracer tool index v1
    salmon index --threads ${task.cpus} -t ${transcriptome} -i salmon_index
    """
}

### quantpe.nf

Defines the QUANTIFICATION process that quantifies the paired-end reads using the index:

'groovy'
nextflow.enable.dsl=2

/*
 * Define the QUANTIFICATION process
 */
process QUANTIFICATION {
    input:
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    tracer tool quant v2
    salmon quant --threads ${task.cpus} --libType=U -i ${salmon_index} -1 ${reads[0]} -2 ${reads[1]} -o ${sample_id}
    tracer end
    """
}

## Tracer Integration

Tracer is integrated into the workflow to track the execution and versions of the processes. The workflow initiates Tracer at the beginning and terminates it at the end. 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Thanks to the Nextflow community for their excellent documentation and support.

## GitPod

[![Open in Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/tracerhub/tracer-usecase) 

---

This README provides a comprehensive guide for users to understand and execute a nextflow RNA-Seq workflow with Tracer enabled tracking. Adjust the repository URL, license details, and acknowledgments as per your specific project details.
