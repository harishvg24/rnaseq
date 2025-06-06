# rnaseq

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code Style: Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)


## Description

This project provides a comprehensive pipeline for RNA sequencing data analysis, from raw reads to differential expression analysis. It automates the process of quality control, alignment, transcript quantification, and statistical analysis, enabling researchers to efficiently analyze RNA-seq data and identify differentially expressed genes. Key features include:

-   Automated end-to-end workflow
-   Support for various aligners and quantification tools (e.g., STAR, Salmon)
-   Integrated quality control using FastQC and MultiQC
-   Differential expression analysis with DESeq2 or edgeR

## Installation

download the setup.sh and make it executable using command chmod +x setup.sh and run the command ./setup.sh for downloading all the required tools
    
 1. Prepare your input data:

    > Input data should be placed in the `data/raw` directory. The pipeline expects FASTQ files for each sample. Describe the expected file naming convention. For example: The pipeline expects paired-end reads with filenames like `sample1_R1.fastq.gz` and `sample1_R2.fastq.gz`.


    data/raw/
    ├── sample1_R1.fastq.gz
    ├── sample1_R2.fastq.gz
    ├── sample2_R1.fastq.gz
    └── sample2_R2.fastq.gz
        The pipeline is configured using a `config.yaml` file. This file specifies the paths to input data, genome index, and other parameters.


    # Example command to execute the main script
    bash rnaseq_pipeline.sh --config config.yaml
    
2.  Interpret the results:

    > Results are stored in the `results` directory. Key output files include:
    >
    > -   `results/differential_expression.csv`: Differential expression analysis results. This file contains a table of genes, their log2 fold changes, p-values, and adjusted p-values.
    > -   `results/qc_report.html`: Quality control report generated by MultiQC. This report provides an overview of the quality of the reads and the alignment.
    > Describe other important output files and their contents. For example:
    > -   `results/alignment/`: Contains the alignment files (BAM/SAM)
    > -   `results/quantification/`: Contains the quantification files (gene counts, TPM values)
