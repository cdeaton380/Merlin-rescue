#!/bin/bash

# Performs quality control checks on fastq files using FASTQC and aggregates the results using MultiQC.

# Directory containing fastq files
fastqDir="/path/to/fastq/files"

# Output directory for FASTQC results
outputDir="/path/to/fastqc/output"

# Run FASTQC on all fastq files in the specified directory
fastqc -o ${outputDir} ${fastqDir}/*.fastq.gz

# Run MultiQC to aggregate the FASTQC reports
multiqc ${outputDir} -o ${outputDir}
