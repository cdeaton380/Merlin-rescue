#!/bin/bash

# Trimming adapters and filtering for quality in RNA-seq data using Cutadapt.

fastqFiles='list_fastq.txt'  # List of fastq file names

# Loop through each fastq file listed and perform trimming and quality filtering
for fastq in $(cat $fastqFiles)
do
    echo "Trimming and filtering sample: $fastq"
    cutadapt -j 10 -q 30,30 -m 20 -a ADAPTER_FWD -A ADAPTER_REV -o "${fastq}_trimmed.1.fastq.gz" -p "${fastq}_trimmed.2.fastq.gz" "${fastq}_R1.fastq.gz" "${fastq}_R2.fastq.gz"
    echo "Completed sample: $fastq"
done
