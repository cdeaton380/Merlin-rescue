#!/bin/bash

# Performs read mapping using HISAT2 and directly pipes the output to Samtools for conversion and sorting.


# List of trimmed fastq files
fastqFiles='list_fastq.txt'

# Loop through each file in the list and perform mapping and post-processing
for fastq in $(cat $fastqFiles)
do
    echo "Mapping and processing sample: $fastq"

    # Mapping reads to the reference genome with HISAT2 --> Use samtools for conversion and sorting
    hisat2 -p 12 -x path_to_hisat2_index -1 "${fastq}_trimmed.1.fastq.gz" -2 "${fastq}_trimmed.2.fastq.gz" | \
    samtools sort -@ 12 -o "${fastq}.sorted.bam"

    echo "DONE Processing sample: $fastq"
done
