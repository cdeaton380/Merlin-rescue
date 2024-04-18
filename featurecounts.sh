#!/bin/bash

# Runs FeatureCounts to count features from BAM files.


# Running FeatureCounts to count features from BAM files
featureCounts -p -t exon -g gene_name --extraAttributes gene_id,gene_type -T 10 -a path_to_annotation_file -o FeatureCounts_Output *.sorted.bam > featurecounts.screen-output.log
