#!/bin/bash

# Using gatk and Picard, remove duplicates for each snow crab pool

export gatk

#deduplicate in parallel - this syntax is specific to the newer versions of Picard and gatk 4.2.5.0
cat crabsortedbams.txt | \
parallel --jobs 15 '/hdd2/gatk-4.2.5.0/gatk --java-options "-Xmx20G" \
MarkDuplicates \
-I {} \
-O {.}.deDup.bam \
-M {.}.deDupMetrics.txt \
--REMOVE_DUPLICATES true
--CREATE_INDEX true'


