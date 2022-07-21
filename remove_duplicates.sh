#!/bin/bash

# Using gatk and Picard, remove duplicates for each snow crab pool

export gatk

#deduplicate in parallel - this syntax is specific to the newer versions of Picard and gatk 4.2.5.0
cat crabsortedbams.txt | \
parallel --jobs 15 '../../gatk-4.2.5.0/gatk --java-options "-Xmx20G" \
MarkDuplicates \
INPUT={} \
OUTPUT={.}.deDup.bam \
M={.}.deDupMetrics.txt \
REMOVE_DUPLICATES=true
ASSUME_SORT_ORDER=coordinate \
CREATE_INDEX=true' #this command seems to be deprecated

#try validating one of them

java -jar ../../../picard.jar ValidateSamFile -I Pool-1.sorted.DeDup.bam -MODE SUMMARY

