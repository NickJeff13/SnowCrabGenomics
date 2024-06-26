#!/bin/bash

projdir=/mnt/sdb/SnowCrabRNA/AllRawReads/

cd $projdir

#export to make parallel happy
export fastp
export projdir
export Sample


#Run fastp on all the files using parallel
#First list all the files
ls -1 | sed -e 's/\_R..fastq.gz$//' | sort | uniq | \
parallel -j 10 \
'fastp -i {}_R1.fastq.gz  -I {}_R2.fastq.gz \
-q 20 \
-o {}_R1.trimmed.fastq.gz \
-O {}_R2.trimmed.fastq.gz \
--thread 12 '

