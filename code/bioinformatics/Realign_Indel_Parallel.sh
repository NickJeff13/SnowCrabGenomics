#!/bin/bash

#After de-duping, need to add read groups to each pool, and index the bams, THEN
#Realign Target in Parallel 
#Followed by Indel realignment in parallel

#go to the de-duplicated bam files
cd /hdd2/Snow_Crab_PoolSeq/Trimmed

#Switch to older java with galternatives
ls *.sorted.deDup.bam | parallel --jobs 15 'java -Xmx16g -jar ../../picard.jar AddOrReplaceReadGroups \
I={} \
O={.}RG.bam \
SORT_ORDER=coordinate \
RGID={.} \
RGLB={.} \
RGPL=illumina \
RGPU={.}.unit RGSM={.}' &&
parallel samtools index ::: *.deDupRG.bam
   
#Run the below line on your bam to see if it contains read groups (RG) now
#java -jar ../../../picard.jar AddOrReplaceReadGroups I=Pool7.sorted.bamdeDup.bam O=Pool7.sorted.DeDupRG.bam RGID=7 RGLB=Pool7 RGPL=illumina RGPU=unit7 RGSM=Pool7
#samtools view -H Pool7.sorted.DeDupRG.bam.bam | grep '^@RG'

#the RealignerTargetCreator program is only in gatk < v.4 so had to use 3.7, which only works on java 1.8, so had to switch java version for now
ls *RG.bam | parallel --jobs 15 ' java -Xmx12g -jar ../../GenomeAnalysisTK.3.7/GenomeAnalysisTK.jar \
 -T RealignerTargetCreator \
 -R ../Genome/SnowCrab/SnowCrabGenome.fasta \
 -I {} \
 -o {.}.intervals'
#the above reference genome needs a Picard dictionary called SnowCrabGenome.dict in this case
&&
#THEN
##Go to the bams
#cd $projdir/align

#get to realigning - this will take a day or two (~40 hours for this project)
ls *RG.bam | \
   parallel --jobs 15 \
 'java -Xmx8G -jar ../../GenomeAnalysisTK.3.7/GenomeAnalysisTK.jar \
 -T IndelRealigner \
 -R ../Genome/SnowCrab/SnowCrabGenome.fasta \
 -I {} \
 -targetIntervals {.}.intervals \
 -o {.}.indelrealigned.bam '
