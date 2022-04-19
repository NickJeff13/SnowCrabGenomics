#!/bin/bash

#align my files to the reference genome using bwa-mem2 and samtools

#index reference genome from NCBI PRJNA602365
#Move to the Snow_Crab_PoolSeq folder

./bwa-mem2 index -p SnowCrabIndex ../Snow_Crab_PoolSeq/Genome/SnowCrabGenome.fasta

#index genome for samtools 

./samtools fqidx /hdd2/Snow_Crab_PoolSeq/Genome/SnowCrabGenome.fasta -o SCSamtoolsIndex

#Create Picard Sequence Dictionary - this is in my reference genome folder
cd /hdd2/Snow_Crab_PoolSeq/Genome/
java -jar ../../picard.jar CreateSequenceDictionary --REFERENCE SnowCrabGenome.fasta --OUTPUT SnowCrabGenome.fasta.dict


#align poolseq reads to reference genome using a for loop
#wanted to use parallel but couldn't figure out the error

cd /hdd2/Snow_Crab_PoolSeq/Trimmed
#added the Read Groups flag in this version (-R) which will hopefully carry over to the de-duped bam files next
#$outfile will add the pool number to the Read Group line, removed the {} which is used with parallel, not a for loop
for f1 in *R1Trimmed.fastq.gz;
  do outfile=${f1%%_R1Trimmed.fastq.gz}"" ;
  echo $outfile\.bam ;
  ../../bwa-mem2/bwa-mem2 mem \
  -t 32 \
  -R "@RG\tID:$outfile\tSM:$outfile\tLB:SnowCrab" \
   ../Genome/SnowCrabGenome.fasta \
  $outfile\_R1Trimmed.fastq.gz  $outfile\_R2Trimmed.fastq.gz\
  | samtools sort -o $outfile\.sorted.bam -T $outfile -@ 32 -m 3G ;
  done


#check the stats and make a text file for each
#samtools flagstat Pool1.sorted.bam > Pool1Results.txt



  
