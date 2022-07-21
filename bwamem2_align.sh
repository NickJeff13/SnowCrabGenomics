#!/bin/bash

#align my files to the reference genome using bwa-mem2 and samtools

#index reference genome from NCBI PRJNA602365
#Move to the Snow_Crab_PoolSeq folder

./bwa-mem2 index -p SnowCrabIndex ../Snow_Crab_PoolSeq/Genome/SnowCrab/SnowCrabGenome.fasta

#index genome for samtools 

./samtools faidx /hdd2/Snow_Crab_PoolSeq/Genome/Csap_genome.fna

#Create Picard Sequence Dictionary - this is in my reference genome folder
cd /hdd2/Snow_Crab_PoolSeq/Genome/
java -jar ../../picard.jar CreateSequenceDictionary --REFERENCE SnowCrabGenome.fasta --OUTPUT SnowCrabGenome.fasta.dict


#align poolseq reads to reference genome using a for loop
#wanted to use parallel but couldn't figure out the error

cd /hdd2/Snow_Crab_PoolSeq/Trimmed
: ' 
added the Read Groups flag in this version (-R) which will hopefully carry over to the de-duped bam files next
$outfile will add the pool number to the Read Group line, removed the {} which is used with parallel, not a for loop
'
for f1 in *R1.trimmed.fastq.gz;
  do outfile=${f1%%_R1.trimmed.fastq.gz}"" ;
  echo $outfile\.bam ;
  ../../bwa-mem2/bwa-mem2 mem \
  -t 32 \
  -R "@RG\tID:$outfile\tSM:$outfile\tLB:SnowCrab" \
   ../Genome/SnowCrab/SnowCrabGenome.fasta \
  $outfile\_R1.trimmed.fastq.gz  $outfile\_R2.trimmed.fastq.gz\
  | samtools sort -o $outfile\.sorted.bam -T $outfile -@ 32 -m 3G ;
  done
&&
samtools flagstat Pool-1.sorted.bam

#check the stats and make a text file for each
#samtools flagstat Pool1.sorted.bam > Pool1Results.txt



  
