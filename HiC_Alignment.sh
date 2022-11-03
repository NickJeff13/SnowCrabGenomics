#!/usr/bin/bash

./bwa-mem2/bwa-mem2 mem \
Snow_Crab_PoolSeq/Genome/SnowCrabGenome.fasta \
qChiOpi1_S_2476_R1_R1_001.fastq.gz \
qChiOpi1_S_2476_R1_R2_001.fastq.gz -t 40 \
| samtools sort -o CrabHiC.sorted.bam -T CrabHiC -@ 32 -m 4G --write-index


#try scaffolding the HiC reads to our indexed genome
cd /hdd2/Snow_Crab_Genome/
#this directory contains my raw HiC files, the aligned bam file, and the indexed reference genome from NCBI
../../../home/mcrg/pin_hic/src/pin_hic_it \
-i 3 \
-x SnowCrabGenome.fasta.fai \
-r SnowCrabGenome.fasta \
-O Scaffolded \
CrabHiC.sorted.bam

': Ouput
[M::cal_asmm] assembly metrics
[M::cal_asmm] Genome Size: 2002919378 bp
[M::cal_asmm] No. scaffolds: 26514
[M::cal_asmm] Maximum: 2536572 bp
[M::cal_asmm] N50: 208145 bp
[M::mk_brks] print average coverage for each 1024 base of the contigs
'

#visualize in JupiterPlot
#Still in Crab Genome folder
../../home/mcrg/JupiterPlot-49e7a27/jupiter \
name=SnowCrab \
ref=SnowCrabGenome.fasta \
fa=Scaffolded/scaffolds_final.fa


