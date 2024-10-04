#!/usr/bin/bash

#align HiC reads to genome from NCBI and create a sorted bam file
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
../../../home/mcrg/JupiterPlot-master/jupiter \
name=SnowCrab \
ref=SnowCrabGenome.fasta \
fa=Scaffolded/scaffolds_final.fa

##################################################
##Try SALSA to see if it's more intuitive!#######
#################################################
#first make our HiC bam alignment file into bed format using bedtools
#assuming in crab fastq folder
../../../home/mcrg/bedtools2/bin/bamToBed -i CrabHiC.sorted.bam > HiCalignment.bed
sort -k 4 HiCalignment.bed > tmp && mv tmp HiCalignment.bed

#create an index of the genome with samtools faidx if not done already
#now run SALSA with the indexed assembly, and the HiC reads (the bed file) and enzyme. m=yes in this code allows for correcting assembly errors with the HiC data
#-e is the enzyme used for restriction, GATC,GANTC are the Arima enzyme sites
python2 ../../../home/mcrg/SALSA-master/run_pipeline.py -a SnowCrabGenome.fasta -l SnowCrabGenome.fasta.fai -b HiCalignment.bed -e GATC,GANTC -o HiCscaffolds -m yes



#next step is to generate HiC contact matrix using juicer_tools
java -jar -Xmx128G juicer_tools_1.22.01.jar pre alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes && mv out.hic.part out.hic
#the out.hic file can be loaded into Juicebox for visualization

#try Juicer
#in home/mcrg/Juicer
./scripts/juicer.sh -d /mnt/sdb/juicer \
-p references/sizes.genome \
-y restriction_sites/SnowCrabGenome_Arima.txt \
-z references/SnowCrabGenome.fasta \
-t 60 \
-D /mnt/sdb/juicer
