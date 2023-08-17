#!/usr/bin/bash

#Run sortmerna on my RNA reads per tissue to filter out bacterial sequences
#Installed sortmerna using conda

#Keep 'ref' databases the same, but run separately for each tissue type

sortmerna --ref /home/mcrg/sortmerna-master/data/rRNA_databases/silva-bac-16s-id90.fasta \
--ref /home/mcrg/sortmerna-master/data/rRNA_databases/silva-arc-16s-id95.fasta \
--ref /home/mcrg/sortmerna-master/data/rRNA_databases/silva-bac-23s-id98.fasta \
--reads Gill/NS.1265.002.NEBNext_dual_i7_G12---NEBNext_dual_i5_G12.SC_RNA_Gill_R1.fastq.gz \
--reads Gill/NS.1265.002.NEBNext_dual_i7_G12---NEBNext_dual_i5_G12.SC_RNA_Gill_R2.fastq.gz \
-workdir /hdd2/SnowCrabRNA/Gill/

