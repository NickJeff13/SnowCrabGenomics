#!/bin/bash

#KMERSGWAS running scripts for Snow Crab PoolSeq data

cd /hdd2/kmerGWAS/external_programs/

#run KMC with canonized kmer counts
./kmc_v3 -t8 -k31 -ci2 @input_files2.txt output_kmc_canon ./ 1> kmc_canon.1 2> kmc_canon.2

#next, run KMC with all kmer counts output
./kmc_v3 -t8 -k31 -ci0 -b  @input_files2.txt output_kmc_all ./ 1> kmc_all.1 2> kmc_all.2

#combine the information from the two KMC runs into one list of kmers
./kmers_add_strand_information -c output_kmc_canon -n output_kmc_all -k 31 -o kmers_with_strand

#remove large kmc files
rm *.kmc*


