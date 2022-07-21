#!/bin/bash
#Run Grenedalf analyses on the generated crab sync file

#1. create a table of snp frequencies
../../../home/mcrg/grenedalf/bin/grenedalf frequency --sync-path /Grenedalf/SnowCrabcounts.sync --omit-invariants --out-dir /Grenedalf/ --verbose --threads 40
#this only took an hour
#Processed 1561135480 genome positions of the input file, and thereof skipped 1189592784 due to being invariant sites.

#2. Calculate Fst in windows across the genome
poolsizes=42,48,48,48,31,37,24,47,47,28,32,31,48,42,39
../../../home/mcrg/grenedalf/bin/grenedalf fst --sync-path Grenedalf/SnowCrabcounts.sync --window-width 20000 --pool-sizes $poolsizes --method spence-nei --omit-na-windows --verbose --threads 40 --sample-name-list poolnames.txt


#3. Allele frequency spectrum heatmap
../../../home/mcrg/grenedalf/bin/grenedalf afs-heatmap --sync-path Grenedalf/SnowCrabcounts.sync --sample-name-list poolnames.txt \
--window-width 20000 --resolution 100 --max-frequency 1 --spectrum-type unfolded --file-prefix SnowCrabAFS --verbose --threads 40

#4. Diversity statistics
../../../home/mcrg/grenedalf/bin/grenedalf diversity --sync-path Grenedalf/SnowCrabcounts.sync --sample-name-list poolnames.txt \
--window-width 20000  --pool-sizes $poolsizes --measure all --min-allele-count 2 --min-coverage 30 --max-coverage 4000 \
--separator-char tab --file-prefix SnowCrabDiversity  --verbose --threads 40
