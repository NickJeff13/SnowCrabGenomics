#subset 30,000 random SNPs for NeEstimator from maf filtered plink files

#use plink to thin a file 
#0.016 is the proportion of 30,000 SNPs of 1.8 million from the pool filtered file
#then recode to plink ped format
plink --bfile snowcrab.poolfiltered --thin 0.016 --recode --out subset_random --allow-extra-chr

#Now use PGDspider to convert ped to genepop for NeEstimator. Load genepop file into NeEstimator and select analyses for Ne
