#!/usr/bin/bash


#First filter the VCF for maf and missing data using vcftools

vcftools --gzvcf allcrub.vcf.gz \
 --maf 0.05 \
 --max-missing 0.9 \
 --recode \
 --out allcrabs.maffiltered.missingfiltered
 

#run plink on snow crab vcf 
VCF=/mnt/sdb/SnowCrab_LCWGS/allcrabs.maffiltered.missingfiltered.vcf.gz
EXON=/mnt/sdb/SnowCrab_LCWGS/allcrub.exon.dp2.recode.vcf.gz

#plink is in our bash source so it will run from anywhere
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids

#now a PCA
plink --vcf $VCF --double-id --allow-extra-chr \
--make-bed --pca --out snowcrab.poolfiltered --mind 0.3
#--set-missing-var-ids --extract snowcrab.prune.in \

#compare with Exon data - dont need the --mind flag as missing inds have been removed for this data
plink --vcf $EXON --double-id --allow-extra-chr \
--make-bed --pca --out exon


# Let's also filter the full SNP dataset with a list of SNPs that overlap with the poolseq data (about 1.8 million SNPs)
#the intersect SNPs file is a tab delim file with Chromosome as column 1 and SNP position as column 2

vcftools --gzvcf vcfs/allcrub.vcf.gz \
--positions Poolseq_LCWGS_IntersectingSNPS.txt \
--out PoolSNPS.only
