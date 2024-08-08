#!/usr/bin/bash

#run plink on snow crab vcf 
VCF=/mnt/sdb/SnowCrab_LCWGS/allcrub.vcf.gz
EXON=/mnt/sdb/SnowCrab_LCWGS/allcrub.exon.dp2.recode.vcf.gz

#plink is in our bash source so it will run from anywhere
plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids

#now a PCA
plink --vcf $VCF --double-id --allow-extra-chr \
--make-bed --pca --out snowcrab.recode --mind
#--set-missing-var-ids --extract snowcrab.prune.in \

#compare with Exon data - dont need the --mind flag as missing inds have been removed for this data
plink --vcf $EXON --double-id --allow-extra-chr \
--make-bed --pca --out exon
