#Look for runs of homozygosity (ROH) using bcftools and our vcf file
#below -G means we're using genotypes, -AF-dflt is a defaul alt allele freq default of 0.4 and '-e -' means use all individuals to estimate allele frequency

bcftools roh -G30 \
--AF-dflt 0.4 \
allcrabs.maffiltered.missingfiltered2.recode.vcf \
-e - \
-o maffiltered_roh.txt
