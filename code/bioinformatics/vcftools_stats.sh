#!/usr/bin/bash

#Run some VCFtools analyses for our snow crab vcfs

#get sequencing depth for individuals
vcftools --gzvcf \
	--depth 
	--out 
	
	
#Heterozygosity and inbreeding coefficient
vcftools --gzvcf  \
	--het
	--out
	
	
	
#Tajima's D - best to run on unfiltered dataset to prevent upward bias when minor alleles are removed
vcftools --gzvcf  \
	--TajimaD
	--out


#Relatedness
vcftools --vcf allcrabs.maffiltered.missingfiltered2.recode.vcf\
	--relatedness
	--out allcrabs.maffiltered
