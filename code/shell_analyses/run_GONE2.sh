#Filter large chromosomes from Crab genome for use in GONE2

#first filter contigs larger than 1 million bases (there aren't many)
awk '$2 > 1000000' SnowCrabGenome.fasta.fai > SnowCrabLargeContigs.fai 

#next create a list or file of the names
 awk '{print $1}' SnowCrabLargeContigs.fai > chromies.txt

#then use these names in plink to filter

plink --vcf allcrub.vcf.gz \
--chr $(<chromies.txt) \
--recode \
--out crab_1Mb_sizefiltered \
--allow-extra-chr \
--double-id

#this led to 1145 individuals with 443529 variants at a genotyping rate of 0.856

#finally, load these into GONE2 so we don't have a problem with too many tiny chromosomes

/home/mcrg/GONE2/gone2 crab_1Mb_sizefiltered.ped -t 40 -g 3 
#in gone, -g 3 indicates low coverage genome sequencing of unphased diploids

#AND IT FAILED - all chromsomes must be larger than 20cM which means this genome is too fragmented to use
