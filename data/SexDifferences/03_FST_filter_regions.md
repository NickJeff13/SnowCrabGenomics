Filter plink files to extract data for region of interest (JACEEZ010007791.1; 31131-31258 bp) - extract data from ALL individuals so need to filter full vcf and update IDs for all samples first
```
./plink --vcf /filepath/snowcrab/angsd_out/allcrub.vcf.gz --maf 0.01 --missing --allow-extra-chr --geno 0.3 --mind 0.3 --update-ids /filepath/snowcrab/plink/update_IDs_pops_2.txt --out /gpfs/fs7/grdi/genarcc/wp3/snowcrab/plink/snowcrab_allsnpsind_maf001 --double-id --make-bed
```
Filter for sex region
```
./plink --bfile /filepath/snowcrab/plink/snowcrab_allsnpsind_maf001  --allow-extra-chr --chr JACEEZ010007791.1 --make-bed -out /filepath/snowcrab/plink/snowcrab_allinds_JACEEZ010007791_sexregion  --from-bp 31131 --to-bp 31258 --recode
```
Filter for top SNP (position 31217 bp) use --recode to get genotype .ped file for plotting
```
./plink --bfile /filepath/snowcrab/plink/snowcrab_allsnpsind_maf001  --allow-extra-chr --chr JACEEZ010007791.1 --make-bed -out /filepath/snowcrab/plink/snowcrab_allinds_JACEEZ010007791_sexSNP_recode  --from-bp 31217 --to-bp 31217 --recode
```
Download the data for R script
```
scp username@address:/filepath/snowcrab/plink/snowcrab_allinds_JACEEZ010007791_sexSNP.ped /Users/Name/Downloads/
```
