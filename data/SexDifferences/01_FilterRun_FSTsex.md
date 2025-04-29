On gpsc, use plink to filter snow crab vcf for maf >0.01 and remove loci missing >30% genotypes and individuals missing >30% data. Filter to keep only individuals that have been sexed (BB and WCB pops). Output file has >9M SNPs
note file path removed working within wp3 directory
```
./plink --vcf /filepath/snowcrab/angsd_out/allcrub.vcf.gz --maf 0.01 --missing --make-bed --allow-extra-chr --geno 0.3 --mind 0.3 --keep /filepath/snowcrab/plink/meta_dat_sexind_snowcrab_keep.txt  --out /filepath/snowcrab/plink/snowcrab_sex_maf001 --double-id
```
On filtered plink file (contains only 175 individuals that have been sexed), update ID names so that sex information is labelled in population/family column
```
./plink --bfile /filepath/snowcrab/plink/snowcrab_sex_maf001  --make-bed --allow-extra-chr --geno 0.3 --mind 0.3 --update-ids /filepath/snowcrab/plink/meta_dat_sexind_snowcrab.txt --out /filepath/snowcrab/plink/snowcrab_sex_maf001_updateID
```
Run FST on plink file with --family to get FST between the males vs. female (coded as pop)
```
./plink --bfile /filepath/snowcrab/plink/snowcrab_sex_maf001_updateID  --allow-extra-chr --fst --family --out /filepath/snowcrab/plink/snowcrab_sex_maf001_updateID_fst
```
Data were downloaded from gpsc to run and plot in R scripts (see other scripts)
```
scp username@address:/filepath/snowcrab/plink/snowcrab_sex_maf001_updateID_fst.fst /Users/Name/Downloads/
```
