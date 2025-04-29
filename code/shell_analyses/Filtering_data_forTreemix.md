Filtering full vcf (9M SNPs) to reduced dataset. Fitlering in plink involved LD pruning, MAF filter, missingness (<30% at SNP and ind), and thinning to <100,000 SNPs.

```
./plink --vcf /filepath/snowcrab/angsd_out/allcrub.vcf.gz --thin-count 100000 --maf 0.05 --missing --indep-pairwise 50 10 0.2 --make-bed --allow-extra-chr --geno 0.3 --mind 0.3 --out /filepath/snowcrab/plink/allcrab_filter_100k --double-id
```

With output bed file - download to local computer
```
scp username@address:/filepath/snowcrab/plink/ /filepath/Downloads/
```

Run PCAdapt to check data (R Script)

Update individual IDs (note had to add SNP names to file; just SNP_1 to SNP_X)
```
./plink --bfile ~/path/snowcrab/allcrab_filter_100k --out ~/path/snowcrab/allcrab_filter_100kupdateid --update-ids ~/path/snowcrab/update_IDs_regions.txt --freq --family --allow-extra-chr --make-bed
```


Gzip allele frequency file for treemix
```
gzip ~/path/snowcrab/allcrab_filter_100kupdateid.frq.strat
```

Move to treemix directory 
run plink2treemix script to convert format 

```
./plink2treemix.py ~/path/snowcrab/allcrab_filter_100kupdateid.frq.strat.gz snowcrab_lcWGS_100k_filter_regions.gz
```

Run treemix. Used 0 to 6 migration events with 1000 bootstraps, no sample size correction (noss), no outgroup specified.
```
for i in {0..6}
do
./treemix -i snowcrab_lcWGS_100k_filter_regions.gz -m $i -o snowcrab_lcWGS_100k_filter_regions.gz.$i -bootstrap -k 1000 -noss > treemix_${i}_log &
done
```

Treemix output files are in 'Treemix_out' folder
