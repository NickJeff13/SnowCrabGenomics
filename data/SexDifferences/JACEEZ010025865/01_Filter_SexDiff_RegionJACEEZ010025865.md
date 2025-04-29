Plink Filter for this region which has larger number of high FST values than other chromosomes
```
./plink --bfile /filepath/snowcrab/plink/snowcrab_allsnpsind_maf001  --allow-extra-chr --chr  JACEEZ010025865.1 --update-ids /filepath/snowcrab/plink/update_IDs_pops_2.txt --make-bed -out /filepath/snowcrab/plink/snowcrab_allinds_JACEEZ010025865_sexregionFST  --from-bp 300000 --to-bp 357000
```
download file and run in R for pcadapt
```
scp username@address:/filepath/snowcrab/plink/snowcrab_allinds_JACEEZ010025865_sexregionFST.* /Users/Name/Downloads/
```
