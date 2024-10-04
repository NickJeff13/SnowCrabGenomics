

#make a binary ped file in plink 1.9
plink --vcf allcrubs.exon.vcf --make-bed --out 

#recode all values in the bim file to zero to make Admixture happy
awk '{for(i=1;i<=NF;i++) $i=0}1' exon.recode.bim > exon.recode2.bim
mvf exon.recode2.bim exon.recode.bim

#run Admixture which is in our PATH
% for K in {1..5}
	do
		admixture --cv exon.bed -j40 $K | teee log${K}.out
	done
