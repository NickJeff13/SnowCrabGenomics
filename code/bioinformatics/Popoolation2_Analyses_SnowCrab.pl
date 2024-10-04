#Start doing some stats on my sync file
#allele frequencies, heterozygosity, and Fst
#in the sync file the allele frequencies are in the format A:T:C:G:N:del

#1. calculate allele frequency differences
perl ../../../../../home/mcrg/popoolation2_1201/snp-frequency-diff.pl --input SnowCrabPopool2_Q25.sync --output-prefix AllPoolsAlleles --min-count 4 --min-coverage 20 --max-coverage 500


#2. Fst for every SNP
perl ../../../../../home/mcrg/popoolation2_1201/fst-sliding.pl --input ZSnowCrabPopool2_Q25.sync --output AllSNPs.fst --suppress-noninformative --min-count 4 --min-coverage 20 \
--max-coverage 500 --min-covered-fraction 1 --window-size 1 --step-size 1 \
--pool-size 47:47:46:48:47:47:47:47:41:39:49:48:47:47:47:47:47:30:30:54:50:40:30
#the window size of 1 here lets us calculate Fst for every SNP


#3. Fst using a sliding window approach of 1000
perl ../../../../../home/mcrg/popoolation2_1201/fst-sliding.pl --input SnowCrabPopool2_Q25.sync --output AllPools_w1000.fst --min-count 4 --min-coverage 20 --max-coverage 500 \
--min-covered-fraction 0.2 --window-size 1000 --step-size 1000 --suppress-noninformative \
--pool-size 47:47:46:48:47:47:47:47:41:35:45:45:45:45:45

#4. Calculate Fisher's exact test
date
perl ../../../../../home/mcrg/popoolation2_1201/fisher-test.pl --input SnowCrabPopool2_Q25.sync --output AllPoolsFisher.fet --suppress-noninformative --min-count 4 --min-coverage 20 --max-coverage 500 --min-covered-fraction 0.2 --window-summary-method geometricmean

#5. Using popoolation here - measure pi for each population mpileup separately
perl ../../../../../home/mcrg/popoolation_1.2.2/Variance-sliding.pl --measure pi \
 --input Pool14.indelrealigned_mpileup --output Pool14.pi --min-count 4 \
 --min-coverage 20 --max-coverage 500 --window-size 10000 --step-size 10000 \
 --pool-size 47 \
 --fastq-type sanger &&
 
 perl ../../../../../home/mcrg/popoolation_1.2.2/Variance-sliding.pl --measure pi \
 --input Pool15.indelrealigned_mpileup --output Pool15.pi --min-count 4 \
 --min-coverage 20 --max-coverage 500 --window-size 10000 --step-size 10000 \
 --pool-size 47 \
 --fastq-type sanger &&
 
 perl ../../../../../home/mcrg/popoolation_1.2.2/Variance-sliding.pl --measure pi \
 --input Pool16.indelrealigned_mpileup --output Pool16.pi --min-count 4 \
 --min-coverage 20 --max-coverage 500 --window-size 10000 --step-size 10000 \
 --pool-size 47 \
 --fastq-type sanger &&
 
 perl ../../../../../home/mcrg/popoolation_1.2.2/Variance-sliding.pl --measure pi \
 --input Pool17.indelrealigned_mpileup --output Pool17.pi --min-count 4 \
 --min-coverage 20 --max-coverage 500 --window-size 10000 --step-size 10000 \
 --pool-size 47 \
 --fastq-type sanger &&
 
 perl ../../../../../home/mcrg/popoolation_1.2.2/Variance-sliding.pl --measure pi \
 --input Pool18.indelrealigned_mpileup --output Pool18.pi --min-count 4 \
 --min-coverage 20 --max-coverage 500 --window-size 10000 --step-size 10000 \
 --pool-size 30 \
 --fastq-type sanger &&
 
 perl ../../../../../home/mcrg/popoolation_1.2.2/Variance-sliding.pl --measure pi \
 --input Pool19.indelrealigned_mpileup --output Pool19.pi --min-count 4 \
 --min-coverage 20 --max-coverage 500 --window-size 10000 --step-size 10000 \
 --pool-size 30 \
 --fastq-type sanger &&
 
 perl ../../../../../home/mcrg/popoolation_1.2.2/Variance-sliding.pl --measure pi \
 --input Pool20.indelrealigned_mpileup --output Pool20.pi --min-count 4 \
 --min-coverage 20 --max-coverage 500 --window-size 10000 --step-size 10000 \
 --pool-size 54 \
 --fastq-type sanger &&
 
 perl ../../../../../home/mcrg/popoolation_1.2.2/Variance-sliding.pl --measure pi \
 --input Pool21.indelrealigned_mpileup --output Pool21.pi --min-count 4 \
 --min-coverage 20 --max-coverage 500 --window-size 10000 --step-size 10000 \
 --pool-size 50 \
 --fastq-type sanger &&
 
 perl ../../../../../home/mcrg/popoolation_1.2.2/Variance-sliding.pl --measure pi \
 --input Pool22.indelrealigned_mpileup --output Pool22.pi --min-count 4 \
 --min-coverage 20 --max-coverage 500 --window-size 10000 --step-size 10000 \
 --pool-size 40 \
 --fastq-type sanger &&
 
 perl ../../../../../home/mcrg/popoolation_1.2.2/Variance-sliding.pl --measure pi \
 --input Pool23.indelrealigned_mpileup --output Pool23.pi --min-count 4 \
 --min-coverage 20 --max-coverage 500 --window-size 10000 --step-size 10000 \
 --pool-size 30 \
 --fastq-type sanger 


###Now visualize the pi outputs
for i in *.pi;
do
perl ../../../../../home/mcrg/popoolation_1.2.2/Visualise-output.pl --input $i --output $i.pdf --ylab pi --chromosomes "Chr01 Chr02 Chr03 Chr04 Chr05 Chr06"
done 

