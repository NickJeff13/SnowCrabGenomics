
cd /hdd2/Snow_Crab_PoolSeq/Trimmed

samtools mpileup -B -f ../Genome/SnowCrab/SnowCrabGenome.fasta Pool-1.sorted.deDupRG.indelrealigned.bam Pool-2.sorted.deDupRG.indelrealigned.bam Pool-3.sorted.deDupRG.indelrealigned.bam Pool-4.sorted.deDupRG.indelrealigned.bam Pool-5.sorted.deDupRG.indelrealigned.bam \
Pool-6.sorted.deDupRG.indelrealigned.bam Pool-7.sorted.deDupRG.indelrealigned.bam Pool-8.sorted.deDupRG.indelrealigned.bam Pool-9.sorted.deDupRG.indelrealigned.bam Pool-10.sorted.deDupRG.indelrealigned.bam \
Pool-11.sorted.deDupRG.indelrealigned.bam Pool-12.sorted.deDupRG.indelrealigned.bam Pool-13.sorted.deDupRG.indelrealigned.bam Pool-14.sorted.deDupRG.indelrealigned.bam Pool-15.sorted.deDupRG.indelrealigned.bam > SnowCrabAllPools.mpileup

#here the -B flag is for base alignment quality but we don't want to filter based on this so use this flag
#-f is the reference genome again - this isn't necessarily needed here but we can use it so we know the actual reference nucleotide instead of just N 

 #Let's also generate one mpileup per pool to look at nucleotide diversity later
 ls *indelrealigned.bam | parallel -j15 'samtools mpileup -B -f ../Genome/SnowCrabGenome.fasta {} > {.}_mpileup'
 
 #Next we will create the sync file in using the new package Grenedalf
 #If we're in the directory where my mpile up is, don't need to specify full path under --pileup-path, just the file name
 
../../../home/mcrg/grenedalf/bin/grenedalf sync-file --pileup-path Trimmed/SnowCrabAllPools.mpileup \
--pileup-min-base-qual 35 \
--pileup-quality-encoding sanger \
--out-dir Grenedalf_Analyses \
--file-prefix SnowCrab35 --verbose --threads 15
#This took about 12 hours to complete but didn't use a ton of ram 

ls *_mpileup | parallel -j 23 'java -ea -Xmx8g -jar ../../../../../home/mcrg/popoolation2_1201/mpileup2sync.jar \
  --input {} --output {.}.sync \
  --fastq-type sanger --min-qual 20 --threads 2'
  
 #for gnu parallel {} is the default input string, and {.} is used to mean we want the output to have the same name as the input string, minus its extension

#Compare Grenedalf output to Popoolation2
java -ea -Xmx80g -jar ../../../../../home/mcrg/popoolation2_1201/mpileup2sync.jar \
  --input SnowCrabAllPools.mpileup --output SnowCrabPopool2.sync \
  --fastq-type sanger --min-qual 35 --threads 24 
