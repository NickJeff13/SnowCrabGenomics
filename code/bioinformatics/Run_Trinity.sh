#Assemeble a Trinity fasta file using all tissue types and reads
#use the full_cleanup option to remove all intermediate files that we don't need
#This assumes I'm in the director with the Trinity executable

./Trinity --seqType fq --max_memory 200G \
--left Brain_R1.fastq,Gill_R1.fastq,Heart_R1.fastq,Hepato_R1.fastq,Leg_R1.fastq \
--right Brain_R2.fastq,Gill_R2.fastq,Heart_R2.fastq,Hepato_R2.fastq,Leg_R2.fastq \
--CPU 40 \
--trimmomatic --quality_trimming_params ILLUMINACLIP:../Trimmomatic-0.39TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 AVGQUAL:30 \
--output SnowCrabtrinity_AssemblyTrimmed \
--full_cleanup \
--verbose

Statistics:
===========
Trinity Version:      Trinity-v2.12.0
Compiler:             GCC
Trinity Parameters:   --seqType fq --max_memory 200G --left Brain_R1.fastq,Gill_R1.fastq,Heart_R1.fastq,Hepato_R1.fastq,Leg_R1.fastq --right Brain_R2.fastq,Gill_R2.fastq,Heart_R2.fastq,Hepato_R2.fastq,Leg_R2.fastq --CPU 40 --trimmomatic --quality_trimming_params ILLUMINACLIP:../Trimmomatic-0.39TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36 AVGQUAL:30 --output SnowCrabtrinity_AssemblyTrimmed --verbose
Paired mode
 Input data
  Left.fasta    10592 MByte
  Right.fasta   10592 MByte
  Number of unique KMERs: 914343069
  Number of reads:        154829456 Output data
  Trinity.fasta 459 MByte

Runtime
=======
Start:       Fri Jun 25 10:06:39 ADT 2021
End:         Sat Jun 26 10:37:27 ADT 2021
Trinity   88248 seconds
  Inchworm (phase 1 - read clustering)  3508 seconds
  Chrysalis (phase 1 - read clustering) 55047 seconds
  Rest (phase 2 - parallel assembly)       29693 seconds
