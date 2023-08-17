#Snow Crab genome assembly using Illumina reads provided by VGP
#begin with filtering and trimming with fastp
fastp -i iChiOpi2-CL11_1.fastq.gz -I iChiOpi2-CL11_2.fastq.gz -o iChiOpi2-CL11_1TRIMMEd.fastq.gz -O iChiOpi2-CL11_2TRIMMED.fastq.gz
'
Read1 before filtering:
total reads: 547653149
total bases: 82695625499
Q20 bases: 78378994465(94.7801%)
Q30 bases: 74205120989(89.7328%)

Read2 before filtering:
total reads: 547653149
total bases: 82695625499
Q20 bases: 77184522984(93.3357%)
Q30 bases: 72336738234(87.4735%)

Read1 after filtering:
total reads: 523279706
total bases: 77074479648
Q20 bases: 74419552779(96.5554%)
Q30 bases: 70893009213(91.9799%)

Read2 after filtering:
total reads: 523279706
total bases: 77073802971
Q20 bases: 73763657104(95.7052%)
Q30 bases: 69637010201(90.3511%)

Filtering result:
reads passed filter: 1046559412
reads failed due to low quality: 48241034
reads failed due to too many N: 38368
reads failed due to too short: 467484
reads with adapter trimmed: 61228590
bases trimmed due to adapters: 3649985174

Duplication rate: 18.7795%

Insert size peak (evaluated by paired-end reads): 151

JSON report: fastp.json
HTML report: fastp.html

fastp -i iChiOpi2-CL11_1.fastq.gz -I iChiOpi2-CL11_2.fastq.gz -o iChiOpi2-CL11_1TRIMMEd.fastq.gz -O iChiOpi2-CL11_2TRIMMED.fastq.gz 
fastp v0.23.2, time used: 3266 seconds
'

#try an assembly of the trimmed fastq's using Clover and compare to assembly in NCBI which I will scaffold with HiC reads
cd /home/clover-2.0
#Start August 9 16:00
./clover -k 40 \
-i1 /mnt/sdb/Snow_Crab_Genome/iChiOpi2-CL11_TRIMMED1.fastq \
-i2 /mnt/sdb/Snow_Crab_Genome/iChiOpi2-CL11_TRIMMED2.fastq \
-o SnowCrabIllum
