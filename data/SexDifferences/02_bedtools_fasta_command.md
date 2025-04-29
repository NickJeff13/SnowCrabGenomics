
Got fasta sequence for chromosome of interest (https://www.ncbi.nlm.nih.gov/nuccore/JACEEZ010007791.1) and uploaded to gpsc
```
scp /Users/Name/Desktop/JACEEZ010007791sequence.fasta username@address:/filepath/snowcrab/plink/
```
installed and used bedtools to extract the sequence from  fasta file
The bed file (bedtools_sexloc.txt) used was in this format to indicate chr and position of sequence:
```
JACEEZ010007791.1	31131	31258
```
bedtools script to getfasta
```
bedtools getfasta -fi /filepath/snowcrab/plink/JACEEZ010007791sequence.fasta  -bed /filepath/snowcrab/plink/bedtools_sexloc.txt -fo /filepath/snowcrab/plink/sex_sequence.txt
```
Download results and blast
```
scp username@address:/filepath/snowcrab/plink/sex_sequence.txt /Users/Name/Downloads/
```
