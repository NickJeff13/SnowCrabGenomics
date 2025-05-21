#Run STAR to align my RNAseq reads to reference crab genome
cd /hdd2/STAR/source/

#create the genome index using STAR, a fasta file, and a gtf file
./STAR --runThreadN 10 --runMode genomeGenerate --genomeDir ../../SnowCrabRNA/STARAlignment --genomeFastaFiles ../../Snow_Crab_PoolSeq/Genome/SnowCrabGenome.fasta --sjdbGTFfile ../../Snow_Crab_PoolSeq/Genome/snowcrabannotation.gtf 

#Now align the reads to our indexed reference
./STAR --runThreadN 10 --genomeDir ../../SnowCrabRNA/STARAlignment --readFilesIn ../../SnowCrabRNA/Leg/Leg_R1.fastq.gz, ../../SnowCrabRNA/Leg/Leg_R2.fastq.gz --readFilesCommand gunzip

#Now align the transcriptome to the indexed genome and annotation file
./STAR --runThreadN 10 --genomeDir ../../SnowCrabRNA/STARAlignment --readFilesIn ../../SnowCrabRNA/ --readFilesCommand gunzip


minimap2 -a -x splice -t 20 -o rnaaln.sam Snow_Crab_Genome/SnowCrabGenome.fasta SnowCrabRNA/AllRawReads/Trinity.fasta 

samtools view -bS rnaaln.sam | samtools sort -o rnaaln.sorted.bam

/home/mcrg/stringtie-3.0.0.Linux_x86_64/stringtie rnaaln.sorted.bam -o sc.transcripts.gtf


cp sc.transcripts.gtf /SnowCrabRNA
cp sc.transcripts.gtf SnowCrabRNA/AllRawReads/

#in the SnowCrabRNA working directory now
/home/mcrg/gffcompare/gffcompare -r ../../Snow_Crab_Genome/snowcrabannotation.gtf -o compare_out sc.transcripts.gtf 

#summary stats for gffcompare of my transcriptome with genome annotations file
'
#= Summary for dataset: sc.transcripts.gtf 
#     Query mRNAs :   43494 in   41153 loci  (40888 multi-exon transcripts)
#            (2000 multi-transcript loci, ~1.1 transcripts per locus)
# Reference mRNAs :   22645 in   21738 loci  (15857 multi-exon)
# Super-loci w/ reference transcripts:     7861
#-----------------| Sensitivity | Precision  |
        Base level:    42.9     |    16.5    |
        Exon level:    35.4     |    20.8    |
      Intron level:    50.4     |    31.1    |
Intron chain level:     4.1     |     1.6    |
  Transcript level:     2.9     |     1.5    |
       Locus level:     3.0     |     1.6    |

     Matching intron chains:     647
       Matching transcripts:     647
              Matching loci:     646

          Missed exons:   40599/91876	( 44.2%)
           Novel exons:  105683/156019	( 67.7%)
        Missed introns:   23614/69865	( 33.8%)
         Novel introns:   68453/113359	( 60.4%)
           Missed loci:   12359/21738	( 56.9%)
            Novel loci:   32639/41153	( 79.3%)
'
