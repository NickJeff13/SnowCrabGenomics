#Run STAR to align my RNAseq reads to reference crab genome
cd /hdd2/STAR/source/

#create the genome index using STAR, a fasta file, and a gtf file
./STAR --runThreadN 10 --runMode genomeGenerate --genomeDir ../../SnowCrabRNA/STARAlignment --genomeFastaFiles ../../Snow_Crab_PoolSeq/Genome/SnowCrabGenome.fasta --sjdbGTFfile ../../Snow_Crab_PoolSeq/Genome/snowcrabannotation.gtf 

#Now align the reads to our indexed reference
./STAR --runThreadN 10 --genomeDir ../../SnowCrabRNA/STARAlignment --readFilesIn ../../SnowCrabRNA/Leg/Leg_R1.fastq.gz, ../../SnowCrabRNA/Leg/Leg_R2.fastq.gz --readFilesCommand gunzip
