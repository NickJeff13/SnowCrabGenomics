#!/usr/bin/bash

##NOTE: INCOMPLETE - could not get the angel_train.py script to work due to some python library issues

#running ANGEL on crab transcriptome

conda activate anaCogent 

#ANGEL scripts are located in the home directory
python dumb_predict.py /mnt/sdb/SnowCrabRNA/AssemblyTrimmed/TranscriptomeTrimmed.fasta crabANGEL \
--min_aa_length 300 \
--cpus 24

'
Total CPU time 16010.96
Longest 500 non-redundant predicted ORFs written to: crabANGEL.nr90.longest_500.cds
Calculating base frequency from /mnt/sdb/SnowCrabRNA/AssemblyTrimmed/TranscriptomeTrimmed.fasta
Calculating hexamer scores from crabANGEL.nr90.longest_500.cds
Scoring predicted ORFs....
Dumb ORF prediction done. Final output written to: crabANGEL.final.cds crabANGEL.final.utr crabANGEL.final.pep
'

#make training dataset
python angel_make_training_set.py crabANGEL crabANGEL.training \
--random \
--cpus 24

#now do the ANGEL classifier training - needs the CDS fasta and UTR fasta
angel_train.py crabANGEL.final.cds crabANGEL.final.utr crabANGEL.classifier.pickle \
--cpus 
