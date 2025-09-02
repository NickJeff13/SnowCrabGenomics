#Updated script for running fastp on raw snow crab fastq files, set as 11 sets of individuals 
#running this script as is in the snowcrab directory will submitt 11 jobs using slurm and the 01_fastp script

for i in {00..11}  ;  do sbatch --export=ALL,set=snowcrab.set$i,paramfile=WGSparams_snowcrab.tsv 01_fastp_parallel_ubu2204.sh ;  done

#use squeue -u nij000 to check job status
#use sacct -u <your_user_name> to look at previous jobs history
#use screen -S <command or session> to keep a screen running even if we disconnect from ssh

#Once all are trimmed we can jellyfish

conda activate jellyfish

ls *.fasta.gz | xargs -n 1 echo gunzip -c > generators

cat generators | \
 parallel -j 16 \
 'jellyfish count -g {} -G 4 -m 21 -s 100M -t 10 -C'