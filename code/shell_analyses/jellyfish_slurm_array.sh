#!/bin/bash

#SBATCH --job-name=jellyfish_parallel
#SBATCH --open-mode=append
#SBATCH --partition=standard
#SBATCH --time=15:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --comment="image=registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04,ssh=true,nsswitch=true"
#SBATCH --account=grdi_genarcc
#SBATCH --array=1-1110

inds.tsv < ls *.fastq.gz

#remove every second line since 2 of each individual using sed
 sed -i '0~2d' inds.tsv

#get conda
source ~/.bashrc

conda activate jellyfish

#define inds for the slurm array
inds=$(sed -n "${SLURM_ARRAY_TASK_ID}p" inds.tsv)

#run jellyfish per array
jellyfish count -m 21 -C -s 2G -o ${inds}_21mer.jf -t 8 <(zcat ${inds}_R1.trimmed.fastq.gz ${inds}_R2.trimmed.fastq.gz)

conda deactivate