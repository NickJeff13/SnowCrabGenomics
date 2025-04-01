##########################################################################################
############Redundancy Analyses of Genotypes and Environment############################
## Written by Nick Jeffery in winter 2025

# Load libraries ----------------------------------------------------------
library(adegenet)
library(hierfstat)
library(vegan)
library(dplyr)
library(lfmm)
library(qvalue)
library(ggplot2)

#We need 1) genotypes per individual, 2) environmental data per site and individual, and 3) GPS coordinates per individual
# We will be using the exome dataset as the full snp dataset will take a long time


# Read in plink genotypes -------------------------------------------------
# these large files are on my local Linux computer so the file paths only work for me, but the data will be available on NCBI and other sources later
# first converted imputed exon snps file to RAW with the following stats: 1072 individuals, 134,087 variants passed QC, genotyping rate is 0.9162

exon.snps <- read.PLINK("/mnt/sdb/SnowCrab_LCWGS/vcfs/rawexon.raw", parallel = TRUE)

# find clusters while we have this loaded to see if there are any
#Note, this took 3 days without finishing so we'll skip it
#grp <- find.clusters(exon.snps, n.pca = 10, method = "kmeans", scale = FALSE)

#exon.snps <- read.table("/mnt/sdb/SnowCrab_LCWGS/vcfs/rawexon.raw", header=T)
dim(exon.snps)
sum(is.na(exon.snps)) # zero NAs so we're good to use this for RDA 

exon.snps@ind.names <- gsub(".realigned.bam", "", exon.snps@ind.names)

exon.snps$loc.names <- c(paste0("Locus","_",1:length(exon.snps$loc.names)))

poplist <- as.factor(c(rep("Mar C",45), rep("Mar E", 45), rep("Mar B",33), rep("Quebec", 33), rep("Lilly Canyon", 22), rep("CMA 3B", 34), rep("CMA 6B", 34), 
                   rep("CMA N5440", 21), rep("CMA 10B",6), "Chaleurs", rep("CMA 10B", 7), "Chaleurs", rep("CMA 10B", 7), "Chaleurs", rep("CMA 10B", 7),
                   "Chaleurs",rep("CMA 10B", 7), "Chaleurs", rep("CMA 4", 7), "Chaleurs", rep("CMA 4", 7), "Chaleurs", rep("CMA 4", 7), "Chaleurs",
                   rep("CMA 4", 7), "Chaleurs",rep("CMA 4", 6), rep("Chaleurs",12), rep("CMA 5A", 33), rep("CMA 8A", 32), rep("Lilly Canyon",9), 
                   rep("Chaleurs", 4), rep("CMA 3N200", 10), rep("NAFO 3L", 33), rep("CMA 3D", 35), rep("CMA 10A", 34), rep("St Marys Bay", 27),
                   rep("Fortune Bay", 35), rep("Trinity Bay", 13), rep("West Cape Breton", 8), rep("Bradelle Bank", 44), rep("CMA 3N200", 15),
                   rep("CMA N5440", 3), "St Marys Bay", rep("CMA N5440", 8), rep("West Cape Breton", 4), rep("Trinity Bay", 20), rep("NENSin", 25), 
                   rep("NAFO 4X", 32), rep("NENSout", 35), rep("CMA 12G", 33), rep("CMA 12C", 35), rep("Laurentian Chan", 35), rep("West Cape Breton", 54), 
                   rep("Mar D", 8), rep("West Cape Breton", 65), rep("Mar D", 26)) )
#the first few reps of West Cape Breton are female, and the last 65 are the males

# quick set of summary stats on the data from adegenet
sum.stat <- summary(exon.snps)