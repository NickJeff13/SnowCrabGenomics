##########################################################################################
############Redundancy Analyses of Genotypes and Environment############################
## Written by Nick Jeffery in winter 2025

# Load libraries ----------------------------------------------------------
library(adegenet)
library(vegan)
library(dplyr)
library(ggplot2)

#We need 1) genotypes per individual, 2) environmental data per site and individual, and 3) GPS coordinates per individual
# We will be using the exome dataset as the full snp dataset will take a long time


# Read in plink genotypes -------------------------------------------------


