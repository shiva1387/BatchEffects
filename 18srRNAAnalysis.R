########################################################
# Analysis of ribosomal sequences-Malaysian algae data #
########################################################
# Author(s): Shiv
# Version: 21042014
# Input: ".txt" file 
# Modified By :Shivshankar Umashankar 
# Functions written here are used for analyzing algal 18srRNA data

#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 

library(gplots)
library(ggplot2)
library(reshape2)
library(vegan)
library(Biostrings)

##### Reading Data

#setwd('../../data/Biochemical_042014/')
Chlorella18s<-readDNAStringSet('Chlorella18s.txt',format="fasta")

#### Calculate distance matrix

distanceMatrix_18s<-stringDist(Chlorella18s, method = "levenshtein", ignoreCase = TRUE, diag = FALSE,
           upper = FALSE, type = "global")