########################
# load mzXML data in R-Arabidopsis TT8 data
# XCMS 
########################
# Author(s): Shiv
# Version: 28082013 
# Input: ".mzXML" files from TT8 agilent data (converted usingMsConverter proteowizard)
# File Location :D:/Shiv_data/Shiv_NUS/PhD_Project/Integrating_Omics/Analysis_R/MetabolomicsSoftwareComparison/Data/TT8_RawData_Metabolomics/
# Output: ".tsv" file containing peaktable
# Modified By :Shivshankar Umashankar 
# Steps followed according to Patti et al VOL.7 NO.3 | 2012 | nature protocols

#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

#############
# User      #
# Specific  #
# Variables #
#############

directory<-"F:/Vinay's Algae Data/Aug2013/Algae/MZXML"; 
# The "PATH", Remember to use "/" and not "/" in the path
# directory contains sub folders each for TT8 and WS 
num_cores<-8; # The number of CPUs available, be aware of memory issue



#############
#   #   #   #
# #   #   # #
#############
ptm<-proc.time() # Process time start

#################
# Load packages #
#################


#################
# Set directory #
#################
setwd(directory) # Sets the path to the folder where the files are located

#################
# Get filenames #
#################
mzxmlfiles<-list.files(directory,pattern=".mzXML",full.names=TRUE) # Get filenames with a given pattern and include subdirectories 


###############################
# Load files and detect peaks # (Be patient)
###############################

set1<-xcmsSet(nSlaves=num_cores,method='centWave',ppm=15,peakwidth=c(5,60), prefilter=c(0,0),snthresh=6) # nSlaves is the number of slaves/cores to be used for parallel peak detection.
#parameters set after discussion with vinay

save(set1,file="D:/Shiv_data/Shiv_NUS/PhD_Project/Integrating_Omics/Analysis_R/MetabolomicsSoftwareComparison/Scripts/Orbitrap/set1_TT8_19042013.Rdata")

save(set1,file="../set1_Algae_29082013.Rdata")

#############################
# Retention time correction #
#############################

set2 <- retcor(set1,method="obiwarp",plottype="deviation")

#save(set2,file="../../../../Scripts/Agilent/set2_TT8_19042013.Rdata")

#################
# Group peaks #
#################
set3 <- group(set2,bw=5,mzwid=0.015,minsamp=1) # bw is the Bandwidth, Setting bw too high will 
#erroneously group peaks of the same m/z value but different rt together. The bw parameter is 
#essentially the maximum time that the features could be separated by.

#########################
# Fill in missing peaks # (Be patient)
#########################
set4 <- fillPeaks(set3) # Retrieves information about the signal intensity where no peaks were detected in some sets

save(set4,file="set4_Algae_29082013_MsLevel1.Rdata")

################# 
# Save peaklist #
################# 
peaklist<-peakTable(set4,filebase="XCMS_Algae_29082013_MsLevel1") # Prints the peaktable as a "tsv" file. "filebase" is the base filename

write.table(peaklist,file="XCMS_Algae_29082013_MsLevel1.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#reporttab <- diffreport(set4, filebase="TT8vsWS")

#################
# Meta XCMS#
################# 

#library(metaXCMS)

##########################
# Report processing time #
##########################
(proc.time()-ptm)[3]# Total processing time (s) 

###########################
##user  system elapsed   ##
#2118.59  213.20 6482.75 ## 
###########################

########################## 