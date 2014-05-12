#########################################################
# Metabolite pathway analysis in R-Malaysian algae data #
#########################################################
# Author(s): Shiv
# Version: 13052014
# Input: ".txt" file 
# Modified By :Shivshankar Umashankar 
# Functions written here are used for analyzing pathway data by performing global test

#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 

library(globaltest)

### Reading Data

setwd('../AlgaeData/results.from.scelse.cluster.211213/')
load('day12_x138_nonzero_sigfeat_s.rda')
day12_nonzero_sigfeat_s_matrix<-day12_nonzero_sigfeat_s[[2]]
load('day4_x138_nonzero_sigfeat_s.rda')


### Splitting rownames, rounding mz and rt values

sigfeat_day12_bc_s<-strsplit(rownames(day12_nonzero_sigfeat_s_matrix[[4]]), "\\@")
#9044
sigfeat_day12_bc_mz_s<-sapply(sigfeat_day12_bc_s , function (x) if(length(x) == 2) x[1] else as.character(NA))
sigfeat_day12_bc_rt_s<-sapply(sigfeat_day12_bc_s , function (x) if(length(x) == 2) x[2] else as.character(NA))
sigfeat_day12_bc_mz_s<-round(as.numeric(sigfeat_day12_bc_mz_s),5)
sigfeat_day12_bc_rt_s<-round(as.numeric(sigfeat_day12_bc_rt_s),0)
sigfeat_day12_bc_mzrt_s<-paste0(sigfeat_day12_bc_mz_s,'@',sigfeat_day12_bc_rt_s)
sigfeat_day12_bc<-data.frame(sigfeat_day12_bc_mz_s,sigfeat_day12_bc_rt_s,day12_nonzero_sigfeat_s_matrix[[4]])

######



setwd('../../data/')
metaboliteData<-read.table("sigfeat_day12_metabolites_identified.txt",header=TRUE,sep='\t')
metaboliteData<-as.data.frame(metaboliteData)
metaboliteData$actualMass<-metaboliteData$Mass+1.005; #To add 1.005 to each compound (as we are searching in the positive mode M+H with 10 ppm error)
metaboliteData$lowerRange<-round(metaboliteData$actualMass-metaboliteData$actualMass*10/10^6,3); #Setting the lower range of 10ppm error to capture metabolites 
metaboliteData$upperRange<-round(metaboliteData$actualMass,3)+0.009; #Setting the upper range of to include variations in third decimal place
metaboliteData<-metaboliteData[order(metaboliteData$actualMass),]

### Mapping metabolites
mappedData<-list()
naValues<-"NA"
cmpMZ <- function(mz, mzlist){which.min(abs(mz-mzlist))} #obtain the mz which is closest to the actual mass

for (i in 1:nrow(metaboliteData))
{
  # STEP1:Obtain values within the range
  interMappedData<- sigfeat_day12_bc[sigfeat_day12_bc$sigfeat_day12_bc_mz_s >= metaboliteData$lowerRange[i] & 
                                     sigfeat_day12_bc$sigfeat_day12_bc_mz_s <= metaboliteData$upperRange[i],]
  #print(paste(i,length(interMappedData)))
  if(!is.na(interMappedData[1,1]))
  {mappedData[[i]]<-cbind(metaboliteData$actualMass[i],metaboliteData$Compound[i],interMappedData)}
  # STEP2: For values which do not map within the tolerance limits,
  #        obtain the mz which is closest to the value
  if(is.na(interMappedData[1,1])) {
   match_mz<-cmpMZ(sigfeat_day12_bc$sigfeat_day12_bc_mz_s,metaboliteData$actualMass[i])
   interMappedData<-sigfeat_day12_bc[match_mz,]
   }
  if(!is.na(interMappedData[1,1]))
  {mappedData[[i]]<-cbind(metaboliteData$actualMass[i],metaboliteData$Compound[i],interMappedData)}
  # STEP3: Check if any mz is left out
  if(is.na(interMappedData[1,1]))
  { print(paste(metaboliteData$lowerRange[i],metaboliteData$upperRange[i]))}
}

mappedData1<-do.call(rbind,mappedData) #Complete list of mapped metabolites and corresponding features
# Will contain mmultiple hits to each metabolite id

### Averaging values of multiple hits

#Averaging by compound ids
mappedMetabolites <- data.table(mappedData1[,c(2,5:122)])
colnames(mappedMetabolites)[1]<-"CompoundName"
mappedMetabolites<-mappedMetabolites[,lapply(.SD, mean),by=CompoundName]
mappedMetabolites<-as.data.frame(t(mappedMetabolites))
colnames(mappedMetabolites)<-as.vector(unlist(mappedMetabolites[1,]))
mappedMetabolites0<-mappedMetabolites[2:119,]
mappedMetabolites<-as.data.frame(apply(mappedMetabolites0,2,as.numeric))
rownames(mappedMetabolites)<-rownames(mappedMetabolites0)
rm(mappedMetabolites0)

#Averaging by sample groups
#extracting sample groups
SampleGroups<-sapply(rownames(mappedMetabolites), function(x) paste(strsplit(x,"_")[[1]][1:2],collapse = "_")) #change d12_test or d4_test accordingly
SampleGroups<-as.vector(SampleGroups)

SampleGroups<-gsub('b1','',SampleGroups)
SampleGroups<-gsub('b2','',SampleGroups)
SampleGroups<-gsub('b3','',SampleGroups)

mappedMetabolites$samplegroup<-as.factor(SampleGroups)
# correct naming of strains
mappedMetabolites$samplegroup<-as.factor(gsub('D12_14','D12_014',mappedMetabolites$samplegroup))
mappedMetabolites$samplegroup<-as.factor(gsub('D12_84','D12_084',mappedMetabolites$samplegroup))
mappedMetabolites$samplegroup<-as.factor(gsub('D12_87','D12_087',mappedMetabolites$samplegroup))
mappedMetabolites$samplegroup<-as.factor(gsub('D12_94','D12_094',mappedMetabolites$samplegroup))
mappedMetabolites<-mappedMetabolites[order(mappedMetabolites$samplegroup),]
mappedMetabolites_21strains<-mappedMetabolites[which(mappedMetabolites$samplegroup!="D12_184"),]
dt <- data.table(mappedMetabolites) # change between mappedMetabolites_21strains and mappedMetabolites to obtain d12_mean_21strains and d12_mean1, respectively.
mappedMetabolites_mean<-dt[,lapply(.SD, mean),by=samplegroup]
mappedMetabolites_mean<-as.data.frame(mappedMetabolites_mean) # 21 or 22 strains
rownames(mappedMetabolites_mean)<-mappedMetabolites_mean[,1]
mappedMetabolites_mean<-mappedMetabolites_mean[,2:655]
# mappedMetabolites_mean_d12_21strains<-mappedMetabolites_mean
# mappedMetabolites_mean_d12<-mappedMetabolites_mean

##### Mantel test

# > veg.dist.d12 <- vegdist(as.matrix(mappedMetabolites_mean),method="euclidean") #d12_mean1[,2:10688]
# > veg.dist.d12.biochem <- vegdist(as.matrix(biochemData_d12),method="euclidean")
# > mantel(veg.dist.d12, veg.dist.d12.biochem)
# 
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = veg.dist.d12, ydis = veg.dist.d12.biochem) 
# 
# Mantel statistic r: 0.2658 
#       Significance: 0.012 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.138 0.189 0.220 0.267 
# 
# Based on 999 permutations

# > mantel(veg.dist.d12, distanceMatrix_18s, method="spear")
# 
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d12, ydis = distanceMatrix_18s, method = "spear") 
# 
# Mantel statistic r: -0.1107 
#       Significance: 0.786 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.160 0.207 0.252 0.291 
# 
# Based on 999 permutations


### Global test analysis

# library(globaltest)
# gt(Y~1, Y~A+B+C, data = X)
# 
# data(exampleX)      # Expression data (40 samples; 1000 genes)
# data(exampleY)      # Clinical outcome for the 40 samples
# pathway1 <- 1:25    # A pathway contains genes 1 to 25
# pathway2 <- 26:50   # another pathway
# gt <- globaltest(exampleX, exampleY, list(pathway1,pathway2))
# gt
