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
nonzero_sigfeat_s_matrix_d12<-day12_nonzero_sigfeat_s[[2]]
nonzero_sigfeat_s_matrix_d12<-nonzero_sigfeat_s_matrix_d12[[4]] #batch corrected data for day12
load('day4_x138_nonzero_sigfeat_s.rda')
nonzero_sigfeat_s_matrix_d4<-day4_nonzero_sigfeat_s[[2]]
nonzero_sigfeat_s_matrix_d4<-nonzero_sigfeat_s_matrix_d4[[7]] #batch corrected data for day4

##################################################
################### Functions ####################
##################################################

### Splitting rownames, rounding mz and rt values
splitMzRt<-function (data_matrix){
  sigfeat_bc_s<-strsplit(rownames(data_matrix), "\\@")
  #9044
  sigfeat_bc_mz_s<-sapply(sigfeat_bc_s , function (x) if(length(x) == 2) x[1] else as.character(NA))
  sigfeat_bc_rt_s<-sapply(sigfeat_bc_s , function (x) if(length(x) == 2) x[2] else as.character(NA))
  sigfeat_bc_mz_s<-round(as.numeric(sigfeat_bc_mz_s),5)
  sigfeat_bc_rt_s<-round(as.numeric(sigfeat_bc_rt_s),0)
  sigfeat_bc_mzrt_s<-paste0(sigfeat_bc_mz_s,'@',sigfeat_bc_rt_s)
  sigfeat_bc<-data.frame(sigfeat_bc_mz_s,sigfeat_bc_rt_s,data_matrix)
  return(sigfeat_bc)
}

### Obtain closest mzrt value
cmpMZ <- function(mz, mzlist){which.min(abs(mz-mzlist))} #obtain the mz which is closest to the actual mass

#extracting sample groups
SampleGroup<-function(data_matrix){

SampleGroups<-sapply(rownames(data_matrix), function(x) paste(strsplit(x,"_")[[1]][1:2],collapse = "_")) 
SampleGroups<-as.vector(SampleGroups)
SampleGroups<-gsub('b1','',SampleGroups)
SampleGroups<-gsub('b2','',SampleGroups)
SampleGroups<-gsub('b3','',SampleGroups)

return(SampleGroups)
}
##################################################
##################################################


setwd('../../data/')
metaboliteData<-read.table("sigfeat_day4_metabolites_identified.txt",header=TRUE,sep='\t')
metaboliteData<-as.data.frame(metaboliteData)
metaboliteData$actualMass<-metaboliteData$Mass+1.005; #To add 1.005 to each compound (as we are searching in the positive mode M+H with 10 ppm error)
metaboliteData$lowerRange<-round(metaboliteData$actualMass-metaboliteData$actualMass*10/10^6,3); #Setting the lower range of 10ppm error to capture metabolites 
metaboliteData$upperRange<-round(metaboliteData$actualMass,3)+0.009; #Setting the upper range of to include variations in third decimal place
metaboliteData<-metaboliteData[order(metaboliteData$actualMass),]

### Mapping metabolites


sigfeat_bc<-splitMzRt(nonzero_sigfeat_s_matrix_d4) #Change to d4 or d12 according to analysis

mappedData<-list()

for (i in 1:nrow(metaboliteData))
{
  # STEP1:Obtain values within the range
  interMappedData<- sigfeat_bc[sigfeat_bc$sigfeat_bc_mz_s >= metaboliteData$lowerRange[i] & 
                                     sigfeat_bc$sigfeat_bc_mz_s <= metaboliteData$upperRange[i],]
  #print(paste(i,length(interMappedData)))
  if(!is.na(interMappedData[1,1]))
  {mappedData[[i]]<-cbind(metaboliteData$actualMass[i],metaboliteData$Compound[i],interMappedData)}
  # STEP2: For values which do not map within the tolerance limits,
  #        obtain the mz which is closest to the value
  if(is.na(interMappedData[1,1])) {
   match_mz<-cmpMZ(sigfeat_bc$sigfeat_bc_mz_s,metaboliteData$actualMass[i])
   interMappedData<-sigfeat_bc[match_mz,]
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
penul_col<-ncol(mappedData1)-1 #as the last column contains p values

mappedMetabolites <- data.table(mappedData1[,c(2,5:penul_col)])
colnames(mappedMetabolites)[1]<-"CompoundName"
mappedMetabolites<-mappedMetabolites[,lapply(.SD, mean),by=CompoundName]
mappedMetabolites<-as.data.frame(t(mappedMetabolites))
colnames(mappedMetabolites)<-as.vector(unlist(mappedMetabolites[1,]))
mappedMetabolites0<-mappedMetabolites[2:nrow(mappedMetabolites),]
mappedMetabolites<-as.data.frame(apply(mappedMetabolites0,2,as.numeric))
rownames(mappedMetabolites)<-rownames(mappedMetabolites0)
rm(mappedMetabolites0)

#Averaging by sample groups

samplegroups<-SampleGroup(mappedMetabolites)

mappedMetabolites$samplegroup<-as.factor(samplegroups)
# correct naming of strains
mappedMetabolites$samplegroup<-as.factor(gsub('D12_14','D12_014',mappedMetabolites$samplegroup))
mappedMetabolites$samplegroup<-as.factor(gsub('D12_84','D12_084',mappedMetabolites$samplegroup))
mappedMetabolites$samplegroup<-as.factor(gsub('D12_87','D12_087',mappedMetabolites$samplegroup))
mappedMetabolites$samplegroup<-as.factor(gsub('D12_94','D12_094',mappedMetabolites$samplegroup))
mappedMetabolites$samplegroup<-as.factor(gsub('D12_94','D12_094',mappedMetabolites$samplegroup))
mappedMetabolites$samplegroup<-as.factor(gsub('D4_14','D4_014',mappedMetabolites$samplegroup))
mappedMetabolites<-mappedMetabolites[order(mappedMetabolites$samplegroup),]
mappedMetabolites_21strains<-mappedMetabolites[which(mappedMetabolites$samplegroup!="D12_184"),]
mappedMetabolites_21strains<-mappedMetabolites[which(mappedMetabolites$samplegroup!="D4_184"),]
dt <- data.table(mappedMetabolites_21strains) # change between mappedMetabolites_21strains and mappedMetabolites to obtain d12_mean_21strains and d12_mean1, respectively.
mappedMetabolites_mean<-dt[,lapply(.SD, mean),by=samplegroup]
mappedMetabolites_mean<-as.data.frame(mappedMetabolites_mean) # 21 or 22 strains
rownames(mappedMetabolites_mean)<-mappedMetabolites_mean[,1]
mappedMetabolites_mean<-mappedMetabolites_mean[,2:ncol(mappedMetabolites_mean)]
# mappedMetabolites_mean_d12_21strains<-mappedMetabolites_mean
#mappedMetabolites_mean_d12<-mappedMetabolites_mean
write.table(round(mappedMetabolites_mean_d12,4),"metabolitesIdentified_day12.txt",quote=FALSE,sep='\t',col.names=NA)
#mappedMetabolites_mean_d4<-mappedMetabolites_mean
write.table(round(mappedMetabolites_mean_d4,4),"metabolitesIdentified_day4.txt",quote=FALSE,sep='\t',col.names=NA)


a<-adonis(mappedMetabolites_mean_d12~ biochemData_d12$biomass+biochemData_d12$lipidProduc+biochemData_d12$totalLipidCon
          +biochemData_d12$totalProtein+biochemData_d12$totalCarbCon, method = "euclidean", perm=999)

a<-adonis(mappedMetabolites_mean_d4~ biochemData_d4$biomass+biochemData_d4$lipidProduc+biochemData_d4$totalLipidCon
          +biochemData_d4$totalProtein+biochemData_d4$totalCarbCon, method = "euclidean", perm=999)

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

#DAY4

# > veg.dist.d4 <- vegdist(as.matrix(mappedMetabolites_mean_d4),method="euclidean")#d4_mean1[,2:10688]
# > veg.dist.d4.biochem <- vegdist(as.matrix(biochemData_d4),method="euclidean")
# > mantel(veg.dist.d4, veg.dist.d4.biochem)
# 
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = veg.dist.d4, ydis = veg.dist.d4.biochem) 
# 
# Mantel statistic r: 0.08611 
#       Significance: 0.241 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.166 0.207 0.251 0.284 
# 
# Based on 999 permutations

# > mantel(veg.dist.d4, distanceMatrix_18s, method="spear")
# 
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d4, ydis = distanceMatrix_18s, method = "spear") 
# 
# Mantel statistic r: 0.1844 
#       Significance: 0.081 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.169 0.218 0.267 0.324 
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

data(QStat)
#QStat.Calculate(QStatX, QStatY, 1000)

PathwayMembership<-as.matrix(c(rep(0,10),rep(1,12)))
testdata<-as.matrix(mappedMetabolites_mean_d12[,1:6])
  
testCalculate<-QStat.Calculate(testdata,PathwayMembership,2)
  
