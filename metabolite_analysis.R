###########################################
# Data Analysis in R-Malaysian algae data #
###########################################
# Author(s): Shiv
# Version: 18052014
# Input: ".txt" file 
# Modified By :Shivshankar Umashankar 
# Functions written here are used for analyzing singletons and non-singletons
# Non-singlteon mz is later used for metabolite identification
# All mz features are mapped using PCDL manager to Metlin metabolites and Metacyc chlorella database
# Patways are identified by mapping them onto kegg chlorella database

#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 

library(data.table)
library(multtest)
library(ggplot2)
library(Biostrings)
library(vegan)

### Reading Data-svd and differential metabolites

setwd('../../AlgaeData/results.from.scelse.cluster.211213/')

load("svd_day4_x138_nonzero.rda")
batch_corrected_mat_d4<-svd_day4_nonzero[[7]]
#write.table(batch_corrected_mat_d4,"batch_corrected_mat_d4.txt",quote=FALSE,sep="\t",col.names=NA)
load("day4_x138_nonzero_sigfeat_s.rda") #loading strain data
day4_nonzero_sigfeat_s_matrix<-day4_nonzero_sigfeat_s[[2]]
batch_corrected_sig_mat_d4<-day4_nonzero_sigfeat_s_matrix[[7]] #as removing after 7 rounds of pc removal, runday effect is removed

## Loading significant features due to RunDay
# load("day4_x138_nonzero_sigfeat_r.rda") #loading runday data
# day4_nonzero_sigfeat_r_matrix<-day4_nonzero_sigfeat_r[[2]]
# batch_corrected_sig_mat_d4_run<-day4_nonzero_sigfeat_r_matrix[[7]] #as removing after 7 rounds of pc removal, runday effect is removed


load("svd_day12_x138_nonzero.rda")
batch_corrected_mat_d12<-svd_day12_nonzero[[4]]
#write.table(batch_corrected_mat_d12,"batch_corrected_mat_d12.txt",quote=FALSE,sep="\t",col.names=NA)
load("day12_x138_nonzero_sigfeat_s.rda") #loading strain data
day12_nonzero_sigfeat_s_matrix<-day12_nonzero_sigfeat_s[[2]]
batch_corrected_sig_mat_d12<-day12_nonzero_sigfeat_s_matrix[[4]] #as removing after 4 rounds of pc removal, runday effect is removed

## Loading significant features due to RunDay
# load("day12_x138_nonzero_sigfeat_r.rda") #loading runday data
# day12_nonzero_sigfeat_r_matrix<-day12_nonzero_sigfeat_r[[2]]
# batch_corrected_sig_mat_d12_run<-day12_nonzero_sigfeat_r_matrix[[4]] #as removing after 4 rounds of pc removal, runday effect is removed


#Ensure all strain ids are in the same format
colnames(batch_corrected_sig_mat_d12)<-as.character(gsub('D12_14','D12_014',colnames(batch_corrected_sig_mat_d12)))
colnames(batch_corrected_sig_mat_d12)<-as.character(gsub('D12_84','D12_084',colnames(batch_corrected_sig_mat_d12)))
colnames(batch_corrected_sig_mat_d12)<-as.character(gsub('D12_87','D12_087',colnames(batch_corrected_sig_mat_d12)))
colnames(batch_corrected_sig_mat_d12)<-as.character(gsub('D12_94','D12_094',colnames(batch_corrected_sig_mat_d12)))
colnames(batch_corrected_sig_mat_d12)<-as.character(gsub('D12_94','D12_094',colnames(batch_corrected_sig_mat_d12)))
colnames(batch_corrected_sig_mat_d4)<-as.character(gsub('D4_14','D4_014',colnames(batch_corrected_sig_mat_d4)))

colnames(batch_corrected_mat_d12)<-as.character(gsub('D12_14','D12_014',colnames(batch_corrected_mat_d12)))
colnames(batch_corrected_mat_d12)<-as.character(gsub('D12_84','D12_084',colnames(batch_corrected_mat_d12)))
colnames(batch_corrected_mat_d12)<-as.character(gsub('D12_87','D12_087',colnames(batch_corrected_mat_d12)))
colnames(batch_corrected_mat_d12)<-as.character(gsub('D12_94','D12_094',colnames(batch_corrected_mat_d12)))
colnames(batch_corrected_mat_d12)<-as.character(gsub('D12_94','D12_094',colnames(batch_corrected_mat_d12)))
colnames(batch_corrected_mat_d4)<-as.character(gsub('D4_14','D4_014',colnames(batch_corrected_mat_d4)))


# Reading singleton data

setwd('../../data/')
isotopeData<-read.table('algae_4setblanks_021213_mzrt_camera.txt',header=TRUE,sep='\t')
#isotopeData$mzrt<-paste0(round(as.numeric(isotopeData$mz),5),'@',isotopeData$rt)
isotopeData$mz<-round(as.numeric(isotopeData$mz),5)
isotopeData$rt<-round(as.numeric(isotopeData$rt),0)
isotopeData$mzrt<-paste0(isotopeData$mz,'@',isotopeData$rt)
nonSingleton<-duplicated(isotopeData$pcgroup) | duplicated(isotopeData$pcgroup, fromLast = TRUE)
mzrt_nonSingleton<-isotopeData[which(nonSingleton==TRUE),]

##### Biochemical data
setwd('Biochemical_042014/')
biochemData0<-read.table("biochemicalData.txt",header=TRUE,sep='\t')
biochemData<-as.data.frame(biochemData0[,3:9])
rownames(biochemData)<-biochemData0[,2]
biochemData_d4<-biochemData[grep('D4', rownames(biochemData)),]
biochemData_d4<-biochemData_d4[order(rownames(biochemData_d4)),]
biochemData_d12<-biochemData[grep('D12', rownames(biochemData)),] 
biochemData_d12<-biochemData_d12[order(rownames(biochemData_d12)),]

#gsub ("a|b|c","",rownames(biochemData_d4))

biochemData_d4$samplegroup<-as.factor(gsub ("a|b|c","",rownames(biochemData_d4)))
biochemData_d4 <- data.table(biochemData_d4)
biochemData_d4<-biochemData_d4[,lapply(.SD, mean),by=samplegroup]
biochemData_d4<-as.data.frame(biochemData_d4)
rownames(biochemData_d4)<-biochemData_d4[,1]
biochemData_d4<-biochemData_d4[,2:8]
biochemData_d4<-as.data.frame(scale(biochemData_d4,center=TRUE, scale=TRUE))

biochemData_d12$samplegroup<-as.factor(gsub ("a|b|c","",rownames(biochemData_d12)))
biochemData_d12 <- data.table(biochemData_d12)
biochemData_d12<-biochemData_d12[,lapply(.SD, mean),by=samplegroup]
biochemData_d12<-as.data.frame(biochemData_d12)
rownames(biochemData_d12)<-biochemData_d12[,1]
biochemData_d12<-biochemData_d12[,2:7] #There is no growth rate associated with day 12
biochemData_d12<-as.data.frame(scale(biochemData_d12,center=TRUE, scale=TRUE))

### Plotting distance matrix
d <- dist(biochemData_d4, method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
png("d4_biochemData_d4.png")
plot (fit,main="d14_biochemData_d4")
dev.off()

##### 18s rRNA data
Chlorella18s<-readDNAStringSet('Chlorella18s.txt',format="fasta")
distanceMatrix_18s<-stringDist(Chlorella18s, method = "levenshtein", ignoreCase = TRUE, diag = FALSE,
                               upper = FALSE, type = "global")
# fit <- hclust(distanceMatrix_18s, method="average") 
# plot (fit)
# dev.off()

#######################################################
################ Functions ############################
#######################################################

combineRoundMzMatrix<-function(batch_corrected_mat)
{
  mat_bc_s<-strsplit(rownames(batch_corrected_mat), "\\@")
  mat_bc_mz_s<-sapply(mat_bc_s , function (x) if(length(x) == 2) x[1] else as.character(NA))
  mat_bc_rt_s<-sapply(mat_bc_s , function (x) if(length(x) == 2) x[2] else as.character(NA))
  mat_bc_mz_s<-round(as.numeric(mat_bc_mz_s),5)
  mat_bc_rt_s<-round(as.numeric(mat_bc_rt_s),0)
  mat_bc_mzrt_s<-paste0(mat_bc_mz_s,'@',mat_bc_rt_s)
  mat_bc<-data.frame(mat_bc_mzrt_s,mat_bc_mz_s,mat_bc_rt_s)
}

#functions to get matching mz and rt
mz_tolerance=10 #ppm
rt_tolerance=5 #secs

cmpMZ <- function(mz, mzlist, mz_cutoff){abs(mz-mzlist) <= mz_cutoff}
#if there is more than one match, change function to pick the closest!
cmpRT <- function(ret, retlist,ret_cutoff){abs(ret-retlist) <= ret_cutoff}

## Matching and extracting non-singletons from dataset

matchAndExtract<-function(mzrtMatrix)
  {
    for(i in 1:nrow(mzrtMatrix)) #getting total length
    {
      ppm_tolerance<-(trunc(mzrtMatrix[i,2])*mz_tolerance)/10^6 #Using mz
      match_mz<-cmpMZ(mzrtMatrix[i,2],mzrt_nonSingleton$mz,ppm_tolerance)
      match_mz_ind<-which(match_mz) #obtaining indices only for the significant features
      match_ret<-cmpRT(mzrtMatrix[i,3],mzrt_nonSingleton$rt[match_mz_ind],rt_tolerance)#Using rt
      mz_rt_match<-match_mz_ind[which(match_ret)]
      matched_value<-mzrt_nonSingleton[mz_rt_match,"mzrt"];
      matched_value=ifelse(exists("matched_value"),matched_value)
      #print(c(i,sigfeat_day4_bc_mz_s[i],sigfeat_day4_bc_rt_s[i],mz_rt_match,matched_value))
      mzrtMatrix[i,4]<-matched_value
    }
  bc_nonsingletons<-mzrtMatrix[!is.na(mzrtMatrix$V4),]
  #writing mz rt to a file
  return (bc_nonsingletons)
}

## Functions associated with mapped metabolites ##

### Obtain closest mzrt value
cmpMZmin <- function(mz, mzlist){which.min(abs(mz-mzlist))} #obtain the mz which is closest to the actual mass

### Splitting rownames, rounding mz and rt values
splitMzRt<-function (data_matrix){
  sigfeat_bc_s<-strsplit(rownames(data_matrix), "\\@")
  sigfeat_bc_mz_s<-sapply(sigfeat_bc_s , function (x) if(length(x) == 2) x[1] else as.character(NA))
  sigfeat_bc_rt_s<-sapply(sigfeat_bc_s , function (x) if(length(x) == 2) x[2] else as.character(NA))
  sigfeat_bc_mz_s<-round(as.numeric(sigfeat_bc_mz_s),5)
  sigfeat_bc_rt_s<-round(as.numeric(sigfeat_bc_rt_s),0)
  sigfeat_bc_mzrt_s<-paste0(sigfeat_bc_mz_s,'@',sigfeat_bc_rt_s)
  sigfeat_bc<-data.frame(sigfeat_bc_mz_s,sigfeat_bc_rt_s,data_matrix)
  return(sigfeat_bc)
}

#extracting sample groups
SampleGroup<-function(data_matrix){
  SampleGroups<-sapply(rownames(data_matrix), function(x) paste(strsplit(x,"_")[[1]][1:2],collapse = "_")) 
  SampleGroups<-as.vector(SampleGroups)
  SampleGroups<-gsub('b1','',SampleGroups)
  SampleGroups<-gsub('b2','',SampleGroups)
  SampleGroups<-gsub('b3','',SampleGroups)
  
  return(SampleGroups)
}

#format metabolite data-add masss
formatMetaboliteData<-function(metaboCpdMassmatrix) { # metaboCpdMassmatrix should contain 2 columns, Mass and Cpd id  
  metaboCpdMassmatrix$actualMass<-metaboCpdMassmatrix$Mass+1.005; #To add 1.005 to each compound (as we are searching in the positive mode M+H with 10 ppm error)
  metaboCpdMassmatrix$lowerRange<-round(metaboCpdMassmatrix$actualMass-metaboCpdMassmatrix$actualMass*10/10^6,3); #Setting the lower range of 10ppm error to capture metabolites 
  metaboCpdMassmatrix$upperRange<-round(metaboCpdMassmatrix$actualMass+metaboCpdMassmatrix$actualMass*10/10^6,3); #Setting the upper range of to include variations in third decimal place
  #metaboCpdMassmatrix$upperRange<-round(metaboCpdMassmatrix$actualMass,3)+0.009; #Setting the upper range of to include variations in third decimal place
  metaboCpdMassmatrix<-metaboCpdMassmatrix[order(metaboCpdMassmatrix$actualMass),]
  return(metaboCpdMassmatrix)
}

#mapping metabolites

mappingMetabolites<-function(massSpecDataMatrix,metaboliteData){  #massSpecDataMatrix is the mz data, metaboliteDat contains cpd id and mass
mappedData<-list()
sigfeat_bc<-splitMzRt(massSpecDataMatrix) #Change to d4 or d12 according to analysis

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
    match_mz<-cmpMZmin(sigfeat_bc$sigfeat_bc_mz_s,metaboliteData$actualMass[i])
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
return (mappedData1)
}

###############################################################
############## Extracting non-singletons ######################

##### DAY12
## significant features
sigfeat_day12_bc_s<-combineRoundMzMatrix(batch_corrected_sig_mat_d12)
significant_nonSingleton_day12<-matchAndExtract(sigfeat_day12_bc_s)
write.table(significant_nonSingleton_day12,"day12_differential_nonsingletonList.txt",quote=FALSE,row.names=FALSE,sep="\t")
## Full data
day12_bc_s<-combineRoundMzMatrix(batch_corrected_mat_d12)
nonSingleton_day12<-matchAndExtract(day12_bc_s)
write.table(nonSingleton_day12,"day12_nonsingletonList.txt",quote=FALSE,row.names=FALSE,sep="\t")

##### DAY4
## significant features
sigfeat_day4_bc_s<-combineRoundMzMatrix(batch_corrected_sig_mat_d4)
significant_nonSingleton_day4<-matchAndExtract(sigfeat_day4_bc_s)
write.table(significant_nonSingleton_day4,"day4_differential_nonsingletonList.txt",quote=FALSE,row.names=FALSE,sep="\t")
## Full data
day4_bc_s<-combineRoundMzMatrix(batch_corrected_mat_d4)
nonSingleton_day4<-matchAndExtract(day4_bc_s)
write.table(nonSingleton_day4,"day4_nonsingletonList.txt",quote=FALSE,row.names=FALSE,sep="\t")

##### Biochemical data comparisons
## significant features
significant_features<-combineRoundMzMatrix(metab_sig)
significant_features_nonSingleton<-matchAndExtract(significant_features)
write.table(significant_features_nonSingleton,"day12_differential_bestVsWorst.txt",quote=FALSE,row.names=FALSE,sep="\t")

################################
# Metabolite pathway analysis  #
################################

metaboliteDataset<-read.table("day12_metabolitesId.txt",header=TRUE,sep='\t')
metaboliteDataset<-as.data.frame(metaboliteDataset)
metaboliteDataset<-formatMetaboliteData(metaboliteDataset)

mappedMetaboliteData<-mappingMetabolites(batch_corrected_mat_d12,metaboliteDataset)

write.table(mappedMetaboliteData,"day12_metabolites_mzrt_compoundId.txt",quote=FALSE,sep='\t')

### Averaging values of multiple hits

#Averaging by compound ids
penul_col<-ncol(mappedMetaboliteData)-1 # -1 to be used only for differential metabolites
#as the last column contains p values

mappedMetabolites <- data.table(mappedMetaboliteData[,c(2,5:penul_col)]) #columns 1,3,4 are meta data
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
#mappedMetabolites_21strains<-mappedMetabolites[which(mappedMetabolites$samplegroup!="D12_184"),]
#mappedMetabolites_21strains<-mappedMetabolites_21strains[which(mappedMetabolites_21strains$samplegroup!="D4_184"),]
dt <- data.table(mappedMetabolites) # change between mappedMetabolites_21strains and mappedMetabolites to obtain d12_mean_21strains and d12_mean1, respectively.
mappedMetabolites_mean<-dt[,lapply(.SD, mean),by=samplegroup]
mappedMetabolites_mean<-as.data.frame(mappedMetabolites_mean) # 21 or 22 strains
rownames(mappedMetabolites_mean)<-mappedMetabolites_mean[,1]
mappedMetabolites_mean<-mappedMetabolites_mean[,2:ncol(mappedMetabolites_mean)]
# mappedMetabolites_mean_d12_21strains<-mappedMetabolites_mean
#mappedMetabolites_mean_d12<-mappedMetabolites_mean
#write.table(round(mappedMetabolites_mean,4),"metabolitesIdentified_day12.txt",quote=FALSE,sep='\t',col.names=NA)
#mappedMetabolites_mean_d4<-mappedMetabolites_mean
#write.table(round(mappedMetabolites_mean,4),"metabolitesIdentified_day4.txt",quote=FALSE,sep='\t',col.names=NA)

### Plotting distance matrix

d <- dist(mappedMetabolites_mean, method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
png("d4_totalMetabolites.png")
plot (fit,main="d4_totalMetabolites_21strains")
dev.off()


##### Analysis of meta data

## Comparing day 12 and day 4 biochem matrix

# mantel(veg.dist.d12.biochem,veg.dist.d4.biochem,method="spear")
# 
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d12.biochem, ydis = veg.dist.d4.biochem,      method = "spear") 
# 
# Mantel statistic r: -0.00405 
# Significance: 0.482 
# 
# Upper quantiles of permutations (null model):
# 90%   95% 97.5%   99% 
# 0.114 0.150 0.184 0.220 
# 
# Based on 999 permutations

## Comparing day12 and day4 metabolite matrices

# mappedMetabolites_mean_d12<-read.table("metabolitesIdentified_day12.txt",header=TRUE,sep='\t',row.names=1)
# mappedMetabolites_mean_d4<-read.table("metabolitesIdentified_day4.txt",header=TRUE,sep='\t',row.names=1)

veg.dist.d12 <- vegdist(as.matrix(mappedMetabolites_mean_d12),method="euclidean") 
veg.dist.d4 <- vegdist(as.matrix(mappedMetabolites_mean_d4),method="euclidean") 
mantel(veg.dist.d12, veg.dist.d4, method="spear")

## Full dataset metabolites (extracted from the full data batch_corrected_mat)

# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d12, ydis = veg.dist.d4, method = "spear") 
# 
# Mantel statistic r: 0.07954 
#       Significance: 0.297 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.187 0.227 0.262 0.319 
# 
# Based on 999 permutations

## Differential metabolites (obtained from the above files, metabolitesIdentified_day12.txt)
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d12, ydis = veg.dist.d4, method = "spear") 
# 
# Mantel statistic r: 0.02115 
# Significance: 0.43 
# 
# Upper quantiles of permutations (null model):
# 90%   95% 97.5%   99% 
# 0.166 0.219 0.253 0.289 
#
# Based on 999 permutations

########### DAY 12

##### FullData
#####
### Permanova

adonis(mappedMetabolites_mean_d12~ biochemData_d12$biomass+biochemData_d12$lipidProduc+biochemData_d12$totalLipidCon
          +biochemData_d12$totalProtein+biochemData_d12$totalCarbCon, method = "euclidean", perm=999)
# Call:
#   adonis(formula = mappedMetabolites_mean_d12 ~ biochemData_d12$biomass +      biochemData_d12$lipidProduc + biochemData_d12$totalLipidCon +      biochemData_d12$totalProtein + biochemData_d12$totalCarbCon,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# biochemData_d12$biomass        1    129.02 129.021 1.76456 0.07947  0.019 *
# biochemData_d12$lipidProduc    1     85.53  85.532 1.16979 0.05268  0.257  
# biochemData_d12$totalLipidCon  1     86.68  86.681 1.18550 0.05339  0.260  
# biochemData_d12$totalProtein   1     80.98  80.984 1.10758 0.04988  0.322  
# biochemData_d12$totalCarbCon   1     71.40  71.39Bio5 0.97644 0.04398  0.497  
# Residuals                     16   1169.88  73.118         0.72060         
# Total                         21   1623.50                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

### Mantel test
veg.dist.d12 <- vegdist(as.matrix(mappedMetabolites_mean_d12),method="euclidean") #d12_mean1[,2:10688]
veg.dist.d12.biochem <- vegdist(as.matrix(biochemData_d12),method="euclidean")
mantel(veg.dist.d12, veg.dist.d12.biochem, method="spear")

# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d12, ydis = veg.dist.d12.biochem, method = "spear") 
# 
# Mantel statistic r: 0.2632 
#       Significance: 0.013 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.150 0.197 0.235 0.288 
# 
# Based on 999 permutations

veg.dist.d12 <- vegdist(as.matrix(mappedMetabolites_mean_d12_21strains),method="euclidean") #d12_mean1[,2:10688]
mantel(veg.dist.d12, distanceMatrix_18s, method="spear")

# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d12, ydis = distanceMatrix_18s, method = "spear") 
# 
# Mantel statistic r: -0.0764 
#       Significance: 0.747 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.166 0.210 0.262 0.310 
# 
# Based on 999 permutations


##### Differential metabolites
#####

# Call:
#   adonis(formula = mappedMetabolites_mean_d12 ~ biochemData_d12$biomass +      biochemData_d12$lipidProduc + biochemData_d12$totalLipidCon +      biochemData_d12$totalProtein + biochemData_d12$totalCarbCon,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# biochemData_d12$biomass        1    118.04 118.039  1.8958 0.08408  0.015 *
#   biochemData_d12$lipidProduc    1     79.25  79.246  1.2727 0.05645  0.220  
# biochemData_d12$totalLipidCon  1     75.40  75.398  1.2109 0.05371  0.224  
# biochemData_d12$totalProtein   1     68.51  68.510  1.1003 0.04880  0.369  
# biochemData_d12$totalCarbCon   1     66.50  66.499  1.0680 0.04737  0.406  
# Residuals                     16    996.23  62.264         0.70960         
# Total                         21   1403.92                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d12, ydis = veg.dist.d12.biochem, method = "spear") 
# 
# Mantel statistic r: 0.2485 
#       Significance: 0.016 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.144 0.185 0.223 0.260 
# 
# Based on 999 permutations

veg.dist.d12 <- vegdist(as.matrix(mappedMetabolites_mean_d12_21strains),method="euclidean") #d12_mean1[,2:10688]
mantel(veg.dist.d12, distanceMatrix_18s, method="spear")
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d12, ydis = distanceMatrix_18s, method = "spear") 
# 
# Mantel statistic r: -0.1124 
#       Significance: 0.814 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.166 0.210 0.255 0.306 
# 
# Based on 999 permutations


########### DAY 4

##### FullData
#####

### Permanova
adonis(mappedMetabolites_mean_d4~ biochemData_d4$biomass+biochemData_d4$lipidProduc+biochemData_d4$totalLipidCon
       +biochemData_d4$totalProtein+biochemData_d4$totalCarbCon +biochemData_d4$growthRate, method = "euclidean", perm=999)
# Call:
#   adonis(formula = mappedMetabolites_mean_d4 ~ biochemData_d4$biomass +      biochemData_d4$lipidProduc + biochemData_d4$totalLipidCon +      biochemData_d4$totalProtein + biochemData_d4$totalCarbCon +      biochemData_d4$growthRate, permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# biochemData_d4$biomass        1     65.97  65.968 1.33942 0.06524  0.135
# biochemData_d4$lipidProduc    1     49.12  49.122 0.99738 0.04858  0.481
# biochemData_d4$totalLipidCon  1     31.44  31.441 0.63838 0.03109  0.945
# biochemData_d4$totalProtein   1     33.04  33.038 0.67081 0.03267  0.900
# biochemData_d4$totalCarbCon   1     55.61  55.606 1.12905 0.05499  0.301
# biochemData_d4$growthRate     1     37.28  37.285 0.75703 0.03687  0.795
# Residuals                    15    738.76  49.251         0.73056       
# Total                        21   1011.22                 1.00000  

### Mantel test
veg.dist.d4 <- vegdist(as.matrix(mappedMetabolites_mean_d4),method="euclidean") #d12_mean1[,2:10688]
veg.dist.d4.biochem <- vegdist(as.matrix(biochemData_d4),method="euclidean")
mantel(veg.dist.d4, veg.dist.d4.biochem, method="spear")
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d4, ydis = veg.dist.d4.biochem, method = "spear") 
# 
# Mantel statistic r: -0.01178 
#       Significance: 0.559 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.148 0.183 0.233 0.262 
# 
# Based on 999 permutations

veg.dist.d4 <- vegdist(as.matrix(mappedMetabolites_mean),method="euclidean") #d12_mean1[,2:10688]
mantel(veg.dist.d4, distanceMatrix_18s, method="spear")
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d4, ydis = distanceMatrix_18s, method = "spear") 
# 
# Mantel statistic r: 0.1081 
#       Significance: 0.225 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.174 0.233 0.276 0.324 
# 
# Based on 999 permutations

##### Differential metabolites
#####

### Permanova
adonis(mappedMetabolites_mean_d4~ biochemData_d4$biomass+biochemData_d4$lipidProduc+biochemData_d4$totalLipidCon
          +biochemData_d4$totalProtein+biochemData_d4$totalCarbCon+biochemData_d4$growthRate, method = "euclidean", perm=999)
# Call:
#   adonis(formula = mappedMetabolites_mean_d4 ~ biochemData_d4$biomass +      biochemData_d4$lipidProduc + biochemData_d4$totalLipidCon +      biochemData_d4$totalProtein + biochemData_d4$totalCarbCon +      biochemData_d4$growthRate, permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# biochemData_d4$biomass        1     48.08  48.085 1.31000 0.06412  0.174
# biochemData_d4$lipidProduc    1     34.04  34.043 0.92746 0.04540  0.515
# biochemData_d4$totalLipidCon  1     23.31  23.312 0.63511 0.03109  0.883
# biochemData_d4$totalProtein   1     23.06  23.059 0.62821 0.03075  0.900
# biochemData_d4$totalCarbCon   1     43.48  43.476 1.18446 0.05798  0.262
# biochemData_d4$growthRate     1     27.32  27.319 0.74428 0.03643  0.751
# Residuals                    15    550.59  36.706         0.73423       
# Total                        21    749.88                 1.00000 

### Mantel test
veg.dist.d4 <- vegdist(as.matrix(mappedMetabolites_mean_d4),method="euclidean") #d12_mean1[,2:10688]
veg.dist.d4.biochem <- vegdist(as.matrix(biochemData_d4),method="euclidean")
mantel(veg.dist.d4, veg.dist.d4.biochem, method="spear")
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d4, ydis = veg.dist.d4.biochem, method = "spear") 
# 
# Mantel statistic r: 0.02485 
#       Significance: 0.398 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.154 0.195 0.233 0.259 
# 
# Based on 999 permutations

veg.dist.d4 <- vegdist(as.matrix(mappedMetabolites_mean),method="euclidean") #d12_mean1[,2:10688]
mantel(veg.dist.d4, distanceMatrix_18s, method="spear")
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d4, ydis = distanceMatrix_18s, method = "spear") 
# 
# Mantel statistic r: 0.1805 
#       Significance: 0.102 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.177 0.237 0.277 0.322 
# 
# Based on 999 permutations

