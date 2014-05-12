###########################################
# Data Analysis in R-Malaysian algae data #
###########################################
# Author(s): Shiv
# Version: 07052014
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

### Reading Data-svd and differential metabolites

setwd('../../AlgaeData/results.from.scelse.cluster.211213/')

load("svd_day4_x138_nonzero.rda")
batch_corrected_mat_d4<-svd_day4_nonzero[[7]]
load("day4_x138_nonzero_sigfeat_s.rda") #loading strain data
day4_nonzero_sigfeat_s_matrix<-day4_nonzero_sigfeat_s[[2]]
batch_corrected_sig_mat_d4<-day4_nonzero_sigfeat_s_matrix[[7]] #as removing after 7 rounds of pc removal, runday effect is removed


load("svd_day12_x138_nonzero.rda")
batch_corrected_mat_d12<-svd_day12_nonzero[[4]]
load("day12_x138_nonzero_sigfeat_s.rda") #loading strain data
day12_nonzero_sigfeat_s_matrix<-day12_nonzero_sigfeat_s[[2]]
batch_corrected_sig_mat_d12<-day12_nonzero_sigfeat_s_matrix[[4]] #as removing after 4 rounds of pc removal, runday effect is removed

# Reading singleton data

setwd('../../data/')
isotopeData<-read.table('algae_4setblanks_021213_mzrt_camera.txt',header=TRUE,sep='\t')
#isotopeData$mzrt<-paste0(round(as.numeric(isotopeData$mz),5),'@',isotopeData$rt)
isotopeData$mz<-round(as.numeric(isotopeData$mz),5)
isotopeData$rt<-round(as.numeric(isotopeData$rt),0)
isotopeData$mzrt<-paste0(isotopeData$mz,'@',isotopeData$rt)
nonSingleton<-duplicated(isotopeData$pcgroup) | duplicated(isotopeData$pcgroup, fromLast = TRUE)
mzrt_nonSingleton<-isotopeData[which(nonSingleton==TRUE),]

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

################################
# Metabolite pathway analysis  #
################################

metaboliteDataset<-read.table("day4_differential_metabolitesId.txt",header=TRUE,sep='\t')
metaboliteDataset<-as.data.frame(metaboliteDataset)
metaboliteDataset<-formatMetaboliteData(metaboliteDataset)

mappedMetaboliteData<-mappingMetabolites(batch_corrected_sig_mat_d4,metaboliteDataset)

### Averaging values of multiple hits

#Averaging by compound ids
penul_col<-ncol(mappedMetaboliteData)-1 #as the last column contains p values

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
mappedMetabolites_21strains<-mappedMetabolites[which(mappedMetabolites$samplegroup!="D12_184"),]
mappedMetabolites_21strains<-mappedMetabolites_21strains[which(mappedMetabolites_21strains$samplegroup!="D4_184"),]
dt <- data.table(mappedMetabolites_21strains) # change between mappedMetabolites_21strains and mappedMetabolites to obtain d12_mean_21strains and d12_mean1, respectively.
mappedMetabolites_mean<-dt[,lapply(.SD, mean),by=samplegroup]
mappedMetabolites_mean<-as.data.frame(mappedMetabolites_mean) # 21 or 22 strains
rownames(mappedMetabolites_mean)<-mappedMetabolites_mean[,1]
mappedMetabolites_mean<-mappedMetabolites_mean[,2:ncol(mappedMetabolites_mean)]
# mappedMetabolites_mean_d12_21strains<-mappedMetabolites_mean
#mappedMetabolites_mean_d12<-mappedMetabolites_mean
#write.table(round(mappedMetabolites_mean_d12,4),"metabolitesIdentified_day12.txt",quote=FALSE,sep='\t',col.names=NA)
#mappedMetabolites_mean_d4<-mappedMetabolites_mean
#write.table(round(mappedMetabolites_mean_d4,4),"metabolitesIdentified_day4.txt",quote=FALSE,sep='\t',col.names=NA)

##### Analysis of meta data

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
# biochemData_d12$totalCarbCon   1     71.40  71.395 0.97644 0.04398  0.497  
# Residuals                     16   1169.88  73.118         0.72060         
# Total                         21   1623.50                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

### Mantel test
veg.dist.d12 <- vegdist(as.matrix(mappedMetabolites_mean_d12_21strains),method="euclidean") #d12_mean1[,2:10688]
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
       +biochemData_d4$totalProtein+biochemData_d4$totalCarbCon, method = "euclidean", perm=999)
# Call:
#   adonis(formula = mappedMetabolites_mean_d4 ~ biochemData_d4$biomass +      biochemData_d4$lipidProduc + biochemData_d4$totalLipidCon +      biochemData_d4$totalProtein + biochemData_d4$totalCarbCon,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# biochemData_d4$biomass        1     90.55  90.554 1.08956 0.05011  0.358  
# biochemData_d4$lipidProduc    1     82.77  82.769 0.99589 0.04580  0.424  
# biochemData_d4$totalLipidCon  1    123.05 123.052 1.48058 0.06810  0.070 .
# biochemData_d4$totalProtein   1     81.48  81.480 0.98039 0.04509  0.482  
# biochemData_d4$totalCarbCon   1     99.39  99.388 1.19585 0.05500  0.241  
# Residuals                    16   1329.76  83.110         0.73589         
# Total                        21   1807.00                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

### Mantel test
veg.dist.d4 <- vegdist(as.matrix(mappedMetabolites_mean_d4),method="euclidean") #d12_mean1[,2:10688]
veg.dist.d4.biochem <- vegdist(as.matrix(biochemData_d4),method="euclidean")
mantel(veg.dist.d4, veg.dist.d4.biochem, method="spear")
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d4, ydis = veg.dist.d4.biochem, method = "spear") 
# 
# Mantel statistic r: 0.2208 
#       Significance: 0.032 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.150 0.195 0.234 0.270 
# 
# Based on 999 permutations

veg.dist.d4 <- vegdist(as.matrix(mappedMetabolites_mean_d4),method="euclidean") #d12_mean1[,2:10688]
mantel(veg.dist.d4, distanceMatrix_18s, method="spear")
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d4, ydis = distanceMatrix_18s, method = "spear") 
# 
# Mantel statistic r: -0.07577 
# Significance: 0.697 
# 
# Upper quantiles of permutations (null model):
# 90%   95% 97.5%   99% 
# 0.171 0.217 0.246 0.286 
# 
# Based on 999 permutations

##### Differential metabolites
#####

### Permanova
adonis(mappedMetabolites_mean_d4~ biochemData_d4$biomass+biochemData_d4$lipidProduc+biochemData_d4$totalLipidCon
          +biochemData_d4$totalProtein+biochemData_d4$totalCarbCon, method = "euclidean", perm=999)
# Call:
#   adonis(formula = mappedMetabolites_mean_d4 ~ biochemData_d4$biomass +      biochemData_d4$lipidProduc + biochemData_d4$totalLipidCon +      biochemData_d4$totalProtein + biochemData_d4$totalCarbCon,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# biochemData_d4$biomass        1     48.08  48.085 1.33128 0.06412  0.183
# biochemData_d4$lipidProduc    1     34.04  34.043 0.94253 0.04540  0.518
# biochemData_d4$totalLipidCon  1     23.31  23.312 0.64543 0.03109  0.856
# biochemData_d4$totalProtein   1     23.06  23.059 0.63841 0.03075  0.876
# biochemData_d4$totalCarbCon   1     43.48  43.476 1.20370 0.05798  0.237
# Residuals                    16    577.90  36.119         0.77066       
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
# Mantel statistic r: 0.03536 
#       Significance: 0.361 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.142 0.196 0.231 0.261 
# 
# Based on 999 permutations

veg.dist.d4 <- vegdist(as.matrix(mappedMetabolites_mean_d4),method="euclidean") #d12_mean1[,2:10688]
mantel(veg.dist.d4, distanceMatrix_18s, method="spear")
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d4, ydis = distanceMatrix_18s, method = "spear") 
# 
# Mantel statistic r: 0.1805 
#       Significance: 0.103 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.182 0.230 0.274 0.318 
# 
# Based on 999 permutations

