########################################################
# Analysis of correlations-Malaysian algae data #
########################################################
# Author(s): Shiv
# Version: 21042014
# Input: ".txt" file 
# Modified By :Shivshankar Umashankar 
# Functions written here are used for analyzing correlations between biochemical and metabolic data
# To try: CCA, distLIM
#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages

library(ade4)
library(data.table)
library(vegan)
library(Biostrings)
library(vegan)

##### Reading in the data

# Metabolomics data
setwd('../../AlgaeData/results.from.scelse.cluster.211213/')

load("svd_day4_x138_nonzero.rda")
batch_corrected_mat_d4<-svd_day4_nonzero[[7]]

load("svd_day12_x138_nonzero.rda")
batch_corrected_mat_d12<-svd_day12_nonzero[[4]]


##### 18s rRNA data

Chlorella18s<-readDNAStringSet('Chlorella18s.txt',format="fasta")

distanceMatrix_18s<-stringDist(Chlorella18s, method = "levenshtein", ignoreCase = TRUE, diag = FALSE,
                               upper = FALSE, type = "global")
fit <- hclust(distanceMatrix_18s, method="average") 
pdf('LevenshteinDistance_22strains.pdf')
plot (fit)
dev.off()

##### Biochemical data

setwd('../../data/Biochemical_042014/')

biochemData0<-read.table("biochemicalData.txt",header=TRUE,sep='\t')
biochemData<-as.data.frame(biochemData0[,3:8])
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
biochemData_d4<-biochemData_d4[,2:7]
biochemData_d4<-as.data.frame(scale(biochemData_d4,center=TRUE, scale=TRUE))

biochemData_d12$samplegroup<-as.factor(gsub ("a|b|c","",rownames(biochemData_d12)))
biochemData_d12 <- data.table(biochemData_d12)
biochemData_d12<-biochemData_d12[,lapply(.SD, mean),by=samplegroup]
biochemData_d12<-as.data.frame(biochemData_d12)
rownames(biochemData_d12)<-biochemData_d12[,1]
biochemData_d12<-biochemData_d12[,2:7]
biochemData_d12<-as.data.frame(scale(biochemData_d12,center=TRUE, scale=TRUE))

##### Obtaining sample groups 

#d12_test<-as.data.frame(t(batch_corrected_mat_d12-min(batch_corrected_mat_d12)))
#d4_test<-as.data.frame(t(batch_corrected_mat_d4-min(batch_corrected_mat_d4)))

SampleGroups<-sapply(rownames(d4_test), function(x) paste(strsplit(x,"_")[[1]][1:2],collapse = "_")) #change d12_test or d4_test accordingly
SampleGroups<-as.vector(SampleGroups)
 
SampleGroups<-gsub('b1','',SampleGroups)
SampleGroups<-gsub('b2','',SampleGroups)
SampleGroups<-gsub('b3','',SampleGroups)

#Metabolomics data- converting into positive values 

d4_test<-as.data.frame(t(batch_corrected_mat_d4-min(batch_corrected_mat_d4)))
d4_test$samplegroup<-as.factor(SampleGroups)
# correct naming of strains
d4_test$samplegroup<-as.factor(gsub('D4_14','D4_014',d4_test$samplegroup))
d4_test<-d4_test[order(d4_test$samplegroup),]
d4_test_21strains<-d4_test[which(d4_test$samplegroup!="D4_184"),]
dt <- data.table(d4_test_21strains) # change between d4_test_21strains and d4_test to obtain d4_mean_21strains and d4_mean1, respectively.
d4_mean<-dt[,lapply(.SD, mean),by=samplegroup]
#d4_mean1<-as.data.frame(d4_mean) #22 strains
d4_mean_21strains<-as.data.frame(d4_mean)

d12_test<-as.data.frame(t(batch_corrected_mat_d12-min(batch_corrected_mat_d12)))
d12_test$samplegroup<-as.factor(SampleGroups)
# correct naming of strains
d12_test$samplegroup<-as.factor(gsub('D12_14','D12_014',d12_test$samplegroup))
d12_test$samplegroup<-as.factor(gsub('D12_84','D12_084',d12_test$samplegroup))
d12_test$samplegroup<-as.factor(gsub('D12_87','D12_087',d12_test$samplegroup))
d12_test$samplegroup<-as.factor(gsub('D12_94','D12_094',d12_test$samplegroup))
d12_test<-d12_test[order(d12_test$samplegroup),]
d12_test_21strains<-d12_test[which(d12_test$samplegroup!="D12_184"),]
dt <- data.table(d12_test_21strains) # change between d12_test_21strains and d12_test to obtain d12_mean_21strains and d12_mean1, respectively.
d12_mean<-dt[,lapply(.SD, mean),by=samplegroup]
#d12_mean1<-as.data.frame(d12_mean) # 22 strains
d12_mean_21strains<-as.data.frame(d12_mean)

### mantel test 

#using vegan
veg.dist.d4 <- vegdist(as.matrix(d4_mean_21strains[,2:13444]),method="euclidean") # Bray-Curtis
veg.dist.d12 <- vegdist(as.matrix(d12_mean_21strains[,2:10688]),method="euclidean")
mantel(veg.dist.d4, veg.dist.d12)
mantel(veg.dist.d4, veg.dist.d12, method="spear")

#using ade4

veg.dist.d4 <- dist(as.matrix(d4_mean1[,2:13443])) # Bray-Curtis
veg.dist.d12 <- dist(as.matrix(d12_mean1[,2:10688]))
mantel.rtest(veg.dist.d4, veg.dist.d12,nrepet=999)

###############################
# Biochemical data
veg.dist.d12 <- vegdist(as.matrix(mappedMetabolites_mean_d12),method="euclidean") #d12_mean1[,2:10688]
veg.dist.d12.biochem <- vegdist(as.matrix(biochemData_d12),method="euclidean")
mantel(veg.dist.d12, veg.dist.d12.biochem)
mantel(veg.dist.d12, veg.dist.d12.biochem, method="spear")
mantel.rtest(veg.dist.d12, veg.dist.d12.biochem,nrepet=999)

# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d12, ydis = veg.dist.d12.biochem, method = "spear") 
# 
# Mantel statistic r: 0.258 
#       Significance: 0.016 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.152 0.197 0.230 0.271 
# 
# Based on 999 permutations


veg.dist.d4 <- vegdist(as.matrix(mappedMetabolites_mean_d4),method="euclidean")#d4_mean1[,2:10688]
veg.dist.d4.biochem <- vegdist(as.matrix(biochemData_d4),method="euclidean")
mantel(veg.dist.d4, veg.dist.d4.biochem)
mantel(veg.dist.d4, veg.dist.d4.biochem, method="spear")
mantel.rtest(veg.dist.d4, veg.dist.d4.biochem,nrepet=999)

# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d4, ydis = veg.dist.d4.biochem, method = "spear") 
# 
# Mantel statistic r: 0.01616 
#       Significance: 0.45 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.150 0.204 0.240 0.282 
# 
# Based on 999 permutations

###############################
# 18s rRNA data

mantel(veg.dist.d12, distanceMatrix_18s, method="spear")
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d12, ydis = distanceMatrix_18s, method = "spear") 
# 
# Mantel statistic r: -0.1005 
#       Significance: 0.767 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.176 0.220 0.258 0.304 
# 
# Based on 999 permutations

mantel.rtest(veg.dist.d12, distanceMatrix_18s,nrepet=999)
# Monte-Carlo test
# Observation: -0.08679852 
# Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
# Based on 999 replicates
# Simulated p-value: 0.722 

mantel(veg.dist.d4, distanceMatrix_18s, method="spear")
# Mantel statistic based on Spearman's rank correlation rho 
# 
# Call:
# mantel(xdis = veg.dist.d4, ydis = distanceMatrix_18s, method = "spear") 
# 
# Mantel statistic r: 0.1316 
#       Significance: 0.193 
# 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.187 0.235 0.281 0.322 
# 
# Based on 999 permutations

mantel.rtest(veg.dist.d4, distanceMatrix_18s,nrepet=999)
# Monte-Carlo test
# Observation: 0.1491506 
# Call: mantelnoneuclid(m1 = m1, m2 = m2, nrepet = nrepet)
# Based on 999 replicates
# Simulated p-value: 0.158 