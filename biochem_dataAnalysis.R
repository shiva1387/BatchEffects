###########################################
# Data Analysis in R-Malaysian algae data #
###########################################
# Author(s): Shiv
# Version: 07052014
# Input: ".txt" file 
# Modified By :Shivshankar Umashankar 
# Functions written here are used for analyzing algal biochemical data

#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 

library(devtools)
#install_github("ggbiplot", "vqv")
library(ggbiplot)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
library(gridExtra)
library(gsubfn)
library(vegan)

### Reading Data

biochemData0<-read.table("biochemicalData.txt",header=TRUE,sep='\t')
biochemData<-as.data.frame(biochemData0[,3:8])
rownames(biochemData)<-biochemData0[,2]
biochemData<-scale(biochemData,center=TRUE,scale=TRUE) # shiv edited on May 7 2014. Added centering and scaling at the start
biochemData_d4<-biochemData[grep('D4', rownames(biochemData)),] 
biochemData_d12<-biochemData[grep('D12', rownames(biochemData)),] 

# Analysis
d<-biochemData_d4[,apply(biochemData_d4, 2, var, na.rm=TRUE) != 0] #removing mz features which have constant variance
fit_d <- prcomp(d) #uses SVD # #shiv edited on May 7 2014. Removed scaling as data has been scaled at the start
scores_pca_d<-as.data.frame(fit_d$x)

# plot(fit_d$x[,1],fit_d$x[,2])
# text(fit_d$x[,1],fit_d$x[,2],rownames(d),pos=3)
  
pcaPlotOutput<-pcaPlot(fit_d,rownames(d),gsub ("a|b|c","",rownames(biochemData_d4)))#biochemData0[,1])
print(pcaPlotOutput)
ggsave("bioChem_pcaPlotOutput_D4.pdf",pcaPlotOutput)
dev.off()

#correlation analysis and plot

correlationPlotOutput<-correlationPlot(biochemData)
print(correlationPlotOutput)
ggsave("bioChem_correlationPlotOutput.pdf",correlationPlotOutput)
dev.off()



# Enlarging dataset for metabolomics data comparison

biochemDataForMetab<-biochemData[rep(1:nrow(biochemData),each=2),] # duplicating each row to avvount for technical replicates
rownames(biochemDataForMetab)<-gsubfn(".", list("a"="b1_r001", "a.1"="b1_r002", "b"="b2_r001","b.1"="b2_r002","c"="b3_r001","c.1"="b3_r002"),rownames(biochemDataForMetab))
rownames(biochemDataForMetab)<-gsub("r001.1","r002",rownames(biochemDataForMetab))

# Day12 
batch_corrected_mat_d12<-svd_day12_nonzero[[4]]
colnames(batch_corrected_mat_d12)<-gsub('D12_14','D12_014',colnames(batch_corrected_mat_d12))
colnames(batch_corrected_mat_d12)<-gsub('D12_84','D12_084',colnames(batch_corrected_mat_d12))
colnames(batch_corrected_mat_d12)<-gsub('D12_87','D12_087',colnames(batch_corrected_mat_d12))
colnames(batch_corrected_mat_d12)<-gsub('D12_94','D12_094',colnames(batch_corrected_mat_d12))

#batch_corrected_mat_d12_sig<-day12_nonzero_sigfeat_s[[2]]
#batch_corrected_mat_d12_sig<-batch_corrected_mat_d12_sig[[4]] 
#batch_corrected_mat_d12_sig<-batch_corrected_mat_d12_sig[,1:118]#removing last column which has wrong pvalues
biochemDataForMetab_d12<-biochemDataForMetab[colnames(batch_corrected_mat_d12),]
biochemDataForMetab_d12<-biochemDataForMetab_d12[complete.cases(biochemDataForMetab_d12),]
biochemDataForMetab_d12$Strain<-gsub("b1_r001|b1_r002|b2_r001|b2_r002|b3_r001|b3_r002","",rownames(biochemDataForMetab_d12))

png("biochem_day12_pairs.png")
pairs(biochemDataForMetab_d12[,1:6],main="DAY 12")
dev.off()

a<-adonis(t(batch_corrected_mat_d12)~ biochemDataForMetab_d12$biomass+biochemDataForMetab_d12$biomassProduc+biochemDataForMetab_d12$lipidProduc+
            biochemDataForMetab_d12$totalLipidCon+biochemDataForMetab_d12$totalProtein+biochemDataForMetab_d12$totalCarbCon+biochemDataForMetab_d12$Strain, 
          method = "euclidean", perm=999)

# Call:
#   adonis(formula = t(batch_corrected_mat_d12) ~ biochemDataForMetab_d12$biomass +      biochemDataForMetab_d12$biomassProduc + biochemDataForMetab_d12$lipidProduc +      biochemDataForMetab_d12$totalLipidCon + biochemDataForMetab_d12$totalProtein +      biochemDataForMetab_d12$totalCarbCon + biochemDataForMetab_d12$Strain,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# biochemDataForMetab_d12$biomass         1      9906  9905.9  8.1659 0.04124  0.001 ***
#   biochemDataForMetab_d12$biomassProduc   1      5444  5443.8  4.4876 0.02266  0.001 ***
#   biochemDataForMetab_d12$lipidProduc     1      4965  4965.4  4.0932 0.02067  0.001 ***
#   biochemDataForMetab_d12$totalLipidCon   1      7444  7444.0  6.1364 0.03099  0.001 ***
#   biochemDataForMetab_d12$totalProtein    1      4914  4914.5  4.0512 0.02046  0.001 ***
#   biochemDataForMetab_d12$totalCarbCon    1      4150  4150.3  3.4213 0.01728  0.001 ***
#   biochemDataForMetab_d12$Strain         21     94223  4486.8  3.6987 0.39223  0.001 ***
#   Residuals                              90    109178  1213.1         0.45448           
# Total                                 117    240225                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

a<-adonis(t(batch_corrected_mat_d12)~ biochemDataForMetab_d12$biomass+biochemDataForMetab_d12$lipidProduc+
            biochemDataForMetab_d12$totalLipidCon+biochemDataForMetab_d12$totalProtein+biochemDataForMetab_d12$totalCarbCon+biochemDataForMetab_d12$Strain, 
          method = "euclidean", perm=999)

# Call:
#   adonis(formula = t(batch_corrected_mat_d12) ~ biochemDataForMetab_d12$biomass +      biochemDataForMetab_d12$lipidProduc + biochemDataForMetab_d12$totalLipidCon +      biochemDataForMetab_d12$totalProtein + biochemDataForMetab_d12$totalCarbCon +      biochemDataForMetab_d12$Strain, permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# biochemDataForMetab_d12$biomass         1      9906  9905.9  8.1567 0.04124  0.001 ***
#   biochemDataForMetab_d12$lipidProduc     1      5552  5551.6  4.5713 0.02311  0.001 ***
#   biochemDataForMetab_d12$totalLipidCon   1      7185  7184.9  5.9161 0.02991  0.001 ***
#   biochemDataForMetab_d12$totalProtein    1      5514  5514.3  4.5405 0.02295  0.001 ***
#   biochemDataForMetab_d12$totalCarbCon    1      4090  4090.2  3.3680 0.01703  0.001 ***
#   biochemDataForMetab_d12$Strain         21     97463  4641.1  3.8216 0.40572  0.001 ***
#   Residuals                              91    110515  1214.5         0.46005           
# Total                                 117    240225                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

permutedData<-batch_corrected_mat_d12[,sample(ncol(batch_corrected_mat_d12))]
a<-adonis(t(permutedData)~ biochemDataForMetab_d12$biomass+biochemDataForMetab_d12$lipidProduc+
            biochemDataForMetab_d12$totalLipidCon+biochemDataForMetab_d12$totalProtein+biochemDataForMetab_d12$totalCarbCon+biochemDataForMetab_d12$Strain, 
          method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(permutedData) ~ biochemDataForMetab_d12$biomass +      biochemDataForMetab_d12$lipidProduc + biochemDataForMetab_d12$totalLipidCon +      biochemDataForMetab_d12$totalProtein + biochemDataForMetab_d12$totalCarbCon +      biochemDataForMetab_d12$Strain, permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# biochemDataForMetab_d12$biomass         1      2433  2432.6 1.18240 0.01013  0.195
# biochemDataForMetab_d12$lipidProduc     1      1711  1710.9 0.83161 0.00712  0.792
# biochemDataForMetab_d12$totalLipidCon   1      2507  2506.7 1.21845 0.01043  0.163
# biochemDataForMetab_d12$totalProtein    1      1918  1918.4 0.93248 0.00799  0.565
# biochemDataForMetab_d12$totalCarbCon    1      2369  2369.3 1.15167 0.00986  0.221
# biochemDataForMetab_d12$Strain         21     42072  2003.4 0.97381 0.17514  0.681
# Residuals                              91    187215  2057.3         0.77933       
# Total                                 117    240225                 1.00000   

# non-batch effect removed data

#ms_data_day12_nonzero

non_batch_corrrected_mat_d12<-log(ms_data_day12_nonzero)
colnames(non_batch_corrrected_mat_d12)<-gsub("_b","b",colnames(non_batch_corrrected_mat_d12))
biochemDataForMetab_d12_nb<-biochemDataForMetab[colnames(non_batch_corrrected_mat_d12),]
biochemDataForMetab_d12_nb<-biochemDataForMetab_d12_nb[complete.cases(biochemDataForMetab_d12_nb),]
biochemDataForMetab_d12_nb$Strain<-gsub("b1_r001|b1_r002|b2_r001|b2_r002|b3_r001|b3_r002","",rownames(biochemDataForMetab_d12_nb))

pairs(biochemDataForMetab_d12_nb[,1:6],main="DAY 12 NB")

#sum(is.na(biochemDataForMetab_d12_nb))

a<-adonis(t(non_batch_corrrected_mat_d12)~ biochemDataForMetab_d12_nb$biomass+biochemDataForMetab_d12_nb$biomassProduc+biochemDataForMetab_d12_nb$lipidProduc+
            biochemDataForMetab_d12_nb$totalLipidCon+biochemDataForMetab_d12_nb$totalProtein+biochemDataForMetab_d12_nb$totalCarbCon+biochemDataForMetab_d12_nb$Strain, 
          method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(non_batch_corrrected_mat_d12) ~ biochemDataForMetab_d12_nb$biomass +      biochemDataForMetab_d12_nb$biomassProduc + biochemDataForMetab_d12_nb$lipidProduc +      biochemDataForMetab_d12_nb$totalLipidCon + biochemDataForMetab_d12_nb$totalProtein +      biochemDataForMetab_d12_nb$totalCarbCon + biochemDataForMetab_d12_nb$Strain,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# biochemDataForMetab_d12_nb$biomass         1     34916   34916  9.9810 0.03384  0.001 ***
#   biochemDataForMetab_d12_nb$biomassProduc   1     14946   14946  4.2724 0.01449  0.001 ***
#   biochemDataForMetab_d12_nb$lipidProduc     1     31238   31238  8.9298 0.03028  0.001 ***
#   biochemDataForMetab_d12_nb$totalLipidCon   1     24841   24841  7.1011 0.02408  0.001 ***
#   biochemDataForMetab_d12_nb$totalProtein    1     33220   33220  9.4962 0.03220  0.001 ***
#   biochemDataForMetab_d12_nb$totalCarbCon    1     19560   19560  5.5916 0.01896  0.001 ***
#   biochemDataForMetab_d12_nb$Strain         21    558131   26578  7.5975 0.54099  0.001 ***
#   Residuals                                 90    314839    3498         0.30517           
# Total                                    117   1031692                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# Day4 
batch_corrected_mat_d4<-svd_day4_nonzero[[7]]#day4_nonzero_sigfeat_s[[2]]
#batch_corrected_mat_d4<-batch_corrected_mat_d4[[7]] #removing last column which has wrong pvalues
#batch_corrected_mat_d4<-batch_corrected_mat_d4[,1:125] #removing last column which has pvalues
colnames(batch_corrected_mat_d4)<-gsub("_b","b",colnames(batch_corrected_mat_d4))
colnames(batch_corrected_mat_d4)<-gsub('D4_14','D4_014',colnames(batch_corrected_mat_d4))

biochemDataForMetab_d4<-biochemDataForMetab[colnames(batch_corrected_mat_d4),]
biochemDataForMetab_d4<-biochemDataForMetab_d4[complete.cases(biochemDataForMetab_d4),]
biochemDataForMetab_d4$Strain<-gsub("b1_r001|b1_r002|b2_r001|b2_r002|b3_r001|b3_r002","",rownames(biochemDataForMetab_d4))

png("biochem_day4_pairs.png")
pairs(biochemDataForMetab_d4[,1:6],main="DAY 4")
dev.off()


#sum(is.na(biochemDataForMetab_d4))

a<-adonis(t(batch_corrected_mat_d4)~ biochemDataForMetab_d4$biomass+biochemDataForMetab_d4$biomassProduc+biochemDataForMetab_d4$lipidProduc+
            biochemDataForMetab_d4$totalLipidCon+biochemDataForMetab_d4$totalProtein+biochemDataForMetab_d4$totalCarbCon+biochemDataForMetab_d4$Strain, 
          method = "euclidean", perm=999)

# Call:
#   adonis(formula = t(batch_corrected_mat_d4) ~ biochemDataForMetab_d4$biomass +      biochemDataForMetab_d4$biomassProduc + biochemDataForMetab_d4$lipidProduc +      biochemDataForMetab_d4$totalLipidCon + biochemDataForMetab_d4$totalProtein +      biochemDataForMetab_d4$totalCarbCon + biochemDataForMetab_d4$Strain,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# biochemDataForMetab_d4$biomass         1      4771  4771.0  3.6367 0.02124  0.001 ***
#   biochemDataForMetab_d4$biomassProduc   1      3792  3792.1  2.8905 0.01688  0.001 ***
#   biochemDataForMetab_d4$lipidProduc     1      3709  3709.4  2.8275 0.01652  0.001 ***
#   biochemDataForMetab_d4$totalLipidCon   1      4962  4962.2  3.7825 0.02209  0.001 ***
#   biochemDataForMetab_d4$totalProtein    1      2520  2520.2  1.9210 0.01122  0.001 ***
#   biochemDataForMetab_d4$totalCarbCon    1      2827  2827.2  2.1550 0.01259  0.001 ***
#   biochemDataForMetab_d4$Strain         21     74771  3560.5  2.7140 0.33290  0.001 ***
#   Residuals                             97    127254  1311.9         0.56656           
# Total                                124    224608                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

a<-adonis(t(batch_corrected_mat_d4)~ biochemDataForMetab_d4$biomass+biochemDataForMetab_d4$lipidProduc+
            biochemDataForMetab_d4$totalLipidCon+biochemDataForMetab_d4$totalProtein+biochemDataForMetab_d4$totalCarbCon+biochemDataForMetab_d4$Strain, 
          method = "euclidean", perm=999)

# Call:
#   adonis(formula = t(batch_corrected_mat_d4) ~ biochemDataForMetab_d4$biomass +      biochemDataForMetab_d4$lipidProduc + biochemDataForMetab_d4$totalLipidCon +      biochemDataForMetab_d4$totalProtein + biochemDataForMetab_d4$totalCarbCon +      biochemDataForMetab_d4$Strain, permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# biochemDataForMetab_d4$biomass         1      4771  4771.0  3.6187 0.02124  0.001 ***
#   biochemDataForMetab_d4$lipidProduc     1      3993  3992.6  3.0283 0.01778  0.001 ***
#   biochemDataForMetab_d4$totalLipidCon   1      3165  3164.7  2.4004 0.01409  0.001 ***
#   biochemDataForMetab_d4$totalProtein    1      2331  2331.3  1.7682 0.01038  0.001 ***
#   biochemDataForMetab_d4$totalCarbCon    1      3492  3491.8  2.6484 0.01555  0.001 ***
#   biochemDataForMetab_d4$Strain         21     77651  3697.6  2.8046 0.34572  0.001 ***
#   Residuals                             98    129206  1318.4         0.57525           
# Total                                124    224608                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

########### CCA

cca.1 <- cca(t(batch_corrected_mat_d4)~ biochemDataForMetab_d4$biomass+biochemDataForMetab_d4$lipidProduc+
               biochemDataForMetab_d4$totalLipidCon+biochemDataForMetab_d4$totalProtein+biochemDataForMetab_d4$totalCarbCon+biochemDataForMetab_d4$Strain)

########## Differential feature identification

### Comparisons

Chlorella<-c("D4_001","D4_006","D4_051","D4_104","D4_14","D4_177","D4_187","D4_207","D4_268","D4_283","D4_322","D4_325","D4_84","D4_87","D4_94")
Parachlorella<-c("D4_245","D4_252","D4_253","D4_254","D4_255","D4_258")
UnIdchlorella<-"D4_184"

AerobicPond<-"D12_051"
SeaBassPond<-c("D12_252","D12_253","D12_2512","D12_255","D12_258")

HighPalmitic<-c("D12_187","D12_84","D12_283","D12_001")
OtherPalmitic<-c("D12_006","D12_14","D12_051","D12_87","D12_94","D12_104","D12_177","D12_184","D12_207",
                 "D12_245","D12_252","D12_253","D12_254","D12_255","D12_258","D12_268","D12_322","D12_325")

ms_data_day12_chlo<-batch_corrected_mat_d12[, grepl(paste(Chlorella,collapse="|"),colnames(batch_corrected_mat_d12))]
ms_data_day12_parchlo<-batch_corrected_mat_d12[, grepl(paste(Parachlorella,collapse="|"),colnames(batch_corrected_mat_d12))]
ms_data_day12_unidchlo<-batch_corrected_mat_d12[, grepl(paste(UnIdchlorella,collapse="|"),colnames(batch_corrected_mat_d12))]

ms_data_day12_aerobic<-batch_corrected_mat_d12[, grepl(paste(AerobicPond,collapse="|"),colnames(batch_corrected_mat_d12))]
ms_data_day12_seabass<-batch_corrected_mat_d12[, grepl(paste(SeaBassPond,collapse="|"),colnames(batch_corrected_mat_d12))]

ms_data_day12_highPalm<-batch_corrected_mat_d12[, grepl(paste(HighPalmitic,collapse="|"),colnames(batch_corrected_mat_d12))]
ms_data_day12_otherPalm<-batch_corrected_mat_d12[, grepl(paste(OtherPalmitic,collapse="|"),colnames(batch_corrected_mat_d12))]


ms_data_day4_chlo<-batch_corrected_mat_d4[, grepl(paste(Chlorella,collapse="|"),colnames(batch_corrected_mat_d4))]
ms_data_day4_parchlo<-batch_corrected_mat_d4[, grepl(paste(Parachlorella,collapse="|"),colnames(batch_corrected_mat_d4))]
ms_data_day4_unidchlo<-batch_corrected_mat_d4[, grepl(paste(UnIdchlorella,collapse="|"),colnames(batch_corrected_mat_d4))]

ms_data_day4_aerobic<-batch_corrected_mat_d4[, grepl(paste(AerobicPond,collapse="|"),colnames(batch_corrected_mat_d4))]
ms_data_day4_seabass<-batch_corrected_mat_d4[, grepl(paste(SeaBassPond,collapse="|"),colnames(batch_corrected_mat_d4))]

### Sample groups

ms_data_day4_chloParchlo<-cbind(ms_data_day4_chlo,ms_data_day4_parchlo)
SampleGroup_chloParchlo<-c(rep("chlorella",ncol(ms_data_day4_chlo)),rep("Parchlorella",ncol(ms_data_day4_parchlo)))

ms_data_day12_aeroSea<-cbind(ms_data_day12_aerobic,ms_data_day12_seabass)
SampleGroup_aeroSea<-c(rep("Aerobic",ncol(ms_data_day12_aerobic)),rep("SeaBass",ncol(ms_data_day12_seabass)))

ms_data_day12_highVsotherPalm<-cbind(ms_data_day12_highPalm,ms_data_day12_otherPalm)
SampleGroup_highVsotherPalm<-c(rep("highPalm",ncol(ms_data_day12_highPalm)),rep("otherPalm",ncol(ms_data_day12_otherPalm)))


SampleGroup<-sapply(colnames(batch_corrected_mat_d12), function(x) paste(strsplit(x,"_")[[1]][1:2],collapse = "_"))
SampleGroup<-as.vector(SampleGroup)
SampleGroup<-gsub('b1','',SampleGroup)
SampleGroup<-gsub('b2','',SampleGroup)
SampleGroup<-gsub('b3','',SampleGroup)
SampleGroup#Sample group without the outlier replicates

######## Permutation based test statistic

classlabel_factor<-as.numeric(as.factor(SampleGroup))-1
data_matrix<-as.matrix(batch_corrected_mat_d12)
dataset_sig_features<-mt.maxT(data_matrix,classlabel_factor,test="f",side="abs",fixed.seed.sampling="y",B=1000,nonpara="n")
# dataset_sig_featureold<-dataset_sig_features
# dataset_sig_features<-as.matrix(dataset_sig_features)
p_values_sig<-sort(dataset_sig_features$rawp) 
id.sig_dataset <- sort(dataset_sig_features[dataset_sig_features$adjp < 0.05,c(1)]) #getting the column which provides index of rows satisying the condition
metab_sig<-cbind(data_matrix[id.sig_dataset,],round(dataset_sig_features$adjp[dataset_sig_features$index %in% id.sig_dataset],5))
# 
# No significant features which were different between chlorella and parachlorella at day12 (f test,FDR) or day 4
# Features are significantly different between aerobic pond Vs SeaBass pond at day 12 -131 features and day4-98 features
# Features are significantly different between high palmitic strains Vs others at day12-109 features and day4-  features 


# #########################################
# # now let's use multtest to get adjusted p-values
# procs <- c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY") 
# adj.all <- mt.rawp2adjp(p_values, procs) 
# adj.sig <- mt.rawp2adjp(p_values_sig, procs)
# 
# #########################################
# # compare number of bonferroni corrected p-values significant at p < 0.05 between the first 500 and full 1000 genes
# length(which(adj.sig$adjp[,8] <= 0.05))
# length(which(adj.all$adjp[,1] <= 0.05))
# #########################################
# # compare results by plotting raw versus adjusted p-values
# plot(adj.all$adjp[,1], adj.all$adjp[,2], pch = 19)  #Bonferroni
# points(adj.all$adjp[,1], adj.all$adjp[,3], col = "red", pch = 19) #Holm
# points(adj.all$adjp[,1], adj.all$adjp[,7], col = "green", pch = 19) #BH
# 
# ##
# mz_rt <- strsplit(rownames(metab_sig), "\\@")
# mz<-sapply(mz_rt , function (x) if(length(x) == 2) x[1] else as.character(NA))
# rt<-sapply(mz_rt , function (x) if(length(x) == 2) x[2] else as.character(NA))
# mz_rt1 <- strsplit(rownames(metab_sig1), "\\@")
# mz1<-sapply(mz_rt1 , function (x) if(length(x) == 2) x[1] else as.character(NA))
# rt1<-sapply(mz_rt1 , function (x) if(length(x) == 2) x[2] else as.character(NA))
# 
# plot(mz,rt,col="#00000033")
# points(mz1,rt1,col="green")


###############################################################################
################ Significant metabolites ## from metab_sig ####################
###############################################################################

##### CAMERA Analysis

isotopeData<-read.table('algae_4setblanks_021213_mzrt_camera.txt',header=TRUE,sep='\t')
#isotopeData$mzrt<-paste0(round(as.numeric(isotopeData$mz),5),'@',isotopeData$rt)
isotopeData$mz<-round(as.numeric(isotopeData$mz),5)
isotopeData$rt<-round(as.numeric(isotopeData$mz),0)
isotopeData$mzrt<-paste0(isotopeData$mz,'@',isotopeData$rt)
nonSingleton<-duplicated(isotopeData$pcgroup) | duplicated(isotopeData$pcgroup, fromLast = TRUE)
mzrt_nonSingleton<-isotopeData[which(nonSingleton==TRUE),]
#23775

sigfeat_day12_mzAndrt<-strsplit(rownames(metab_sig), "\\@")
sigfeat_day12_mz<-sapply(sigfeat_day12_mzAndrt , function (x) if(length(x) == 2) x[1] else as.character(NA))
sigfeat_day12_rt<-sapply(sigfeat_day12_mzAndrt , function (x) if(length(x) == 2) x[2] else as.character(NA))
sigfeat_day12_mz<-round(as.numeric(sigfeat_day12_mz),5)
sigfeat_day12_rt<-round(as.numeric(sigfeat_day12_rt),0)
sigfeat_day12_mzrt<-paste0(sigfeat_day12_mz,'@',sigfeat_day12_rt)
sigfeat_day12<-data.frame(sigfeat_day12_mzrt,sigfeat_day12_mz,sigfeat_day12_rt)

#functions to get matching mz and rt
cmpMZ <- function(mz, mzlist, mz_cutoff){abs(mz-mzlist) <= mz_cutoff}
#if there is more than one match, change function to pick the closest!
cmpRT <- function(ret, retlist,ret_cutoff){abs(ret-retlist) <= ret_cutoff}

###day12
mz_tolerance=10 #ppm
rt_tolerance=5 #secs
for(i in 1:length(sigfeat_day12))
{
  
  ppm_tolerance<-(trunc(sigfeat_day12_mz[i])*mz_tolerance)/10^6 
  match_mz<-cmpMZ(sigfeat_day12_mz[i],mzrt_nonSingleton$mz,ppm_tolerance)
  match_mz_ind<-which(match_mz) #obtaining indices only for the significant features
  match_ret<-cmpRT(sigfeat_day12_rt[i],mzrt_nonSingleton$rt[match_mz_ind],rt_tolerance)
  mz_rt_match<-match_mz_ind[which(match_ret)]
  matched_value<-mzrt_nonSingleton[mz_rt_match,"mzrt"];
  matched_value=ifelse(exists("matched_value"),matched_value)
  #print(c(i,sigfeat_day12_bc_mz_s[i],sigfeat_day12_bc_rt_s[i],mz_rt_match,matched_value))
  sigfeat_day12[i,4]<-matched_value
}
# write.table(sigfeat_day12,"HighpalmVsotherDay12.txt",quote=FALSE,sep='\t')
# sigfeat_day12_nonsingletons<-sigfeat_day12[!is.na(sigfeat_day12$V4),]
# sigfeat_day12_nonsingletonList_mz<-strsplit(sigfeat_day12_bc_nonsingletons$V4, "\\@")
#writing mz rt to a file


#########################################################
####################### Functions #######################
#########################################################
# pca plot

pcaPlot<-function(pca_summary,dataLabel,dataGroups){
  g <- ggbiplot(pca_summary, obs.scale = 1, var.scale = 1, ellipse = FALSE, groups=dataGroups,labels=dataLabel,
                circle = TRUE)
  g <- g + scale_color_discrete(name = '')
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
  return(g)
}

#function to compute correlation plot

correlationPlot<-function(env_parameters){
combos <- combn(ncol(env_parameters),2) # combinations without repetitions

#combos <- expand.grid(rep(list(1:ncol(env_parameters)), 2 )) # combinations with repetitions
#combos <- as.matrix(combos)
#combos <- t(combos) # transpose matrix

mat1 <- adply(combos, 2, function(x) {
  #test <- chisq.test(env_parameters[, x[1]], env_parameters[, x[2]])
  test<-cor.test(env_parameters[, x[1]], env_parameters[, x[2]], use="all.obs", method="kendall", exact=FALSE) 
  
  out <- data.frame("Row" = colnames(env_parameters)[x[1]]
                    , "Column" = colnames(env_parameters[x[2]])
                    , "p.value" = round(test$p.value, 3)
                    ,  "statistic"= test$statistic
                    ,  "tau.value" = round(test$estimate,3)
  )
  return(out)
})

tau.valu.mat<-mat1
# mat 2 is the same table but formated to show pvalue of corrrelations
mat2 <- adply(combos, 2, function(x) {
  #test <- chisq.test(env_parameters[, x[1]], env_parameters[, x[2]])
  test<-cor.test(env_parameters[, x[1]], env_parameters[, x[2]], use="all.obs", method="kendall", exact=FALSE) 
  
  out <- data.frame("Row" = colnames(env_parameters)[x[1]]
                    , "Column" = colnames(env_parameters[x[2]])
                    , "tau.value" = round(test$estimate,3)
                    , "statistic"= test$statistic
                    , "p.value" = round(test$p.value, 3)
  )
  return(out)
})

p.valu.mat<-mat2

### ploting using ggplot2

mat1$stars <- cut(mat1$p.value, breaks=c(-Inf, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
mat1$tauval.scal<-cut(mat1$tau.value,breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),include.lowest=TRUE, label=c("(-0.75,-1)","(-0.5,-0.75)","(-0.25,-0.5)","(0,-0.25)","(0,0.25)","(0.25,0.5)","(0.5,0.75)","(0.75,1)")) # the label for the legend

meas <- as.character(unique(mat1$Column))
#pdf('env_parameters_sediment.pdf',height=16,width=16)

p<-ggplot(mat1, aes(Row, Column)) +coord_flip()+
  geom_tile(aes(fill=mat1$tauval.scal),colour="white") +  
  scale_fill_brewer(palette = "RdYlGn",name="Correlation") + geom_text(aes(label=mat1$stars), color="black", size=5.5) + 
  labs(y=NULL, x=NULL, fill="p.value") +  scale_y_discrete(limits=meas[length(meas):1]) + #flip the x axis 
  theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1,size=18),axis.text.y=element_text(size=18),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black')) #draws x and y axis line

p1<-p+theme(legend.text = element_text(size = 15),legend.title = element_text(size = 20)) #change the font size for legend
return(p1)
}