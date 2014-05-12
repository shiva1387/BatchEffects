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
library(Biostrings)
library(psych)

### Reading Data

biochemData0<-read.table("biochemicalData.txt",header=TRUE,sep='\t')
biochemData<-as.data.frame(biochemData0[,3:8])
rownames(biochemData)<-biochemData0[,2]
#biochemData<-as.data.frame(scale(biochemData,center=TRUE,scale=TRUE)) # shiv edited on May 7 2014. Added centering and scaling at the start
biochemData_d4<-biochemData[grep('D4', rownames(biochemData)),] 
biochemData_d4<-as.data.frame(scale(biochemData_d4,center=TRUE, scale=TRUE))
biochemData_d12<-biochemData[grep('D12', rownames(biochemData)),] 
biochemData_d12<-as.data.frame(scale(biochemData_d12,center=TRUE, scale=TRUE))

biochemData1<-as.data.frame(scale(biochemData,center=TRUE,scale=TRUE)) #scaling and centering the complete matrix
biochemData1<-biochemData-min(biochemData)

# Analysis
d<-biochemData1[,apply(biochemData1, 2, var, na.rm=TRUE) != 0] #removing mz features which have constant variance
fit_d <- prcomp(d) #uses SVD # #shiv edited on May 7 2014. Removed scaling as data has been scaled at the start
scores_pca_d<-as.data.frame(fit_d$x)

# plot(fit_d$x[,1],fit_d$x[,2])
# text(fit_d$x[,1],fit_d$x[,2],rownames(d),pos=3)
  
pcaPlotOutput<-pcaPlot(fit_d,rownames(d),biochemData0[,1])#gsub ("a|b|c","",rownames(biochemData_d12)))#biochemData1[,1])
#print(pcaPlotOutput)
ggsave("bioChem_pcaPlotOutput.pdf",pcaPlotOutput)
dev.off()

#correlation analysis and plot

correlationPlotOutput<-correlationPlot(biochemData1)
print(correlationPlotOutput)
ggsave("bioChem_correlationPlotOutput.pdf",p1)
dev.off()

#Pairs plot

png("biochem_day4_pairs.png")
pairs(biochemData_d4,main="DAY 4")
dev.off()

# Enlarging dataset for metabolomics data comparison

biochemDataForMetab<-biochemData[rep(1:nrow(biochemData1),each=2),] # duplicating each row to avvount for technical replicates
rownames(biochemDataForMetab)<-gsubfn(".", list("a"="b1_r001", "a.1"="b1_r002", "b"="b2_r001","b.1"="b2_r002","c"="b3_r001","c.1"="b3_r002"),rownames(biochemDataForMetab))
rownames(biochemDataForMetab)<-gsub("r001.1","r002",rownames(biochemDataForMetab))

# Day12 
batch_corrected_mat_d12<-svd_day12_nonzero[[4]]
batch_corrected_mat_d12<-batch_corrected_mat_d12-min(batch_corrected_mat_d12)
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

a<-adonis(t(batch_corrected_mat_d12)~ biochemDataForMetab_d12$biomass+biochemDataForMetab_d12$biomassProduc+biochemDataForMetab_d12$lipidProduc+
            biochemDataForMetab_d12$totalLipidCon+biochemDataForMetab_d12$totalProtein+biochemDataForMetab_d12$totalCarbCon+biochemDataForMetab_d12$Strain, 
          method = "euclidean", perm=999)

# Call:
#   adonis(formula = t(batch_corrected_mat_d12) ~ biochemDataForMetab_d12$biomass +      biochemDataForMetab_d12$biomassProduc + biochemDataForMetab_d12$lipidProduc +      biochemDataForMetab_d12$totalLipidCon + biochemDataForMetab_d12$totalProtein +      biochemDataForMetab_d12$totalCarbCon + biochemDataForMetab_d12$Strain,      permutations = 999, method = "bray") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
# biochemDataForMetab_d12$biomass         1  0.005050 0.0050503 10.9537 0.04364  0.001 ***
#   biochemDataForMetab_d12$biomassProduc   1  0.002892 0.0028916  6.2717 0.02499  0.001 ***
#   biochemDataForMetab_d12$lipidProduc     1  0.002903 0.0029026  6.2956 0.02508  0.001 ***
#   biochemDataForMetab_d12$totalLipidCon   1  0.004190 0.0041902  9.0881 0.03621  0.001 ***
#   biochemDataForMetab_d12$totalProtein    1  0.003249 0.0032493  7.0475 0.02808  0.001 ***
#   biochemDataForMetab_d12$totalCarbCon    1  0.002412 0.0024119  5.2312 0.02084  0.001 ***
#   biochemDataForMetab_d12$Strain         21  0.053525 0.0025488  5.5281 0.46255  0.001 ***
#   Residuals                              90  0.041495 0.0004611         0.35860           
# Total                                 117  0.115716                   1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

a<-adonis(t(batch_corrected_mat_d12)~ biochemDataForMetab_d12$biomass+biochemDataForMetab_d12$lipidProduc+
            biochemDataForMetab_d12$totalLipidCon+biochemDataForMetab_d12$totalProtein+biochemDataForMetab_d12$totalCarbCon+biochemDataForMetab_d12$Strain, 
          method = "bray", perm=999)

# Call:
#   adonis(formula = t(batch_corrected_mat_d12) ~ biochemDataForMetab_d12$biomass +      biochemDataForMetab_d12$lipidProduc + biochemDataForMetab_d12$totalLipidCon +      biochemDataForMetab_d12$totalProtein + biochemDataForMetab_d12$totalCarbCon +      biochemDataForMetab_d12$Strain, permutations = 999, method = "bray") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
# biochemDataForMetab_d12$biomass         1  0.005050 0.0050503 10.9302 0.04364  0.001 ***
#   biochemDataForMetab_d12$lipidProduc     1  0.003198 0.0031982  6.9218 0.02764  0.001 ***
#   biochemDataForMetab_d12$totalLipidCon   1  0.003937 0.0039371  8.5209 0.03402  0.001 ***
#   biochemDataForMetab_d12$totalProtein    1  0.003569 0.0035693  7.7250 0.03085  0.001 ***
#   biochemDataForMetab_d12$totalCarbCon    1  0.002367 0.0023668  5.1224 0.02045  0.001 ***
#   biochemDataForMetab_d12$Strain         21  0.055548 0.0026451  5.7248 0.48004  0.001 ***
#   Residuals                              91  0.042047 0.0004621         0.36336           
# Total                                 117  0.115716                   1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

permutedData<-batch_corrected_mat_d12[,sample(ncol(batch_corrected_mat_d12))]
a<-adonis(t(permutedData)~ biochemDataForMetab_d12$biomass+biochemDataForMetab_d12$lipidProduc+
            biochemDataForMetab_d12$totalLipidCon+biochemDataForMetab_d12$totalProtein+biochemDataForMetab_d12$totalCarbCon+biochemDataForMetab_d12$Strain, 
          method = "bray", perm=999)
# Call:
#   adonis(formula = t(permutedData) ~ biochemDataForMetab_d12$biomass +      biochemDataForMetab_d12$lipidProduc + biochemDataForMetab_d12$totalLipidCon +      biochemDataForMetab_d12$totalProtein + biochemDataForMetab_d12$totalCarbCon +      biochemDataForMetab_d12$Strain, permutations = 999, method = "bray") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)
# biochemDataForMetab_d12$biomass         1  0.001204 0.00120376 1.21556 0.01040  0.197
# biochemDataForMetab_d12$lipidProduc     1  0.001167 0.00116682 1.17825 0.01008  0.240
# biochemDataForMetab_d12$totalLipidCon   1  0.000628 0.00062755 0.63370 0.00542  0.973
# biochemDataForMetab_d12$totalProtein    1  0.000767 0.00076654 0.77405 0.00662  0.805
# biochemDataForMetab_d12$totalCarbCon    1  0.000789 0.00078910 0.79683 0.00682  0.786
# biochemDataForMetab_d12$Strain         21  0.021046 0.00100217 1.01199 0.18187  0.398
# Residuals                              91  0.090117 0.00099029         0.77877       
# Total                                 117  0.115716                    1.00000  

# non-batch effect removed data

#ms_data_day12_nonzero

non_batch_corrrected_mat_d12<-log(ms_data_day12_nonzero)
colnames(non_batch_corrrected_mat_d12)<-gsub("_b","b",colnames(non_batch_corrrected_mat_d12))
biochemDataForMetab_d12_nb<-biochemDataForMetab[colnames(non_batch_corrrected_mat_d12),]
biochemDataForMetab_d12_nb<-biochemDataForMetab_d12_nb[complete.cases(biochemDataForMetab_d12_nb),]
biochemDataForMetab_d12_nb$Strain<-gsub("b1_r001|b1_r002|b2_r001|b2_r002|b3_r001|b3_r002","",rownames(biochemDataForMetab_d12_nb))

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
batch_corrected_mat_d4<-batch_corrected_mat_d4-min(batch_corrected_mat_d4)
#batch_corrected_mat_d4<-batch_corrected_mat_d4[[7]] #removing last column which has wrong pvalues
#batch_corrected_mat_d4<-batch_corrected_mat_d4[,1:125] #removing last column which has pvalues
colnames(batch_corrected_mat_d4)<-gsub("_b","b",colnames(batch_corrected_mat_d4))
colnames(batch_corrected_mat_d4)<-gsub('D4_14','D4_014',colnames(batch_corrected_mat_d4))

biochemDataForMetab_d4<-biochemDataForMetab[colnames(batch_corrected_mat_d4),]
biochemDataForMetab_d4<-biochemDataForMetab_d4[complete.cases(biochemDataForMetab_d4),]
biochemDataForMetab_d4$Strain<-gsub("b1_r001|b1_r002|b2_r001|b2_r002|b3_r001|b3_r002","",rownames(biochemDataForMetab_d4))

#sum(is.na(biochemDataForMetab_d4))

a<-adonis(t(batch_corrected_mat_d4)~ biochemDataForMetab_d4$biomass+biochemDataForMetab_d4$biomassProduc+biochemDataForMetab_d4$lipidProduc+
            biochemDataForMetab_d4$totalLipidCon+biochemDataForMetab_d4$totalProtein+biochemDataForMetab_d4$totalCarbCon+biochemDataForMetab_d4$Strain, 
          method = "bray", perm=999)

# Call:
#   adonis(formula = t(batch_corrected_mat_d4) ~ biochemDataForMetab_d4$biomass +      biochemDataForMetab_d4$biomassProduc + biochemDataForMetab_d4$lipidProduc +      biochemDataForMetab_d4$totalLipidCon + biochemDataForMetab_d4$totalProtein +      biochemDataForMetab_d4$totalCarbCon + biochemDataForMetab_d4$Strain,      permutations = 999, method = "bray") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)    
# biochemDataForMetab_d4$biomass         1  0.002207 0.00220749  4.8525 0.02440  0.001 ***
#   biochemDataForMetab_d4$biomassProduc   1  0.001792 0.00179246  3.9402 0.01981  0.001 ***
#   biochemDataForMetab_d4$lipidProduc     1  0.001703 0.00170255  3.7425 0.01882  0.001 ***
#   biochemDataForMetab_d4$totalLipidCon   1  0.002365 0.00236546  5.1997 0.02614  0.001 ***
#   biochemDataForMetab_d4$totalProtein    1  0.001229 0.00122885  2.7012 0.01358  0.001 ***
#   biochemDataForMetab_d4$totalCarbCon    1  0.001276 0.00127556  2.8039 0.01410  0.001 ***
#   biochemDataForMetab_d4$Strain         21  0.035782 0.00170388  3.7455 0.39546  0.001 ***
#   Residuals                             97  0.044127 0.00045492         0.48770           
# Total                                124  0.090481                    1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

a<-adonis(t(batch_corrected_mat_d4)~ biochemDataForMetab_d4$biomass+biochemDataForMetab_d4$lipidProduc+
            biochemDataForMetab_d4$totalLipidCon+biochemDataForMetab_d4$totalProtein+biochemDataForMetab_d4$totalCarbCon+biochemDataForMetab_d4$Strain, 
          method = "bray", perm=999)

# Call:
#   adonis(formula = t(batch_corrected_mat_d4) ~ biochemDataForMetab_d4$biomass +      biochemDataForMetab_d4$lipidProduc + biochemDataForMetab_d4$totalLipidCon +      biochemDataForMetab_d4$totalProtein + biochemDataForMetab_d4$totalCarbCon +      biochemDataForMetab_d4$Strain, permutations = 999, method = "bray") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
# biochemDataForMetab_d4$biomass         1  0.002207 0.0022075  4.8209 0.02440  0.001 ***
#   biochemDataForMetab_d4$lipidProduc     1  0.001787 0.0017867  3.9020 0.01975  0.001 ***
#   biochemDataForMetab_d4$totalLipidCon   1  0.001682 0.0016820  3.6733 0.01859  0.001 ***
#   biochemDataForMetab_d4$totalProtein    1  0.001161 0.0011610  2.5356 0.01283  0.001 ***
#   biochemDataForMetab_d4$totalCarbCon    1  0.001543 0.0015427  3.3690 0.01705  0.001 ***
#   biochemDataForMetab_d4$Strain         21  0.037228 0.0017727  3.8715 0.41144  0.001 ***
#   Residuals                             98  0.044874 0.0004579         0.49595           
# Total                                124  0.090481                   1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

########### CCA

cca.1 <- cca(t(batch_corrected_mat_d4)~ biochemDataForMetab_d4$biomass+biochemDataForMetab_d4$lipidProduc+
               biochemDataForMetab_d4$totalLipidCon+biochemDataForMetab_d4$totalProtein+biochemDataForMetab_d4$totalCarbCon+biochemDataForMetab_d4$Strain)

########## Differential feature identification

### Comparisons

Chlorella<-c("D12_001","D12_006","D12_051","D12_104","D12_14","D12_177","D12_187","D12_207","D12_268","D12_283","D12_322","D12_325","D12_84","D12_87","D12_94")
Parachlorella<-c("D12_245","D12_252","D12_253","D12_254","D12_255","D12_258")
UnIdchlorella<-"D12_184"

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
SampleGroup

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
#updated on 12-May-2014 by shiv to include p.adjust method. Thus, using 'corr.test' from psych package.

correlationPlot<-function(env_parameters){
  env_parameters<-biochemData1
  combos <- combn(ncol(env_parameters),2) # combinations without repetitions
  
  #combos <- expand.grid(rep(list(1:ncol(env_parameters)), 2 )) # combinations with repetitions
  #combos <- as.matrix(combos)
  #combos <- t(combos) # transpose matrix
  
  mat1 <- adply(combos, 2, function(x) {
    #test <- chisq.test(env_parameters[, x[1]], env_parameters[, x[2]])
    test<-corr.test(as.data.frame(env_parameters[, x[1]]), as.data.frame(env_parameters[, x[2]]), use="all.obs", method="kendall", adjust="bonferroni") 
    
    out <- data.frame("Row" = colnames(env_parameters)[x[1]]
                      , "Column" = colnames(env_parameters[x[2]])
                      , "p.value" = as.numeric(round(test$p, 3))
                      ,  "tau.value" = as.numeric(round(test$r,3))
    )
    return(out)
  })
  
  tau.valu.mat<-mat1
  # mat 2 is the same table but formated to show pvalue of corrrelations
  mat2 <- adply(combos, 2, function(x) {
    #test <- chisq.test(env_parameters[, x[1]], env_parameters[, x[2]])
    test<-corr.test(as.data.frame(env_parameters[, x[1]]), as.data.frame(env_parameters[, x[2]]), use="all.obs", method="kendall", adjust="bonferroni") 
    
    out <- data.frame("Row" = colnames(env_parameters)[x[1]]
                      , "Column" = colnames(env_parameters[x[2]])
                      , "tau.value" = as.numeric(round(test$r,3))
                      , "p.value" = as.numeric(round(test$p, 3))
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