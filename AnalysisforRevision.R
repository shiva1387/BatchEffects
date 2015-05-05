# Data Analysis in R-for addressing reviewers comments#
#######################################################
# Author(s): Shiv
# Version: 26032015 
# Input: ".txt" file from Kirwan et al(2014), Dataset 8A from MTBLS79 (MetaboLights Accession) 

#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 
#source("http://bioconductor.org/biocLite.R")
library(ggplot2)
library(data.table)
library(stringr)
library(vegan)
library(gplots)
library(raster)
library(RColorBrewer)
library(reshape2)
library(preprocessCore)
library(multtest)

#############
# User      #
# Specific  #
# Variables #
#############

directory<- "../data/kirwan2014/"
#Contains path for .tsv file
# The "PATH", Remember to use "/" and not "/" in the path
setwd(directory)
getwd()

#################
# Get filenames #
#################

mzfilename<-"Kirwan_Dataset08a_SFPM_PQN.txt"
ms_data_total<-read.table(mzfilename,sep="\t",header=T,check.names=FALSE,row.names=1)
ms_data_total<-as.data.frame(t(ms_data_total)) #samples as columns and mz as rows
str(ms_data_total)

metaData<-read.table("Kirwan_MetaInfo08a_SFPM_PQN.txt",sep="\t",header=T,check.names=FALSE)

#######################################
###### Exploratory data analysis ######
#######################################

## Functions

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

ScaleData<-function(data_matrix){
  processed_data<-scale(data_matrix,center=T,scale=T)
  colnames(processed_data)<-colnames(data_matrix)
  rownames(processed_data)<-rownames(data_matrix)
  return(processed_data)
}
normalizeData<-function(data_matrix){
  processed_data<-normalize.quantiles(as.matrix(data_matrix),copy=TRUE)
  colnames(processed_data)<-colnames(data_matrix)
  rownames(processed_data)<-rownames(data_matrix)
  return(processed_data)
}
compute_pca<-function(dataset){#,preprocess_method) {
  dataset<-as.matrix(log1p(dataset)) #log transform the data using natural log
processed_data<-dataset #as the data has already been normalized by kirwan
  #   if(preprocess_method=="norm") {
#     processed_data<-normalize.quantiles(as.matrix(dataset),copy=TRUE)
#   }   else  {
#     processed_data<-scale(dataset,center=T,scale=T)
#     processed_data<-processed_data-min(processed_data)
#   }
  pca_results <- princomp(processed_data,cor=F,scores=T) ### IMP: choose quantile normalized or scaled data
  return(pca_results)
}

# for only a single dataframe such as those of uncorrected data
compute_perm_ftest<-function(dataset,classlabel_factor) {
  data_matrix<-as.matrix(dataset)
  dataset_sig_features<-mt.maxT(data_matrix,classlabel_factor,test="f",side="abs",fixed.seed.sampling="y",B=1000,nonpara="n")
  p.values<-dataset_sig_features
  id.sig_dataset <- sort(dataset_sig_features[dataset_sig_features$adjp < 0.05,c(1)])
  metab_sig<-cbind(data_matrix[id.sig_dataset,],round(dataset_sig_features$adjp[dataset_sig_features$index %in% id.sig_dataset],5))
  return(list(p.values,metab_sig))
}

# for only a single dataframe such as those of uncorrected data
compute_r2_from_fvalue_singleDataFrame<-function(dataset,classlabel) {
  classlabel_factor<-as.numeric(as.factor(classlabel))-1
  fval_teststat<-as.data.frame(matrix(NA,nrow(dataset),1))
  r2_teststat<-as.data.frame(matrix(NA,nrow(dataset),1))
  data_matrix<-as.matrix(dataset)
  dataset_fvalue<-mt.teststat(data_matrix,classlabel_factor,test="f",nonpara="n")
  dataset_r2<-r2.from.Fstat(dataset_fvalue,length(unique(classlabel)),ncol(data_matrix))
  fval_teststat[,1]<-dataset_fvalue
  r2_teststat[,1]<-dataset_r2
  return (list(fval_teststat,r2_teststat))
}

compute_linearModel_associations<-function(results.from.pca,batch) {
  lm_pca_scores<-apply(results.from.pca$loadings,2, function(x) {
    lm_val<-lm(x~ as.factor(batch) )#
    lm_cor<-summary(lm_val)
    p.val_model<-anova(lm_val)$'Pr(>F)'[1]
    fvalue_model<-anova(lm_val)$'F value'[1] # modified by shiv to return r2 for each model term on 09-Jan 2015
    r2_model<-anova(lm_val)$'Sum Sq'/sum(anova(lm_val)$'Sum Sq') #includes the r2 term for residuals
    r2_model<-r2_model[1] #includes the r2 term only for runday and runday/strain
    return(list(lm_cor$r.squared,p.val_model,fvalue_model,r2_model))
  })
}

extract.variables.pc.associations<-function(linearmodel_list,variable.extract) {
  if(variable.extract=="r2model") { #WARNING:code implictly assumes r2 is in the first column and p.val in the second
    return (sapply(linearmodel_list, function(x){as.numeric(x[1])}))
  } else if (variable.extract=="pvalue"){
    return (sapply(linearmodel_list, function(x){x[2]}))
  } else if (variable.extract=="fvalue"){
    return (sapply(linearmodel_list, function(x){x[3]}))
  } else {
    return (sapply(linearmodel_list, function(x){x[4]}))
  }
}


#### PCA
ms_data_samples<- ms_data_total[, -grep('QC', colnames(ms_data_total))] #removing QC for calculating associations with biological samples only

dataForPlotting<-ms_data_total# ms_data_samples #ms_data_total #used for batches
metainfoForPlotting<-metaData

ms_data_princomp<-compute_pca(dataForPlotting) 
residual_variance<-ms_data_princomp$sdev^2/sum(ms_data_princomp$sdev^2)

plot(ms_data_princomp$loadings[,2]~ms_data_princomp$loadings[,1])
text(ms_data_princomp$loadings[,2]~ms_data_princomp$loadings[,1], labels = colnames(dataForPlotting), cex=0.6, pos=4)

### PCA ggplot
forPlot<-data.frame(PCaxisA = ms_data_princomp$loadings[,1],PCaxisB = ms_data_princomp$loadings[,3], 
                    dataType=metainfoForPlotting$Sample_Type, batch=metainfoForPlotting$Batches_forBatchCorrectionAlgorithm)

pch_types<-c(15, 16, 17, 18)
pch_values<-sample(pch_types, 4)

##plotting
plot1<- ggplot(data=forPlot, aes(x=PCaxisA, y=PCaxisB, colour= factor(batch), shape = factor(dataType))) + geom_point(size=4)
#+facet_wrap(dataType~GrowthStage_f)#for samples
plot2<- plot1+ scale_shape_manual('GrowthStage',values=pch_values) +  ## For single datatype, shapes based on timepoint
  theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1,size=12),axis.text.y=element_text(size=12),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black'))+ 
  xlab(paste0("PC 1 loadings","\n","Variation exp= ",round(residual_variance[1]*100,2),"%")) + 
  ylab(paste0("PC 3 loadings","\n","Variation exp= ",round(residual_variance[3]*100,2),"%")) +
  ggtitle("PCA Samples - Kirwan et al")


pdf("pca_pc13_samples_Kirwan.pdf",height=12,width=12)
plot2
dev.off()

### Calculating PC associations

#Strain
ms_data_samples<- ms_data_total[, -grep('QC', colnames(ms_data_total))] #removing QC for calculating associations with biological samples only
ms_data_princomp_strain<-compute_pca(ms_data_samples) 
residual_variance_strain<-ms_data_princomp$sdev^2/sum(ms_data_princomp$sdev^2)

metainfoForPlotting_samples<-metainfoForPlotting[-grep('QC',metainfoForPlotting$SampleName),]

lm_pca_loadings_strain<-compute_linearModel_associations(ms_data_princomp_strain, metainfoForPlotting_samples$Class)
lm_pca_loadings_r.sq.model_strain<-extract.variables.pc.associations(lm_pca_loadings_strain,"r2model")
lm_pca_loadings_pval_strain<-extract.variables.pc.associations(lm_pca_loadings_strain,"pvalue")
lm_pca_loadings_pval_strain<-do.call(rbind.data.frame, lm_pca_loadings_pval_strain)
lm_pca_loadings_pval_strain<-as.data.frame(sapply(lm_pca_loadings_pval_strain, function(x) p.adjust(as.numeric(x),method="BH"))) # FDR correction

lm_pca_loadings_fval_strain<-extract.variables.pc.associations(lm_pca_loadings_strain,"fvalue")
lm_pca_loadings_fval_strain<-do.call(rbind.data.frame, lm_pca_loadings_fval_strain)
lm_pca_loadings_r2_strain<-extract.variables.pc.associations(lm_pca_loadings_strain,"r2")
lm_pca_loadings_r2_strain<-do.call(rbind.data.frame, lm_pca_loadings_r2_strain)

r2.pval.dataset_strain<-cbind(as.numeric(paste0(1:length(lm_pca_loadings_r.sq.model_strain))),lm_pca_loadings_pval_strain, 
                              lm_pca_loadings_r2_strain)

colnames(r2.pval.dataset_strain)<-c("PrincipalComponents","P-value(strain)","R2(batch)")
write.table(r2.pval.dataset_strain,"PC_associations_strain.txt",quote=FALSE,row.names=FALSE,sep="\t")

#Batch
#lm_pca_loadings<-compute_linearModel_associations(ms_data_princomp, metainfoForPlotting_samples$Batches_forBatchCorrectionAlgorithm) #no QC
lm_pca_loadings<-compute_linearModel_associations(ms_data_princomp, metainfoForPlotting$Batches_forBatchCorrectionAlgorithm) #with QC

lm_pca_loadings_r.sq.model<-extract.variables.pc.associations(lm_pca_loadings,"r2model")
lm_pca_loadings_pval<-extract.variables.pc.associations(lm_pca_loadings,"pvalue")
lm_pca_loadings_pval<-do.call(rbind.data.frame, lm_pca_loadings_pval)
lm_pca_loadings_pval<-as.data.frame(sapply(lm_pca_loadings_pval, function(x) p.adjust(as.numeric(x),method="BH"))) # FDR correction

lm_pca_loadings_fval<-extract.variables.pc.associations(lm_pca_loadings,"fvalue")
lm_pca_loadings_fval<-do.call(rbind.data.frame, lm_pca_loadings_fval)
lm_pca_loadings_r2<-extract.variables.pc.associations(lm_pca_loadings,"r2")
lm_pca_loadings_r2<-do.call(rbind.data.frame, lm_pca_loadings_r2)

r2.pval.dataset<-cbind(as.numeric(paste0(1:length(lm_pca_loadings_r.sq.model))),lm_pca_loadings_pval, lm_pca_loadings_r2)

colnames(r2.pval.dataset)<-c("PrincipalComponents","P-value(batch)","R2(batch)")
write.table(r2.pval.dataset,"PC_associations.txt",quote=FALSE,row.names=FALSE,sep="\t")
#### Add a column called Significant to mark significant R2 values

#Strain
r2.pval.dataset_forPlot_strain<-r2.pval.dataset_strain
r2.pval.dataset_forPlot_strain$"Sig-R2(strain)"<-rep("Sig-R2(strain)",nrow(r2.pval.dataset_forPlot_strain))
r2.pval.dataset_forPlot_strain[r2.pval.dataset_forPlot_strain[,2]> 0.05,4]<-"NS-R2(strain)"
## Nasty workaround to obtain dataset in a format compatible with ggplot
r2.pval.dataset_forPlot_dataType_strain<-r2.pval.dataset_forPlot_strain[,c(1,3,4)]
r2.pval.dataset_forPlot_dataType_strain$"data"<-rep("Strain",nrow(r2.pval.dataset_forPlot_dataType_strain))
colnames(r2.pval.dataset_forPlot_dataType_strain)<-c("PrincipalComponents","R2","SigR2","data")

#Batch
r2.pval.dataset_forPlot<-r2.pval.dataset
r2.pval.dataset_forPlot$"Sig-R2(batch)"<-rep("Sig-R2(batch)",nrow(r2.pval.dataset_forPlot))
r2.pval.dataset_forPlot[r2.pval.dataset_forPlot[,2]> 0.05,4]<-"NS-R2(batch)"

## Nasty workaround to obtain dataset in a format compatible with ggplot
r2.pval.dataset_forPlot_dataType<-r2.pval.dataset_forPlot[,c(1,3,4)]
r2.pval.dataset_forPlot_dataType$"data"<-rep("RunDay",nrow(r2.pval.dataset_forPlot_dataType))
colnames(r2.pval.dataset_forPlot_dataType)<-c("PrincipalComponents","R2","SigR2","data")

# r2.pval.dataset_forPlot_ResidualVariance<-data.frame(as.numeric(paste0(1:length(lm_pca_loadings_r.sq.model))),
#                                                      round(residual_variance*100,2),
#                                                      rep("ResidualVariance",length(lm_pca_loadings_r.sq.model)))
# r2.pval.dataset_forPlot_ResidualVariance$"data"<-rep("ResidualVariance",nrow(r2.pval.dataset_forPlot_ResidualVariance))
# colnames(r2.pval.dataset_forPlot_ResidualVariance)<-c("PrincipalComponents","R2","SigR2","data")
# 
# r2.pval.dataset_forPlot_merged<-rbind(r2.pval.dataset_forPlot_dataType,r2.pval.dataset_forPlot_ResidualVariance)

r2.pval.dataset_forPlot_merged<-rbind(r2.pval.dataset_forPlot_dataType,r2.pval.dataset_forPlot_dataType_strain)

###################### 

# cols = gg_color_hue(8)

# r2.pval.dataset_forPlot_merged$"SigR2" = factor(r2.pval.dataset_forPlot_merged$"SigR2", 
#                                                 levels=c("Sig-R2(batch)","NS-R2(batch)",
#                                                          "ResidualVariance"))
# r2.pval.dataset_forPlot_merged$data = factor(r2.pval.dataset_forPlot_merged$data, 
#                                               levels=c("R2","ResidualVariance"))

r2.pval.dataset_forPlot_merged$"SigR2" = factor(r2.pval.dataset_forPlot_merged$"SigR2", 
                                                 levels=c("Sig-R2(strain)","NS-R2(strain)",
                                                          "Sig-R2(batch)","NS-R2(batch)"))
  
#Plotting strain and batch associations
set.seed(1)
# ggplot
plot1<- ggplot(data=r2.pval.dataset_forPlot_merged, aes(x=log(PrincipalComponents), y=R2,  group=data,
                                                        colour= SigR2, shape = SigR2)) +
  geom_point(size=2) +  geom_line(size=0.2) +   scale_colour_manual(values=c("#7CAE00", "#7CAE00","#F8766D", "#F8766D")) + 
  scale_shape_manual(values=c(16,1,17,2)) + #facet_grid(data~.,scales="free_y") +
  theme_bw() + theme(axis.text.x=element_text(hjust = 1,size=8),axis.text.y=element_text(size=8),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     #panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black'))

#ggsave("PC_associations_kirwan_batch_strain_noQC_log.pdf",plot1,height=8,width=10)
ggsave("PC_associations_kirwan_batchwithQC_strain_log.pdf",plot1,height=8,width=10)

##Only plotting batch and run day
# ## ggplot
# set.seed(1)
# # ggplot 
# plot1<- ggplot(data=r2.pval.dataset_forPlot_merged, aes(x=PrincipalComponents, y=R2,  group=data,
#                                                         colour= SigR2, shape = SigR2)) +
#   geom_point(size=2) +  geom_line(size=0.2) + facet_grid(data~.,scales="free_y") +
#   theme_bw() + theme(axis.text.x=element_text(hjust = 1,size=8),axis.text.y=element_text(size=8),
#                      panel.grid.major.x = element_blank(), # to x remove gridlines
#                      panel.grid.major.y = element_blank(), # to y remove gridlines
#                      #panel.border = element_blank(),  # remove top and right border
#                      panel.background = element_blank(),
#                      axis.line = element_line(color = 'black'))
# 
# ggsave("PC_associations_kirwan.pdf",plot1,height=8,width=10)

## SVD

svd_dataset<-svd(ms_data_total)

# Removing PC1
svd_dataset$d1<-svd_dataset$d
svd_dataset$d1[1]<-0 #specifying the pc component to be made zero
dataset_rmbatch1<-svd_dataset$u %*% diag(svd_dataset$d1) %*% t(svd_dataset$v) #Removing the variation in start_pc_comp.
dataset_rmbatch1<-dataset_rmbatch1-min(dataset_rmbatch1) #To ensure values are above zero
rownames(dataset_rmbatch1)<-rownames(ms_data_total)
colnames(dataset_rmbatch1)<-colnames(ms_data_total)

#Removing PC1,3
svd_dataset$d1<-svd_dataset$d
svd_dataset$d1[c(1,3)]<-0 #specifying the pc component to be made zero
dataset_rmbatch2<-svd_dataset$u %*% diag(svd_dataset$d1) %*% t(svd_dataset$v) #Removing the variation in start_pc_comp.
dataset_rmbatch2<-dataset_rmbatch2-min(dataset_rmbatch2)
rownames(dataset_rmbatch2)<-rownames(ms_data_total)
colnames(dataset_rmbatch2)<-colnames(ms_data_total)

#Removing PC1,3,4
svd_dataset$d1<-svd_dataset$d
svd_dataset$d1[c(1,3,4)]<-0 #specifying the pc component to be made zero
dataset_rmbatch3<-svd_dataset$u %*% diag(svd_dataset$d1) %*% t(svd_dataset$v) #Removing the variation in start_pc_comp.
dataset_rmbatch3<-dataset_rmbatch3-min(dataset_rmbatch3)
rownames(dataset_rmbatch3)<-rownames(ms_data_total)
colnames(dataset_rmbatch3)<-colnames(ms_data_total)

#Removing PC1,3,4,6,7
svd_dataset$d1<-svd_dataset$d
svd_dataset$d1[c(1,3,4,6,7)]<-0 #specifying the pc component to be made zero
dataset_rmbatch5<-svd_dataset$u %*% diag(svd_dataset$d1) %*% t(svd_dataset$v) #Removing the variation in start_pc_comp.
dataset_rmbatch5<-dataset_rmbatch5-min(dataset_rmbatch5)
rownames(dataset_rmbatch5)<-rownames(ms_data_total)
colnames(dataset_rmbatch5)<-colnames(ms_data_total)


#### Differential analysis

#permutational analysis
classlabel<-as.numeric(metaData$Batches_forBatchCorrectionAlgorithm)-1
ms_data_total_rna_perm<-compute_perm_ftest(dataset_rmbatch5, classlabel)
ms_data_total_rna_perm_pval<-ms_data_total_rna_perm[[1]]
ms_data_total_rna_perm_pval$BHadj<-p.adjust(as.numeric(ms_data_total_rna_perm_pval$rawp),method="BH")
ms_data_total_rna_perm_matrix<-ms_data_total_rna_perm[[2]]
dim(ms_data_total_rna_perm_matrix)
# 814 173 #uncorrected
# 756 173 #removing PC1
# 649 173 #removing PC1,3
# 515 173 #removing PC1,3,4
# 163 173 #removing PC1,3,4,6,7

##t-test
dataset_rmbatch_samples<- dataset_rmbatch5[, -grep('QC', colnames(dataset_rmbatch5))] #removing Blanks

# Comparing C with S
dataset_rmbatch_samples_C <- dataset_rmbatch_samples[, grep('_C', colnames(dataset_rmbatch_samples))] 
dataset_rmbatch_samples_S <- dataset_rmbatch_samples[, grep('_S', colnames(dataset_rmbatch_samples))] 
dataset_rmbatch_samples_CS<-cbind(dataset_rmbatch_samples_C,dataset_rmbatch_samples_S)

a.list <- apply(dataset_rmbatch_samples_CS,1,function(x){t.test(x[1:66],x[67:134])$p.value})
a.list <-p.adjust(as.numeric(a.list),method="BH")
id.sig <- which(a.list < 0.05 );
metab.sig_CS<-dataset_rmbatch_samples_CS[id.sig,]
nrow(metab.sig_CS)
#1709 #1629 after BH correction#uncorrected
#1591 #removing PC1 after BH correction
#1609 #removing PC1,3 after BH correction
#1600 #removing PC1,3,4 after BH correction
#1635 #removing PC1,3,4,6,7 after BH correction

# > a<-intersect(rownames(ms_data_total_rna_perm_matrix),rownames(metab.sig_CS))
# > length(a)
# [1] 52 #52 features common between batch associated variables and biologically significant variables
#[1] 42 #when using BH corrected t.test results

##Comparing RSD (median) values
rsd<-function(x){
  rsd_value<-(sd(x)/mean(as.numeric(x)))*100
  return(rsd_value)
}

# For QC's
dataToCompare<-ms_data_total # To ensure all values are in the positive range
ms_data_qc<- dataToCompare[, grep('QC', colnames(dataToCompare))] #removing Blanks
dt <- data.table(t(ms_data_qc))
dt$batch<-metaData[grep('QC',metaData$SampleName),4]
ms_data_qc_rsd<-t(as.data.frame(dt[,lapply(.SD,rsd),by=batch]))
ms_data_qc_rsd[is.na(ms_data_qc_rsd)] <- 0
ms_data_qc_rsd1<-apply(ms_data_qc_rsd[2:nrow(ms_data_qc_rsd),],2,median)#first row is batch information
ms_data_qc_rsd1<-apply(ms_data_qc_rsd,2,median)#first row is batch information

###QC's
#9.462764 10.028034  7.656542 10.317175  9.025295 11.751496  9.231841  8.801779 #Uncorrected
#0.008832475 0.009711961 0.007457190 0.008953825 0.008534387 0.011082341 0.009639492 0.009325934 #removing PC1
#0.01094398 0.01137419 0.00948960 0.01134566 0.01109268 0.01323674 0.01145353 0.01131700 #removing PC1,3
#0.01094398 0.01137419 0.00948960 0.01134566 0.01109268 0.01323674 0.01145353 0.01131700 #removing PC1,3,4
#0.011937755 0.012861325 0.009836365 0.011690704 0.010437756 0.014580289 0.011612906 0.014955328 #removing PC1,3,4,6,7

#For Biological samples

dataToCompare<-ms_data_total# To ensure all values are in the positive range)
ms_data_bio<- dataToCompare[, grep('_S01', colnames(dataToCompare))] #Selecting cows
#ms_data_bio<- dataToCompare[, grep('_C01', colnames(dataToCompare))] #Selecting sheep
dt <- data.table(t(ms_data_bio))
dt$batch<-metaData[grep('_S01',metaData$SampleName),5]
ms_data_bio_rsd<-as.numeric(dt[,lapply(.SD,rsd),by=batch])
ms_data_bio_rsd[is.na(ms_data_bio_rsd)] <- 0
ms_data_bio_rsd1<-median(ms_data_bio_rsd[2:2489])
ms_data_bio_rsd1

#Cow1
#19.0 #Uncorrected
#0.02 #removing PC1,3,4,6,7
#Cow2
#19.16 #Uncorrected
#0.19 #removing PC1,3,4,6,7
#Cow3
#19.55 #Uncorrected
#0.02 #removing PC1,3,4,6,7
#Cow4
#17.94 #Uncorrected
#0.02 #removing PC1,3,4,6,7
#Cow5
#19.53 #Uncorrected
#0.02 #removing PC1,3,4,6,7
#Cow6
#14.42 #Uncorrected
#0.02 #removing PC1,3,4,6,7
#Cow7
#21.92 #Uncorrected
#0.02 #removing PC1,3,4,6,7
#Cow8
#23.4 #Uncorrected
#0.02 #removing PC1,3,4,6,7
#Cow9
#18.2 #Uncorrected
#0.02 #removing PC1,3,4,6,7
#Cow10
#22.58 #Uncorrected
#0.02 #removing PC1,3,4,6,7


#Sheep

#Sheep1
#18.28 #Uncorrected
#0.01 #removing PC1,3,4,6,7
#Sheep2
#19.60 #Uncorrected
#0. #02removing PC1,3,4,6,7
#Sheep3
#21.24 #Uncorrected
#0.19 #removing PC1,3,4,6,7
#Sheep4
#19.16 #Uncorrected
#0.02 #removing PC1,3,4,6,7
#Sheep5
#21.65 #Uncorrected
#0.02 #removing PC1,3,4,6,7
#Sheep6
#21.07 #Uncorrected
#0.01 #removing PC1,3,4,6,7
#Sheep7
#19.20 #Uncorrected
#0.02 #removing PC1,3,4,6,7
#Sheep8
#22.93 #Uncorrected
#0.03 #removing PC1,3,4,6,7
#Sheep9
#20.12 #Uncorrected
#0.03 #removing PC1,3,4,6,7
#Sheep10
#19.38 #Uncorrected
#0.03 #removing PC1,3,4,6,7
