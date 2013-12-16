###########################################
# Data Analysis in R-Malaysian algae data #
###########################################
# Author(s): Shiv
# Version: 16122013 
# Input: ".tsv" file from XCMS 
# Software: XCMS
# Modified By :Shivshankar Umashankar
# Modified server version for xcms file "algae_4setblanks_021213.tsv" produced using xcms_1.38
##############################################
# Data,software used for obaining input data #
##############################################
# The files, algae_417reps_sorted.txt,Algae_22strains_indreps_blanks_matrix.tsv,algae_4setblanks_060913.tsv are the same.
# 
# algae_417reps_sorted.txt and algae_4setblanks_060913 are the same,the difference is in algae_417reps_sorted.txt, 
# the columns are sorted based on column names
# 
# The difference between Algae_22strains_indreps_blanks_matrix and the other 2 files is, 
# Algae_22strains_indreps_blanks_matrix contains only the 22 strains + blanks and matrix. 
# It does not contain any duplicate columns.
# 
# After numberous xcms parameter optimization, the below parameters provided the best resolution and features. 
# This was determined by looking at the mz/rt list for the blanks and sample strain and comparing the counts/TIC with the raw
# MS file using Mass Hunter(Agilent software to display raw data from the instrument). 
# Features were checked to see if they were real or noise picked up by XCMS.
# The known features were checked for accuracy in peak characteristics such as peak width, retention time.
# This files were processed using the R code below.
# 
# MZXML File location: smbl.nus.edu.sg/data/metagenome_data/Desktop/Shiv_algae/Algae/MZXML/algae_all
# 
# #The version used was R version 3.0.1 (2013-05-16) -- "Good Sport" and xcms_1.38.0
# 
#### R code for XCMS ##########
############################## Rerunning analysis on the same set. this time storing intermediate steps for CAMERA
#file generated as algae_4setblanks_021213
# library(xcms)
# set1<-xcmsSet(nSlaves=44,method='centWave',ppm=30,peakwidth=c(5,60), prefilter=c(0,0),snthresh=6)
# save(set1,file="set1_algae_4setblanks_021213.rda")
# set2 <- group(set1,bw=5,mzwid=0.015,minsamp=1,minfrac=0) 
# set3 <- retcor(set2,method="obiwarp",plottype="none")
# set4 <- group(set3,bw=5,mzwid=0.015,minsamp=1,minfrac=0)
# save(set4,file="set4_algae_4setblanks_021213.rda")
# set5 <- fillPeaks(set4) 
# save(set5,file="set5_algae_4setblanks_021213.rda")
# peaklist<-peakTable(set5,filebase="algae_4setblanks_021213")

#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 

library(ggplot2)
library(data.table)
library(stringr)
library(vegan)
library(ape)
library(gplots)
library(raster)
library(gridBase)
library(EMT)
library(reshape2)
library(affy)
library(preprocessCore)
library(multtest)
library(compiler)

#############
# User      #
# Specific  #
# Variables #
#############

directory<- "F:/Vinay's Algae Data/Aug2013/Data"
#Contains path for .tsv file
# The "PATH", Remember to use "/" and not "/" in the path
setwd(directory)
getwd()

#################
# Get filenames #
######################## reading in the mass spec data 

mzfilename<-"algae_4setblanks_021213.tsv"
ms_data_total<-read.table(mzfilename,sep="\t",header=T,check.names=FALSE,row.names=1)
str(ms_data_total)

#to avoid memory issues, the file after removing duplicates and outlier replicates is stored as algae_22strains_243sample_021213.tsv
#this file contains 243 samples

################Formatting column names
ori_names_ms_data_total<-names(ms_data_total)
names(ms_data_total)<-ori_names_ms_data_total
names(ms_data_total)<-gsub('\\[', '', names(ms_data_total))
names(ms_data_total)<-gsub('\\]', '', names(ms_data_total))
a<-names(ms_data_total)
b<-gsub('\\(raw)', '', a)
c<-gsub('\\, ','\\-', b)
names(ms_data_total)<-c
rm(a,b,c)

strains_22<-c("D12_001","D12_006","D12_14","D12_051","D12_84","D12_87","D12_94","D12_104","D12_177","D12_184","D12_187","D12_207","D12_245","D12_252","D12_253","D12_254","D12_255","D12_258","D12_268","D12_283","D12_322","D12_325",
              "D4_001","D4_006","D4_14","D4_051","D4_084","D4_087","D4_094","D4_104","D4_177","D4_184","D4_187","D4_207","D4_245","D4_252","D4_253","D4_254","D4_255","D4_258","D4_268","D4_283","D4_322","D4_325")

duplicate_removal<-c("D12_001b1b.r001","D12_001b1b.r002","D12_001b2b.r001","D12_001b2b.r002","D12_001b3b.r001","D12_001b3b.r002",
                     "D4_001b1.r001","D4_001b1.r002","D4_001b2.r001","D4_001b2.r002","D4_001b3.r001","D4_001b3.r002")
ms_data_22<-ms_data_total[, grepl(paste(strains_22,collapse="|"),names(ms_data_total))]
ms_data_22<-ms_data_22[, !names(ms_data_22) %in% duplicate_removal] 
ms_data_22<-ms_data_22[,35:dim(ms_data_22)[2]]
names(ms_data_22)<-gsub('\\.', '\\_', names(ms_data_22))

#Only for full data -22 strains
#ms_data_total_strains<-ms_data_total[, -grep('ACN', names(ms_data_total))] #removing blanks
#ms_data_total_strains<-ms_data_total_strains[, -grep('matri_', names(ms_data_total_strains))] #removing matrix

metadata<-read.table("Algae_indreps_22strains_blanks_matrix_metadata.txt",header=TRUE,row.name=1,sep="\t")
metadata<-as.data.frame(metadata)
metadata_strains<-unique(metadata)
metadata<-t(metadata)
RunDay<-metadata[3,]
GrowthStage<-metadata[1,]

# Biochemical paramters for 11 strains in day 12 group

biochem_para<-read.table("biochemicalparam_day12_replicates.txt",header=TRUE,row.name=1,sep="\t")
biochem_para<-as.data.frame(biochem_para)
biochem_para<-biochem_para[order(biochem_para[,1]),]

# Assigning sample groups

SampleGroups<-sapply(names(ms_data_22), function(x) paste(strsplit(x,"_")[[1]][1:2],collapse = "_"))
SampleGroups<-as.vector(SampleGroups)

SampleGroups<-gsub('b1','',SampleGroups)
SampleGroups<-gsub('b2','',SampleGroups)
SampleGroups<-gsub('b3','',SampleGroups)

#write.table(SampleGroups,"SampleGroups.txt",row.name=FALSE,quote=F)
#SampleGroups1<-read.table("SampleGroups.txt",header=F)

#no_of_info_cols<-1
#data_start_index<-no_of_info_cols+1;
#data_end_index<-dim(ms_data_total)[2];

ms_data<-round(ms_data_22,0) # removing decimal places from intensities. This helps reducing size dataset using the logic that intensities are not so accurate!
#ms_data<-ms_data_total[,1:311]
names(ms_data)
dim(ms_data)
SampleGroup<-SampleGroups

#Measure of variation between replicates of each strain #based on boxplot

groups<-unique(SampleGroup)

groups_with_outliers<-0
outlier_replicates<-0
a<-1

for(i in 1:length(groups))
{
  ms_data_grp1<-ms_data[, grep(groups[i], names(ms_data))]
  ms_data_grp1[ms_data_grp1 < 1e-3]<- 1
  ms_data_grp1<-log(ms_data_grp1)
  boxplot_replicates<-boxplot(ms_data_grp1,plot=FALSE)
  if(ncol(ms_data_grp1)==6)
  {
    cv_group<-c(round(cv(ms_data_grp1[,1]), 2),round(cv(ms_data_grp1[,2]), 2),round(cv(ms_data_grp1[,3]), 2),
                round(cv(ms_data_grp1[,4]), 2),round(cv(ms_data_grp1[,5]), 2),round(cv(ms_data_grp1[,6]), 2))
    boxplot_cv<-boxplot(cv_group,plot=FALSE)
  }
  if(ncol(ms_data_grp1)==4)
  {
    cv_group<-c(round(cv(ms_data_grp1[,1]), 2),round(cv(ms_data_grp1[,2]), 2),
                round(cv(ms_data_grp1[,3]), 2),round(cv(ms_data_grp1[,4]), 2))
    boxplot_cv<-boxplot(cv_group,plot=FALSE)
  }
  if(ncol(ms_data_grp1)==7)
  {
    cv_group<-c(round(cv(ms_data_grp1[,1]), 2),round(cv(ms_data_grp1[,2]), 2),round(cv(ms_data_grp1[,3]), 2),
                round(cv(ms_data_grp1[,4]), 2),round(cv(ms_data_grp1[,5]), 2),round(cv(ms_data_grp1[,6]), 2),
                round(cv(ms_data_grp1[,7]), 2))
    boxplot_cv<-boxplot(cv_group,plot=FALSE)
  }
  
  if(length(boxplot_cv$out)>0)
  {
    groups_with_outliers[a]<-groups[i]
    outlier_replicates[a]<-list(names(ms_data_grp1)[which(cv_group%in%boxplot_cv$out)])
    a<-a+1
  }
}
outlier_replicates<-unlist(outlier_replicates)
# 
# > groups_with_outliers
# [1] "23rd_ACN" "D12_001"  "D12_006"  "D12_14"   "D12_051"  "D12_104"  "D12_207"  "D12_255"  "D12_283"  "D4_094"   "D4_104"   "D4_245"   "D4_252"   "D4_254"   "D4_268"  
# [16] "D4_325" 

# > outlier_replicates
#[1] "D12_001b3_r001" "D12_006b1_r002" "D12_051b3_r002" "D12_14b2_r002"  "D12_177b3_r002"
#[6] "D12_187b1_r001" "D12_207b1_r001" "D12_207b2_r002" "D12_252b1_r001" "D12_252b3_r002"
#[11] "D12_255b1_r001" "D12_258b1_r001" "D4_094_b1_r002" "D4_104_b3_r002" "D4_14b1_r001"  
#[16] "D4_245b1_r001"  "D4_254b1_r001"  "D4_268_b1_r002" "D4_325_b1_r002"

strains_with_biochem<-c("D12_253","D12_051","D12_252","D12_187","D12_255","D12_001","D12_177","D12_258","D12_254","D12_184","D12_87")


###### Data: inclusion and exclusion, removing replicates(samples with large variance) from the dataset 
ms_data<-ms_data[, !names(ms_data) %in% outlier_replicates] #without the outlier replicates
rownames(ms_data)<-paste0(ms_data_total$mz,"@",ms_data_total$rt) #rownames(ms_data_total)
# a<-cbind(ms_data_total$mz,ms_data_total$rt,ms_data)
# write.table(a,"algae_22strains_243sample_021213.tsv",quote=FALSE,sep="\t",row.name=FALSE)

#### ALL FURTHER ANALYSIS TO BE PERFORMED ON STRAINS AFTER REMOIER OUTLIERS ####

SampleGroup<-sapply(names(ms_data), function(x) paste(strsplit(x,"_")[[1]][1:2],collapse = "_"))
SampleGroup<-as.vector(SampleGroup)
SampleGroup<-gsub('b1','',SampleGroup)
SampleGroup<-gsub('b2','',SampleGroup)
SampleGroup<-gsub('b3','',SampleGroup)
SampleGroup#Sample group without the outlier replicates

SampleGroup_day4<-SampleGroup[grep('D4', SampleGroup)] 
SampleGroup_day12<-SampleGroup[grep('D12',SampleGroup)] 

metadata_rem<-metadata[, !colnames(metadata) %in% outlier_replicates]
RunDay<-metadata_rem[3,]
GrowthStage<-metadata_rem[1,]
#RunDay_rem<-as.vector(RunDay_rem)
GrowthStage<-as.vector(GrowthStage)

GrowthStage<-gsub('D4', 1, GrowthStage)
GrowthStage<-gsub('D12', 2, GrowthStage)

#strains alone 24:268

GrowthStage_strains<-GrowthStage[24:268]
GrowthStage_strains<-gsub('D4', 1, GrowthStage_strains)
GrowthStage_strains<-gsub('D12', 2, GrowthStage_strains)

RunDay_strains<-RunDay[names(ms_data)]
RunDay_strains<-as.numeric(RunDay_strains)
names(RunDay_strains)<-names(ms_data)

ms_data_day4<-ms_data[, grep('D4', names(ms_data))] 
ms_data_day12<-ms_data[, grep('D12', names(ms_data))] 
ms_data_day12_bc<-ms_data[, grepl(paste(strains_with_biochem,collapse="|"),names(ms_data))] #11 strains with biochemical param

GrowthStage_day4<-GrowthStage[grep('D4', names(ms_data),perl=TRUE, value=TRUE)]
GrowthStage_day12<-GrowthStage[grep('D12', names(ms_data),perl=TRUE, value=TRUE)] 

RunDay_day4<-RunDay[grep('D4',names(RunDay_strains),perl=TRUE, value=TRUE)]
RunDay_day12<-RunDay[grep('D12', names(RunDay_strains),perl=TRUE, value=TRUE)]

### Extracting features which are non-zero across all strains

ms_data_day4_zero_nonzero<-ms_data_day4
ms_data_day4_zero_nonzero[ms_data_day4_zero_nonzero<1e-3]<-NA
ms_data_day4_nonzero<- ms_data_day4_zero_nonzero[complete.cases(ms_data_day4_zero_nonzero),]
write.table(ms_data_day4_nonzero,"ms_data_day4_nonzero.txt",quote=FALSE,sep="\t",col.names=NA)
#ms_data_day4_nonzero msd4<-read.table("ms_data_day4_nonzero.txt",sep="\t",row.names=1,header=TRUE) #to read in files directly

ms_data_day12_zero_nonzero<-ms_data_day12
ms_data_day12_zero_nonzero[ms_data_day12_zero_nonzero<1e-3]<-NA
ms_data_day12_nonzero<- ms_data_day12_zero_nonzero[complete.cases(ms_data_day12_zero_nonzero),]
write.table(ms_data_day12_nonzero,"ms_data_day12_nonzero.txt",quote=FALSE,sep="\t",col.names=NA)
#ms_data_day12_nonzero<-read.table("ms_data_day12_nonzero.txt",sep="\t",row.names=1,header=TRUE) #to read in files directly
  
ms_data_day12_nonzero_bc<-ms_data_day12_nonzero[, grepl(paste(strains_with_biochem,collapse="|"),names(ms_data_day12_nonzero))] #11 strains with biochemical param
ms_data_day12_nonzero_bc<-ms_data_day12_nonzero_bc[,order(names(ms_data_day12_nonzero_bc))]

########################################
###### Multivariate data analysis ######
########################################
###############################################
## R functions- batch effect removal ##########
###############################################
#function to compute PCA and perform linear models against runday and strain
#PCA using princomp on mz by samp

compute_pca<-function(dataset,preprocess_method) {
  dataset<-as.matrix(log(dataset)) #log transform the data using natural log
  if(preprocess_method=="norm") {
    processed_data<-normalize.quantiles(as.matrix(dataset),copy=TRUE)
  }   else  {
    processed_data<-scale(dataset,center=T,scale=T)
    processed_data<-processed_data-min(processed_data)
  }
  pca_results <- princomp(processed_data,cor=F,scores=T) ### IMP: choose quantile normalized or scaled data
  return(pca_results)
}

#function to compute linear model 
compute_linearModel<-function(results.from.pca,dependent.factor) { #dependent.factor is either RunDay_day4 or 12 (or) SampleGroup_day4 or 12 
  lm_pca_scores<-apply(results.from.pca$loadings,2, function(x) {
    lm_val<-lm(x~ as.factor(dependent.factor)) 
    lm_cor<-summary(lm_val)
    p.val<-anova(lm_val)$'Pr(>F)'[1]
    return(list(lm_cor$r.squared,p.val))
  })
}

#function to extract r2 value from list containg r2 and p.val returned from linear model
compute.r2.pval<-function(linearmodel_list,r2.pval) {
  if(r2.pval=="r2") { #WARNING:code implictly assumes r2 is in the first column and p.val in the second
    return (sapply(linearmodel_list, function(x){as.numeric(x[1])}))
  } else{
    return (sapply(linearmodel_list, function(x){as.numeric(x[2])}))
  }
}

#function to compute svd
compute_svd<-function(dataset,preprocess_method,start_pc_comp,end_pc_comp,recursive_pcs) {
  dataset<-as.matrix(log(dataset)) #log transform the data using natural log
  if(preprocess_method=="norm") {
    processed_data<-normalize.quantiles(as.matrix(dataset),copy=TRUE)
  }   else  {
    processed_data<-scale(dataset,center=T,scale=T)
    processed_data<-processed_data-min(processed_data)
  }
  
  if(recursive_pcs=="recursive") {
    svd_dataset<-svd(processed_data)
    if(end_pc_comp){end_pc<-end_pc_comp} else{end_pc<-ncol(svd_dataset)}
    pc_combi_list<-sapply(start_pc_comp:end_pc,function(x) seq(start_pc_comp:x))
    dataset_rmbatch<-vector("list", length(pc_combi_list)) #list() #
    for(i in 1:length(pc_combi_list))
    {
      svd_dataset$d1<-svd_dataset$d
      svd_dataset$d1[c(unlist(pc_combi_list[i]))]<-0 #Removing the variation caused by runday(id using pca) where the multiple r2 correlation is above 0.5
      dataset_rmbatch1<-svd_dataset$u %*% diag(svd_dataset$d1) %*% t(svd_dataset$v)
      rownames(dataset_rmbatch1)<-rownames(dataset)
      colnames(dataset_rmbatch1)<-names(dataset)
      dataset_rmbatch[[i]]<-dataset_rmbatch1
    }       
  } else{
    dataset_rmbatch<-list()
    svd_dataset<-svd(processed_data)
    svd_dataset$d1<-svd_dataset$d
    svd_dataset$d1[start_pc_comp]<-0 
    end_pc_comp#not used
    dataset_rmbatch1<-svd_dataset$u %*% diag(svd_dataset$d1) %*% t(svd_dataset$v)
    rownames(dataset_rmbatch1)<-rownames(dataset)
    colnames(dataset_rmbatch1)<-names(dataset)
    dataset_rmbatch[[1]]<-dataset_rmbatch1
  }
  return(dataset_rmbatch)
}

# compute permutative f test statistics to difentify signigicant features
compute_perm_ftest<-function(dataset,classlabel) {
  classlabel_factor<-as.numeric(as.factor(classlabel))-1
  if(length(dataset)> 1 ) {
    p.values<-vector("list", length(dataset)) ### Change to p.values<-vector("list", length(dataset)) # make it faster
    sig_metab_dataset<-vector("list", length(dataset)) 
    for(i in 1:length(dataset))
    { data_matrix<-as.matrix(dataset[[i]])
      dataset_sig_features<-mt.maxT(data_matrix,classlabel_factor,test="f",side="abs",fixed.seed.sampling="y",B=100000,nonpara="n")
      p.values[[i]]<-dataset_sig_features
      id.sig_dataset <- sort(dataset_sig_features[dataset_sig_features$adjp < 0.05,c(1)])
      metab_sig<-cbind(data_matrix[id.sig_dataset,],round(dataset_sig_features$adjp[id.sig_dataset],5))
      sig_metab_dataset[[i]]<-metab_sig
    }
  } else { # only a single dataset
    classlabel_factor<-as.numeric(as.factor(classlabel))-1
    p.values<-list()
    sig_metab_dataset<-list()
    data_matrix<-as.matrix(dataset[[1]])
    dataset_sig_features<-mt.maxT(data_matrix,classlabel_factor,test="f",side="abs",fixed.seed.sampling="y",B=100000,nonpara="n")
    p.values[[1]]<-dataset_sig_features
    id.sig_dataset <- sort(dataset_sig_features[dataset_sig_features$adjp < 0.05,c(1)])
    metab_sig<-cbind(data_matrix[id.sig_dataset,],round(dataset_sig_features$adjp[id.sig_dataset],5))
    sig_metab_dataset[[1]]<-metab_sig
  }
  return(list(p.values,sig_metab_dataset))
}

compute_perm_ftest_faster<-cmpfun(compute_perm_ftest)

#computer list of p.val
compute_pval_list<-function(sigfeat_pval) {
  listofpval<-vector("list", length(sigfeat_pval)) #list()
  for(i in 1:length(sigfeat_pval))
  {
    listofpval[i]<-sigfeat_pval[[i]][4]
  }
  sigfeat_pval_pvaldf<-do.call(cbind.data.frame, listofpval)
}

#compute no of sig features
compute_no_features<-function(sig_feat) {
  no_of_sig_features<-rep(NA,length(sig_feat))
  for(i in 1:length(sig_feat))
  {
    no_of_sig_features[i]<-nrow(sig_feat[[i]][1])
  }
  return(no_of_sig_features)
}

###############################################
#day4

fit_day4_mzbysam_princomp<-compute_pca(ms_data_day4_nonzero,"scale")
fit_day4_mzbysam_pc_variation<-fit_day4_mzbysam_princomp$sdev^2/sum(fit_day4_mzbysam_princomp$sdev^2)
lm_pca_scores_strain_day4_nonzero_loadings<-compute_linearModel(fit_day4_mzbysam_princomp,SampleGroup_day4)
lm_pca_scores_strain_day4_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_strain_day4_nonzero_loadings,"r2")
lm_pca_scores_strain_day4_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_strain_day4_nonzero_loadings,"pval")

lm_pca_scores_runday4_nonzero_loadings<-compute_linearModel(fit_day4_mzbysam_princomp,RunDay_day4)
lm_pca_scores_runday4_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_runday4_nonzero_loadings,"r2")
lm_pca_scores_runday4_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_runday4_nonzero_loadings,"pval")

svd_day4_nonzero <- compute_svd(ms_data_day4_nonzero,"scale",1,ncol(ms_data_day4_nonzero),"recursive")
save(svd_day4_nonzero,file='svd_day4_x138_nonzero.rda')
#str(svd_day4_nonzero)
day4_nonzero_sigfeat_r<-compute_perm_ftest_faster(svd_day4_nonzero,RunDay_day4)
day4_nonzero_sigfeat_r_pval<-day4_nonzero_sigfeat_r[[1]]
day4_nonzero_sigfeat_r_matrix<-day4_nonzero_sigfeat_r[[2]]
day4_nonzero_sigfeat_r_pvaldf<-compute_pval_list(day4_nonzero_sigfeat_r_pval)
day4_no_nonzero_sigfeat_r<-compute_no_features(day4_nonzero_sigfeat_r_matrix)
save(day4_nonzero_sigfeat_r,file='day4_x138_nonzero_sigfeat_r.rda')

day4_nonzero_sigfeat_s<-compute_perm_ftest_faster(svd_day4_nonzero,SampleGroup_day4)
day4_nonzero_sigfeat_s_pval<-day4_nonzero_sigfeat_s[[1]]
day4_nonzero_sigfeat_s_matrix<-day4_nonzero_sigfeat_s[[2]]
day4_nonzero_sigfeat_s_pvaldf<-compute_pval_list(day4_nonzero_sigfeat_s_pval)
day4_no_nonzero_sigfeat_s<-compute_no_features(day4_nonzero_sigfeat_s_matrix)
save(day4_nonzero_sigfeat_s,file='day4_x138_nonzero_sigfeat_s.rda')

###################
#### day12
###################

fit_day12_mzbysam_princomp<-compute_pca(ms_data_day12_nonzero,"scale")
fit_day12_mzbysam_pc_variation<-fit_day12_mzbysam_princomp$sdev^2/sum(fit_day12_mzbysam_princomp$sdev^2)
lm_pca_scores_strain_day12_nonzero_loadings<-compute_linearModel(fit_day12_mzbysam_princomp,SampleGroup_day12)
lm_pca_scores_strain_day12_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_strain_day12_nonzero_loadings,"r2")
lm_pca_scores_strain_day12_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_strain_day12_nonzero_loadings,"pval")

lm_pca_scores_runday12_nonzero_loadings<-compute_linearModel(fit_day12_mzbysam_princomp,RunDay_day12)
lm_pca_scores_runday12_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_runday12_nonzero_loadings,"r2")
lm_pca_scores_runday12_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_runday12_nonzero_loadings,"pval")

svd_day12_nonzero <- compute_svd(ms_data_day12_nonzero,"scale",1,ncol(ms_data_day12_nonzero),"recursive")
#svd_day12_nonzero <- compute_svd(ms_data_day12_nonzero,"scale",0,0,"nonrecursive")
save(svd_day12_nonzero,file='svd_day12_x138_nonzero.rda')
#str(svd_day12_nonzero)
day12_nonzero_sigfeat_r<-compute_perm_ftest_faster(svd_day12_nonzero,RunDay_day12)
day12_nonzero_sigfeat_r_pval<-day12_nonzero_sigfeat_r[[1]]
day12_nonzero_sigfeat_r_matrix<-day12_nonzero_sigfeat_r[[2]]
day12_nonzero_sigfeat_r_pvaldf<-compute_pval_list(day12_nonzero_sigfeat_r_pval)
day12_no_nonzero_sigfeat_r<-compute_no_features(day12_nonzero_sigfeat_r_matrix)
save(day12_nonzero_sigfeat_r,file='day12_x138_nonzero_sigfeat_r.rda')

day12_nonzero_sigfeat_s<-compute_perm_ftest_faster(svd_day12_nonzero,SampleGroup_day12)
day12_nonzero_sigfeat_s_pval<-day12_nonzero_sigfeat_s[[1]]
day12_nonzero_sigfeat_s_matrix<-day12_nonzero_sigfeat_s[[2]]
day12_nonzero_sigfeat_s_pvaldf<-compute_pval_list(day12_nonzero_sigfeat_s_pval)
day12_no_nonzero_sigfeat_s<-compute_no_features(day12_nonzero_sigfeat_s_matrix)
save(day12_nonzero_sigfeat_s,file='day12_x138_nonzero_sigfeat_s.rda')


png("pval_factor.png",height=800,width=800)
par(mfrow=c(2,2))
matplot(day4_nonzero_sigfeat_r_pvaldf,type='l',ylab="pval",xlab="day4 Runday")
matplot(day12_nonzero_sigfeat_r_pvaldf,type='l',ylab="pval",xlab="day12 Runday")
matplot(day4_nonzero_sigfeat_s_pvaldf,type='l',ylab="pval",xlab="day4 Strain")
matplot(day12_nonzero_sigfeat_s_pvaldf,type='l',ylab="pval",xlab="day12 Strain")
dev.off()

png("pval_sigfeat.png",height=800,width=800)
par(mfrow=c(2,2))
plot(1:124,day4_no_nonzero_sigfeat_r,type='l',col='blue',xlab="day4 Runday")
plot(1:121,day12_no_nonzero_sigfeat_r,type='l',col='red',xlab="day12 Runday")
plot(1:124,day4_no_nonzero_sigfeat_s,type='l',col='blue',xlab="day4 Strain")
plot(1:121,day12_no_nonzero_sigfeat_s,type='l',col='red',xlab="day12 Strain")
dev.off()



