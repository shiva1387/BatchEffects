###########################################
# Data Analysis in R-Malaysian algae data #
###########################################
# Author(s): Shiv
# Version: 01122013 
# Input: ".tsv" file from XCMS 
# Software: XCMS
# Modified By :Shivshankar Umashankar 
# Modular version of AlgaeDataAnalysis.R
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
# #The version used was R version 3.0.1 (2013-05-16) -- "Good Sport" and xcms_1.36.0
# # Interestingly, the no of features chnaged only for the big matrix of 417 reps. There is no difference between the features detected 
# # between xcms_1.36 and xcms_1.38 for the blank or tt8 matrix 
#### R code for XCMS ##########
# library(xcms)
# 
# set1<-xcmsSet(nSlaves=44,method='centWave',ppm=30,peakwidth=c(5,60), prefilter=c(0,0),snthresh=6) 
# set2 <- group(set1,bw=5,mzwid=0.015,minsamp=1,minfrac=0) 
# set3 <- retcor(set2,method="obiwarp",plottype="none")
# set4 <- group(set3,bw=5,mzwid=0.015,minsamp=1,minfrac=0)
# set5 <- fillPeaks(set4) 
# peaklist<-peakTable(set5,filebase="algae_4setblanks_060913")

##############################
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
library(xcms)
library(stringr)
library(vegan)
library(ape)
library(gplots)
library(raster)
library(gridBase)
library(EMT)
library(RColorBrewer)
library(agricolae)
library(reshape2)
library(affy)
library(preprocessCore)
library(multtest)

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
#################


#mzxmlfiles<-list.files(directory,pattern=".mzXML",full.names=TRUE) # Get filenames with a given pattern and include subdirectories 

####### reading in the mass spec data 

mzfilename<-"Algae_22strains_indreps_blanks_matrix.tsv"
ms_data_total<-read.table(mzfilename,sep="\t",header=T,check.names=FALSE,row.names=1)
str(ms_data_total)

################Formatting column names

names(ms_data_total)<-gsub('\\[', '', names(ms_data_total))
names(ms_data_total)<-gsub('\\]', '', names(ms_data_total))
a<-names(ms_data_total)
b<-gsub('\\(raw)', '', a)
c<-gsub('\\, ','\\-', b)
names(ms_data_total)<-c
rm(a,b,c)

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

SampleGroups<-sapply(names(ms_data_total), function(x) paste(strsplit(x,"_")[[1]][1:2],collapse = "_"))
SampleGroups<-as.vector(SampleGroups)

SampleGroups<-gsub('b1','',SampleGroups)
SampleGroups<-gsub('b2','',SampleGroups)
SampleGroups<-gsub('b3','',SampleGroups)

#write.table(SampleGroups0,"SampleGroups.txt",row.name=FALSE,quote=F)
#SampleGroups1<-read.table("SampleGroups.txt",header=F)

#no_of_info_cols<-1
#data_start_index<-no_of_info_cols+1;
#data_end_index<-dim(ms_data_total)[2];

ms_data<-ms_data_total[,1:311]
names(ms_data)
dim(ms_data)
SampleGroup<-SampleGroups[1:311]


#ms_data<-log(ms_data)

#name_list <- strsplit(SampleGroup, "_")
#GrowthStage<-sapply(name_list , function (x) if(length(x) == 2) x[1] else as.character(NA))
#StrainId<-sapply(name_list , function (x) if(length(x) == 2) x[2] else as.character(NA))
#StrainName<-as.factor(names(ms_data))
#StrainName0<-names(ms_data)

#Measure of variation between replicates of each strain #based on boxplot

groups<-unique(SampleGroup)

groups_with_outliers<-0
outlier_replicates<-0
a<-1

for(i in 1:length(groups))
{
  ms_data_grp1<-ms_data[, grep(groups[i], names(ms_data))]
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
# [1] "23rd_ACN_Blank01" "D12_001b3_r001"   "D12_006b1_r002"   "D12_14b2_r002"    "D12_14b3_r001"    "D12_051b3_r002"   "D12_104b1_r001"   "D12_207b2_r002"   "D12_255b1_r001"  
# [10] "D12_283b1_r002"   "D4_094_b1_r002"   "D4_104_b3_r002"   "D4_245b1_r001"    "D4_252b2_r001"    "D4_254b1_r001"    "D4_254b1_r002"    "D4_268_b1_r002"   "D4_325_b1_r002"  

strains_with_biochem<-c("D12_253","D12_051","D12_252","D12_187","D12_255","D12_001","D12_177","D12_258","D12_254","D12_184","D12_87")

###### Data: inclusion and exclusion, removing replicates(samples with large variance) from the dataset 

#rem <- c('D4_094_b1_r002','D4_104_b3_r002','D4_245b1_r001','D4_254b1_r001','D4_325_b1_r002')

ms_data<-ms_data[, !names(ms_data) %in% outlier_replicates] #without the outlier replicates
rownames(ms_data)<-rownames(ms_data_total)
#write.table(ms_data,"ms_data_wo_outliers.tsv",sep='\t',quote=FALSE,col.names=NA)

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

RunDay_strains<-as.numeric(RunDay[24:268])
names(RunDay_strains)<-colnames(ms_data)[24:268]

ms_data_day4<-ms_data[, grep('D4', names(ms_data))] 
ms_data_day12<-ms_data[, grep('D12', names(ms_data))] 
ms_data_day12_bc<-ms_data[, grepl(paste(strains_with_biochem,collapse="|"),names(ms_data))] #11 strains with biochemical param

GrowthStage_day4<-GrowthStage[grep('D4', names(GrowthStage),perl=TRUE, value=TRUE)]
GrowthStage_day12<-GrowthStage[grep('D12', names(GrowthStage),perl=TRUE, value=TRUE)] 

RunDay_day4<-RunDay[grep('D4',names(RunDay_strains),perl=TRUE, value=TRUE)]
RunDay_day12<-RunDay[grep('D12', names(RunDay_strains),perl=TRUE, value=TRUE)]

zeroes_in_column_rem<-z1[, !colnames(z1) %in% outlier_replicates]
zeroes_in_column_day4<-z1[grep('D4', colnames(z1))] 
zeroes_in_column_day12<-z1[grep('D12', colnames(z1))]

mz_grp_zero_day4<-mz_grp_zero[, grep('D4', names(mz_grp_zero))] 
mz_grp_zero_day12<-mz_grp_zero[, grep('D12', names(mz_grp_zero))]


# converting to data.table format

ms_data_table<-data.table(t(log(ms_data))) #log transformed
ms_data_table<-cbind(ms_data_table,as.factor(SampleGroup))
lapply(ms_data_table,class) #checking col classes

ms_data_day12_zero[is.na(ms_data_day12_zero)]<-0
ms_data_day12_zero_table<-data.table(t(log(ms_data_day12_zero)))
ms_data_day12_zero_table<-cbind(ms_data_day12_zero_table,as.factor(SampleGroup_day12))
lapply(ms_data_day12_zero_table,class) #checking col classes

ms_data_day4_zero[is.na(ms_data_day4_zero)]<-0
ms_data_day4_zero_table<-data.table(t(log(ms_data_day4_zero)))
ms_data_day4_zero_table<-cbind(ms_data_day4_zero_table,as.factor(SampleGroup_day4))
lapply(ms_data_day4_zero_table,class) #checking col classes

#######################################
###### Exploratory data analysis ######
#######################################

##########################
# Histogram of mz and rt #
##########################

mz_rt <- strsplit(rownames(ms_data_total), "\\@")
mz<-sapply(mz_rt , function (x) if(length(x) == 2) x[1] else as.character(NA))
rt<-sapply(mz_rt , function (x) if(length(x) == 2) x[2] else as.character(NA))

png("algae_417reps_retcor_MPP_sorted.png")
plot(rt,mz,pch=1,main="algae_417reps_retcor_MPP_sorted",ylab="m/z",xlab="rt")
dev.off()

### Measure of variation ########

png('boxplot.png',width=8000)
boxplot(log(ms_data),range=0,xlab="Samples",ylab="Log2 feature counts",las=1,cex.axis=0.7)
dev.off()

### Visual representation of runday Vs strains ####

runday_strains<-read.table("algae_metainfor_runday_strains.txt",header=TRUE,sep="\t")
acast(tmp, x~y, value.var="z")
runday_strains_matrix<-acast(runday_strains, SampleName~RunDay, value.var="Value")
runday_strains_matrix[is.na(runday_strains_matrix)]<-0
pdf("RunDayVsStrains.pdf",width=16,height=16)
pheatmap(runday_strains_matrix,cluster_rows=FALSE,cluster_cols=FALSE,scale="none",breaks=c(12,16,18,22,25),col=brewer.pal(4,"Accent"),fontsize=12)
dev.off()

########################################################################
###### Systems level analysis of non-detected features in ms data ######
########################################################################

zero_count_column<-function(x){ d=as.vector(x); sum(d < 1+1e-3) } ### Counting zeros in each column # here zero is given as 1
total_feature_strain<-function(x) length(x)                       ### Counting number of features for each m/z feature in each group
zero_count<-function(x) sum(x < 1+1e-3)                             ### Counting number of zeroes for each m/z feature in each group

feature_present_count<-function(x) {                              ### Counting number of features observed for each m/z feature in each group
  no_of_zero<-sum(x < 1+1e-3) 
  no_feature<-length(x)-no_of_zero
  return(no_feature)
}

zero_percent_strain<-function(x) {                              ### Counting number of features observed for each m/z feature in each group
  no_of_zero<-sum(x < 1+1e-3) 
  zero_percent<-round((no_of_zero/length(x)),2)
  return(zero_percent)
}

all_zero<-function(x){ d=as.vector(x); sum(d == 1) } 

#column numbers
#blanks: 1:23
#day12: 24:144
#day4: 145:268
#matrix: 269:293

zeroes_in_column<-lapply(ms_data,zero_count_column) 
zeroes_in_column<-as.numeric(zeroes_in_column)
z1<-t(zeroes_in_column)
names(z1)<-colnames(ms_data)


zero_bins<-c(0,500,1000,1500,2000,2500,3000,3500,8000,14000)


png("zeros_per_column.png",res=150,width=1200,height=1200)
par(mfrow=c(2,2))
hist(z1[1:23],breaks=zero_bins, xaxt='n',main="Blanks",xlab="No of zeros",ylab="No of samples:n-23",freq=TRUE)
axis(side=1, at=zero_bins, labels=zero_bins)
hist(z1[145:268],breaks=zero_bins, xaxt='n',main="Day4",xlab="No of zeros",ylab="No of samples:n-124",freq=TRUE)
axis(side=1, at=zero_bins, labels=zero_bins)
hist(z1[269:293],breaks=zero_bins, xaxt='n',main="Matrix",xlab="No of zeros",ylab="No of samples:n-25",freq=TRUE)
axis(side=1, at=zero_bins, labels=zero_bins)
hist(z1[24:144],breaks=zero_bins, xaxt='n',main="Day12",xlab="No of zeros",ylab="No of samples:n-121",freq=TRUE)
axis(side=1, at=zero_bins, labels=zero_bins)
dev.off()


mz_feature_total2<-t(ms_data_table[,lapply(.SD,total_feature_strain),by=V2])
mz_feature_total<-as.data.frame(mz_feature_total2[2:nrow(mz_feature_total2),])
colnames(mz_feature_total)<-unique(SampleGroup)

write.table(mz_feature_total,"mz_features_total.txt",quote=FALSE,sep="\t",col.names=NA)
mz_grp_total<-read.table("mz_features_total.txt",sep='\t',header=TRUE,row.name=1)

mz_feature_observed2<-t(ms_data_table[,lapply(.SD,feature_present_count),by=V2])
mz_feature_observed<-as.data.frame(mz_feature_observed2[2:nrow(mz_feature_observed2),])
colnames(mz_feature_observed)<-unique(SampleGroup)

write.table(mz_feature_observed,"mz_features_observed.txt",quote=FALSE,sep="\t",col.names=NA)
mz_grp_observed<-read.table("mz_features_observed.txt",sep='\t',header=TRUE,row.name=1)


zero_mz_feature2<-t(ms_data_table[,lapply(.SD,zero_count),by=V2])
zero_mz_feature<-as.data.frame(zero_mz_feature2[2:nrow(zero_mz_feature2),])
colnames(zero_mz_feature)<-unique(SampleGroup)
write.table(zero_mz_feature,"mz_features_zeroes.txt",quote=FALSE,sep="\t",col.names=NA)
mz_grp_zero<-read.table("mz_features_zeroes.txt",sep='\t',header=TRUE,row.name=1)


ms_data_d4_zerocount2<-t(ms_data_day4_zero_table[,lapply(.SD,zero_count),by=V2])
ms_data_d4_zerocount<-as.data.frame(ms_data_d4_zerocount2[2:nrow(ms_data_d4_zerocount2),])
colnames(ms_data_d4_zerocount)<-unique(SampleGroup_day4)
write.table(ms_data_d4_zerocount,"ms_data_d4_zerocount.txt",quote=FALSE,sep="\t",col.names=NA)
ms_data_d4_zerocount<-read.table("ms_data_d4_zerocount.txt",sep='\t',header=TRUE,row.name=1)

ms_data_d12_zerocount2<-t(ms_data_day12_zero_table[,lapply(.SD,zero_count),by=V2])
ms_data_d12_zerocount<-as.data.frame(ms_data_d12_zerocount2[2:nrow(ms_data_d12_zerocount2),])
colnames(ms_data_d12_zerocount)<-unique(SampleGroup_day12)
write.table(ms_data_d12_zerocount,"ms_data_d12_zerocount.txt",quote=FALSE,sep="\t",col.names=NA)
ms_data_d12_zerocount<-read.table("ms_data_d12_zerocount.txt",sep='\t',header=TRUE,row.name=1)

zero_percent_mz_feature2<-t(ms_data_table[,lapply(.SD,zero_percent_strain),by=V2])
zero_percent_mz_feature<-as.data.frame(zero_percent_mz_feature2[2:nrow(zero_percent_mz_feature2),])
colnames(zero_percent_mz_feature)<-unique(SampleGroup)
write.table(zero_percent_mz_feature,"mz_features_zero_percent.txt",quote=FALSE,sep="\t",col.names=NA)
mz_grp_zero_percent<-read.table("mz_features_zero_percent.txt",sep='\t',header=TRUE,row.name=1)

no_features_all_zero<-lapply(mz_grp_zero_percent,all_zero) 
names(no_features_all_zero)<-unique(SampleGroup)

png("nof-features-allblanks.png",res=150,width=1200,height=1200)
par(mfrow=c(2,2))
plot(1:4,no_features_all_zero[1:4], xaxt = "n",xlab='Blanks', ylab='no of features-all reps zero')
axis(1, at=1:4, labels=names(no_features_all_zero)[1:4])
plot(1:22,no_features_all_zero[27:48], xaxt = "n",xlab='Day 4 Strains', ylab='no of features-all reps zero')
axis(1, at=1:22, labels=names(no_features_all_zero)[5:26],cex.axis=0.6)
plot(1:4,no_features_all_zero[49:52], xaxt = "n",xlab='Matrix', ylab='no of features-all reps zero')
axis(1, at=1:4, labels=names(no_features_all_zero)[49:52])
plot(1:22,no_features_all_zero[5:26], xaxt = "n",xlab='Day 12 Strains', ylab='no of features-all reps zero')
axis(1, at=1:22, labels=names(no_features_all_zero)[5:26],cex.axis=0.6)
dev.off()

#day12 is columns 5:26 and day 4 is columns 27:48

### Extracting features which are non-zero across all strains
###

ms_data_day4_zero_nonzero<-ms_data_day4
ms_data_day4_zero_nonzero[ms_data_day4_zero_nonzero<1+1e-3]<-NA
ms_data_day4_nonzero<- ms_data_day4_zero_nonzero[complete.cases(ms_data_day4_zero_nonzero),]
ms_data_day4_zero<- ms_data_day4_zero_nonzero[!complete.cases(ms_data_day4_zero_nonzero),]
ms_data_day4_nonzero[is.na(ms_data_day4_nonzero)]<-0
write.table(ms_data_day4_nonzero,"ms_data_day4_nonzero.txt",quote=FALSE,sep="\t",col.names=NA)

ms_data_day12_zero_nonzero<-ms_data_day12
ms_data_day12_zero_nonzero[ms_data_day12_zero_nonzero<1+1e-3]<-NA
ms_data_day12_nonzero<- ms_data_day12_zero_nonzero[complete.cases(ms_data_day12_zero_nonzero),]
ms_data_day12_zero<- ms_data_day12_zero_nonzero[!complete.cases(ms_data_day12_zero_nonzero),]
ms_data_day12_nonzero[is.na(ms_data_day12_nonzero)]<-0
write.table(ms_data_day12_nonzero,"ms_data_day12_nonzero.txt",quote=FALSE,sep="\t",col.names=NA)

ms_data_day12_nonzero_bc<-ms_data_day12_nonzero[, grepl(paste(strains_with_biochem,collapse="|"),names(ms_data_day12_nonzero))] #11 strains with biochemical param
ms_data_day12_nonzero_bc<-ms_data_day12_nonzero_bc[,order(names(ms_data_day12_nonzero_bc))]

mz_grp_zero<-data.frame(cbind(mz,rt,mz_grp_zero))

for (i in seq(3,ncol(mz_grp_zero),2))
{
  mz<-as.vector(mz_grp_zero$mz)
  rt<-as.vector(mz_grp_zero$rt)
  ## now plot scatter plot
  png(paste("noofZeros/", i, ".png", sep = ""),height=800,width=1200)
  par(mfrow=c(2,2))
  plot(mz,mz_grp_zero[,i],main=names(mz_grp_zero)[i],xlab="mz",ylab="No of zeros")
  hist(mz_grp_zero[,i],main=names(mz_grp_zero)[i],xlab="no of zeroes",ylab="Features", breaks=6)
  plot(mz,mz_grp_zero[,i+1],main=names(mz_grp_zero)[i+1],xlab="mz",ylab="No of zeros")
  hist(mz_grp_zero[,i+1],main=names(mz_grp_zero)[i+1],xlab="no of zeroes",ylab="Features", breaks=6)
  dev.off()
}
graphics.off()

# plotting %of zeroes for each feature

mz_grp_zero_percent<-data.frame(cbind(mz,rt,mz_grp_zero_percent))

for (i in seq(3,ncol(mz_grp_zero_percent),4))
{
  mz<-as.vector(mz_grp_zero_percent$mz)
  rt<-as.vector(mz_grp_zero_percent$rt)
  ## now plot scatter plot
  png(paste("percentofzeroes/", i, ".png", sep = ""),height=800,width=1200)
  par(mfrow=c(2,2))
  hist(mz_grp_zero_percent[,i],main=names(mz_grp_zero)[i],xlab="% of replicates which are zeros",ylab="no of Features", breaks=6)
  hist(mz_grp_zero_percent[,i+1],main=names(mz_grp_zero)[i+1],xlab="% of replicates which are zeros",ylab="no of Features", breaks=6)
  hist(mz_grp_zero_percent[,i],main=names(mz_grp_zero)[i+2],xlab="% of replicates which are zeros",ylab="no of Features", breaks=6)
  hist(mz_grp_zero_percent[,i],main=names(mz_grp_zero)[i+3],xlab="% of replicates which are zeros",ylab="no of Features", breaks=6)
  dev.off()
}
graphics.off()

#Mean-Vs-Zeros

for (i in seq(3,ncol(mz_grp_mean),2))
{
  mz<-as.vector(mz_grp_mean$mz)
  rt<-as.vector(mz_grp_mean$rt)
  ## now plot scatter plot
  png(paste("Mean-Vs-Zeros/", i, ".png", sep = ""),height=800,width=1200)
  par(mfrow=c(2,2))
  
  cor1<-cor(mz_grp_zero[,i],mz_grp_mean[,i], use="all.obs", method="pearson")
  plot(mz_grp_zero[,i],mz_grp_mean[,i],main=names(mz_grp_zero)[i],xlab="No of zeroes",ylab="mean")
  text(3, 10, paste("R2 =", round(cor1, 3)),cex=1)
  hist(mz_grp_mean[,i],main=names(mz_grp_zero)[i],xlab="mean",ylab="No of Features")
  
  cor2<-cor(mz_grp_zero[,i+1],mz_grp_mean[,i+1], use="all.obs", method="pearson")  
  plot(mz_grp_zero[,i+1],mz_grp_mean[,i+1],main=names(mz_grp_zero)[i+1],xlab="No of zeroes",ylab="mean")
  text(3, 10, paste("R2 =", round(cor2, 3)),cex=1)
  hist(mz_grp_mean[,i+1],main=names(mz_grp_zero)[i+1],xlab="mean",ylab="No of Features")
  
  dev.off()
}
graphics.off()


############################################################################
###### End-Systems level analysis of non-detected features in ms data ######
############################################################################

## Generating boxplot and cv/mean for each sample group

groups<-unique(SampleGroup)
no<-1

no_of_features<-as.integer(nrow(ms_data_total))

for(i in 1:length(groups))
{
  png(paste(groups[i],".png"),width=2000,height=800)
  ms_data_grp1<-ms_data[, grep(groups[i], names(ms_data))]
  par(mar=c(10,2,4,2),mfrow=c(1,2))
  ms_data_grp1<-log(ms_data_grp1)
  boxplot(ms_data_grp1,range=0,ylab="Log2 feature counts", main=paste0(groups[i],"\n","For each col:mean,cv,no of zeroes"),xaxt="n",las=1)
  if(ncol(ms_data_grp1)==6)
  {
    axis(side=1,cex.axis=1,at=seq(1,6,1),as.character(c(paste0(names(ms_data_grp1[1]),"\n",round(mean(ms_data_grp1[,1]), 2),",",round(cv(ms_data_grp1[,1]), 2),"\n",zeroes_in_column[no]),
                                                        paste0(names(ms_data_grp1[2]),"\n",round(mean(ms_data_grp1[,2]), 2),",",round(cv(ms_data_grp1[,2]), 2),"\n",zeroes_in_column[no+1]),
                                                        paste0(names(ms_data_grp1[3]),"\n",round(mean(ms_data_grp1[,3]), 2),",",round(cv(ms_data_grp1[,3]), 2),"\n",zeroes_in_column[no+2]),
                                                        paste0(names(ms_data_grp1[4]),"\n",round(mean(ms_data_grp1[,4]), 2),",",round(cv(ms_data_grp1[,4]), 2),"\n",zeroes_in_column[no+3]),
                                                        paste0(names(ms_data_grp1[5]),"\n",round(mean(ms_data_grp1[,5]), 2),",",round(cv(ms_data_grp1[,5]), 2),"\n",zeroes_in_column[no+4]),
                                                        paste0(names(ms_data_grp1[6]),"\n",round(mean(ms_data_grp1[,6]), 2),",",round(cv(ms_data_grp1[,6]), 2),"\n",zeroes_in_column[no+5])   
    )))
    
    
    
    ### Plotting colored bar charts explaining % of values in each col
    
    c = c(((no_of_features-zeroes_in_column[no])/no_of_features)*100,(zeroes_in_column[no]/no_of_features)*100,((no_of_features-zeroes_in_column[no+1])/no_of_features)*100,(zeroes_in_column[no+1]/no_of_features)*100,
          ((no_of_features-zeroes_in_column[no+2])/no_of_features)*100,(zeroes_in_column[no+2]/no_of_features)*100,((no_of_features-zeroes_in_column[no+3])/no_of_features)*100,(zeroes_in_column[no+3]/no_of_features)*100,
          ((no_of_features-zeroes_in_column[no+4])/no_of_features)*100,(zeroes_in_column[no+4]/no_of_features)*100,((no_of_features-zeroes_in_column[no+5])/no_of_features)*100,(zeroes_in_column[no+5]/no_of_features)*100)
    a = c("Present","Zero","Present","Zero","Present","Zero","Present","Zero","Present","Zero","Present","Zero")
    b = c(names(ms_data_grp1[1]),names(ms_data_grp1[1]),names(ms_data_grp1[2]),names(ms_data_grp1[2]),names(ms_data_grp1[3]),names(ms_data_grp1[3]),
          names(ms_data_grp1[4]),names(ms_data_grp1[4]),names(ms_data_grp1[5]),names(ms_data_grp1[5]),names(ms_data_grp1[6]),names(ms_data_grp1[6]))
    
    dat = data.frame(Group=a, Member=b, Percentage=c)
    
    ## the last one is the current plot
    plot.new()              ## suggested by @Josh
    vps <- baseViewports()
    pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
    vp1 <-plotViewport(c(1,1,0,1)) ## create new vp with margins, you play with this values 
    require(ggplot2)
    
    p <-  ggplot(dat, aes(x = factor(Member),  y = Percentage, fill=Group))+ xlab("Replicates") + geom_bar(stat = "identity")+ labs(title= "Features Pres Vs Abs\n")+
      ## some setting in the title to get something near to the other plots
      theme(plot.title = element_text(size = rel(2),face ='bold'),axis.text.x = element_text(size=15))+scale_fill_manual(values = c('red','green'))
    print(p,vp = vp1)        ## suggested by @bpatiste
    
    no<-no+6
    
  }
  else if(ncol(ms_data_grp1)==4)
  {
    axis(side=1,at=seq(1,4,1),as.character(c(paste0(names(ms_data_grp1[1]),"\n",round(mean(ms_data_grp1[,1]), 2),",",round(cv(ms_data_grp1[,1]), 2),"\n",zeroes_in_column[no]),
                                             paste0(names(ms_data_grp1[2]),"\n",round(mean(ms_data_grp1[,2]), 2),",",round(cv(ms_data_grp1[,2]), 2),"\n",zeroes_in_column[no+1]),
                                             paste0(names(ms_data_grp1[3]),"\n",round(mean(ms_data_grp1[,3]), 2),",",round(cv(ms_data_grp1[,3]), 2),"\n",zeroes_in_column[no+2]),
                                             paste0(names(ms_data_grp1[4]),"\n",round(mean(ms_data_grp1[,4]), 2),",",round(cv(ms_data_grp1[,4]), 2),"\n",zeroes_in_column[no+3])   
    )))
    
    ### Plotting colored bar charts explaining % of values in each col
    
    c = c(((no_of_features-zeroes_in_column[no])/no_of_features)*100,(zeroes_in_column[no]/no_of_features)*100,((no_of_features-zeroes_in_column[no+1])/no_of_features)*100,(zeroes_in_column[no+1]/no_of_features)*100,
          ((no_of_features-zeroes_in_column[no+2])/no_of_features)*100,(zeroes_in_column[no+2]/no_of_features)*100,((no_of_features-zeroes_in_column[no+3])/no_of_features)*100,(zeroes_in_column[no+3]/no_of_features)*100)
    a = c("Present","Zero","Present","Zero","Present","Zero","Present","Zero")
    b = c(names(ms_data_grp1[1]),names(ms_data_grp1[1]),names(ms_data_grp1[2]),names(ms_data_grp1[2]),
          names(ms_data_grp1[3]),names(ms_data_grp1[3]),names(ms_data_grp1[4]),names(ms_data_grp1[4]))
    
    dat = data.frame(Group=a, Member=b, Percentage=c)
    
    ## the last one is the current plot
    plot.new()              ## suggested by @Josh
    vps <- baseViewports()
    pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
    vp1 <-plotViewport(c(1,1,0,1)) ## create new vp with margins, you play with this values 
    require(ggplot2)
    
    p <-  ggplot(dat, aes(x = factor(Member),  y = Percentage, fill=Group))+ xlab("Replicates") + geom_bar(stat = "identity")+ labs(title= "Features Pres Vs Abs\n")+
      ## some setting in the title to get something near to the other plots
      theme(plot.title = element_text(size = rel(1.4),face ='bold'),axis.text.x = element_text(size=15))+scale_fill_manual(values = c('red','green'))
    print(p,vp = vp1)        ## suggested by @bpatiste
    
    no<-no+4
    
  }
  else if(ncol(ms_data_grp1)==7)
  {
    axis(side=1,at=seq(1,7,1),as.character(c(paste0(names(ms_data_grp1[1]),"\n",round(mean(ms_data_grp1[,1]), 2),",",round(cv(ms_data_grp1[,1]), 2),"\n",zeroes_in_column[no]),
                                             paste0(names(ms_data_grp1[2]),"\n",round(mean(ms_data_grp1[,2]), 2),",",round(cv(ms_data_grp1[,2]), 2),"\n",zeroes_in_column[no+1]),
                                             paste0(names(ms_data_grp1[3]),"\n",round(mean(ms_data_grp1[,3]), 2),",",round(cv(ms_data_grp1[,3]), 2),"\n",zeroes_in_column[no+2]),
                                             paste0(names(ms_data_grp1[4]),"\n",round(mean(ms_data_grp1[,4]), 2),",",round(cv(ms_data_grp1[,4]), 2),"\n",zeroes_in_column[no+3]),
                                             paste0(names(ms_data_grp1[5]),"\n",round(mean(ms_data_grp1[,5]), 2),",",round(cv(ms_data_grp1[,5]), 2),"\n",zeroes_in_column[no+4]),
                                             paste0(names(ms_data_grp1[6]),"\n",round(mean(ms_data_grp1[,6]), 2),",",round(cv(ms_data_grp1[,6]), 2),"\n",zeroes_in_column[no+5]),
                                             paste0(names(ms_data_grp1[7]),"\n",round(mean(ms_data_grp1[,7]), 2),",",round(cv(ms_data_grp1[,7]), 2),"\n",zeroes_in_column[no+6])
    )))
    
    ### Plotting colored bar charts explaining % of values in each col
    
    c = c(((no_of_features-zeroes_in_column[no])/no_of_features)*100,(zeroes_in_column[no]/no_of_features)*100,((no_of_features-zeroes_in_column[no+1])/no_of_features)*100,(zeroes_in_column[no+1]/no_of_features)*100,
          ((no_of_features-zeroes_in_column[no+2])/no_of_features)*100,(zeroes_in_column[no+2]/no_of_features)*100,((no_of_features-zeroes_in_column[no+3])/no_of_features)*100,(zeroes_in_column[no+3]/no_of_features)*100,
          ((no_of_features-zeroes_in_column[no+4])/no_of_features)*100,(zeroes_in_column[no+4]/no_of_features)*100,((no_of_features-zeroes_in_column[no+5])/no_of_features)*100,(zeroes_in_column[no+5]/no_of_features)*100,
          ((no_of_features-zeroes_in_column[no+5])/no_of_features)*100,(zeroes_in_column[no+5]/no_of_features)*100)
    a = c("Present","Zero","Present","Zero","Present","Zero","Present","Zero","Present","Zero","Present","Zero","Present","Zero")
    b = c(names(ms_data_grp1[1]),names(ms_data_grp1[1]),names(ms_data_grp1[2]),names(ms_data_grp1[2]),names(ms_data_grp1[3]),names(ms_data_grp1[3]),
          names(ms_data_grp1[4]),names(ms_data_grp1[4]),names(ms_data_grp1[5]),names(ms_data_grp1[5]),names(ms_data_grp1[6]),names(ms_data_grp1[6]),
          names(ms_data_grp1[7]),names(ms_data_grp1[7]))
    
    dat = data.frame(Group=a, Member=b, Percentage=c)
    
    ## the last one is the current plot
    plot.new()              ## suggested by @Josh
    vps <- baseViewports()
    pushViewport(vps$figure) ##   I am in the space of the autocorrelation plot
    vp1 <-plotViewport(c(1,1,0,1)) ## create new vp with margins, you play with this values 
    require(ggplot2)
    
    p <-  ggplot(dat, aes(x = factor(Member),  y = Percentage, fill=Group))+ xlab("Replicates") + geom_bar(stat = "identity")+ labs(title= "Features Pres Vs Abs\n")+
      ## some setting in the title to get something near to the other plots
      theme(plot.title = element_text(size = rel(1.4),face ='bold'),axis.text.x = element_text(size=15))+scale_fill_manual(values = c('red','green'))
    print(p,vp = vp1)        ## suggested by @bpatiste
    
    no<-no+7
  }
  dev.off()
  graphics.off()
}
rm(no)
#plot(t(log1p(ms_data))~SampleGroup,range=0,xlab="Samples",ylab="Log2 feature counts",las=1,col=SampleGroup)

# plotting coefficient of variations among samples(cols)

#Coefficient of variation
png("cv_samples.png",width=4000)
plot(1:ncol(ms_data),apply(ms_data,2,cv)) # For each column 
dev.off()

# cv for each feature in each sample group

cv1 <- function(x) (sd(x)/mean(x)) * 100
cv_mz_feature2<-t(ms_data_tst[,lapply(.SD,cv1),by=V2])
#lapply(cv_mz_feature,class)
cv_mz_feature<-as.data.frame(cv_mz_feature2[2:nrow(cv_mz_feature2),])
colnames(cv_mz_feature)<-unique(SampleGroup)
write.table(cv_mz_feature,"mz_features_cv.txt",quote=FALSE,sep="\t")
mz_grp_cv<-read.table("mz_features_cv.txt",sep='\t',header=TRUE,row.name=1)
mz_grp_cv[is.na(mz_grp_cv)]<-0

# mean for each feature in each sample group

mean_mz_feature2<-t(ms_data_tst[,lapply(.SD,mean),by=V2])
mean_mz_feature<-as.data.frame(mean_mz_feature2[2:nrow(mean_mz_feature2),])
colnames(mean_mz_feature)<-unique(SampleGroup)
write.table(mean_mz_feature,"mz_features_mean.txt",quote=FALSE,sep="\t")
mz_grp_mean<-read.table("mz_features_mean.txt",sep='\t',header=TRUE,row.name=1)
mz_grp_mean[is.na(mz_grp_mean)]<-0

# mean for each feature(excluding zeroes) in each sample group

#set.seed(1234) 
#testvec <- sample(0:10, 100, replace=TRUE) 
mean_wo_zero <- function(x) {
  no_of_zero<-sum(x < 1e-3)
  mean1<-sum(x)/(length(x)-no_of_zero)
  return(mean1)}

#mean(testvec)
#[1] 4.31
#mean(testvec[testvec != 0]) 
#[1] 4.842697
# mean_wo_zero(testvec)
# [1] 4.842697

mean_non_zero_mz_feature2<-t(ms_data_tst[,lapply(.SD,mean_wo_zero),by=V2])
mean_non_zero_mz_feature<-as.data.frame(mean_non_zero_mz_feature2[2:nrow(mean_non_zero_mz_feature2),])
colnames(mean_non_zero_mz_feature)<-unique(SampleGroup)
write.table(mean_non_zero_mz_feature,"mean_non_zero_mz_feature.txt",quote=FALSE,sep="\t")
mz_grp_mean_non_zero<-read.table("mean_non_zero_mz_feature.txt",sep='\t',header=TRUE,row.name=1)
mz_grp_mean_non_zero[is.na(mz_grp_mean_non_zero)]<-0

####

mz_grp_mean<-data.frame(cbind(mz,rt,mz_grp_mean))
mz_grp_cv<-data.frame(cbind(mz,rt,mz_grp_cv))

for (i in seq(3,ncol(mz_grp_mean),4))
{
  mz<-as.vector(mz_grp_mean$mz)
  rt<-as.vector(mz_grp_mean$rt)
  rt<-as.numeric(rt)
  print (i)
  y<-mz_grp_mean[,i]
  rt_mean<-cbind(as.character(rt),as.numeric(y))
  rt_mean<-data.frame(rt,mz_grp_mean[,i],mz_grp_mean[,i+1],mz_grp_mean[,i+2],mz_grp_mean[,i+3])
  rt_mean<-rt_mean[order(rt_mean[,1]),]
  
  ## now plot scatter plot
  png(paste("mz/", i, "_mz.png", sep = ""),height=800,width=1200)
  par(mfrow=c(2,2))
  plot(mz,mz_grp_mean[,i],main=names(mz_grp_mean)[i],xlab="mz",ylab="mean of log(counts)")
  plot(mz,mz_grp_mean[,i+1],main=names(mz_grp_mean)[i+1],xlab="mz",ylab="mean of log(counts)")
  plot(mz,mz_grp_mean[,i+2],main=names(mz_grp_mean)[i+2],xlab="mz",ylab="mean of log(counts)")
  plot(mz,mz_grp_mean[,i+3],main=names(mz_grp_mean)[i+3],xlab="mz",ylab="mean of log(counts)")
  dev.off()
  
  ## now plot histogram
  # hist(mz_grp_mean[,i],main=names(mz_grp_mean)[i],ylab="number of features",xlab="mean of log(counts)")
  # hist(mz_grp_mean[,i+1],main=names(mz_grp_mean)[i+1],ylab="number of features",xlab="mean of log(counts)")
  # hist(mz_grp_mean[,i+2],main=names(mz_grp_mean)[i+2],ylab="number of features",xlab="mean of log(counts)")
  # hist(mz_grp_mean[,i+3],main=names(mz_grp_mean)[i+3],ylab="number of features",xlab="mean of log(counts)")
  
  png(paste("rt/", i, "_rt.png", sep = ""),height=800,width=1200)
  par(mfrow=c(2,2))
  plot(rt_mean[,1],rt_mean[,2],main=names(mz_grp_mean)[i],xlab="rt",ylab="mean of log(counts)")
  plot(rt_mean[,1],rt_mean[,3],main=names(mz_grp_mean)[i+1],xlab="rt",ylab="mean of log(counts)")
  plot(rt_mean[,1],rt_mean[,4],main=names(mz_grp_mean)[i+2],xlab="rt",ylab="mean of log(counts)")
  plot(rt_mean[,1],rt_mean[,5],main=names(mz_grp_mean)[i+3],xlab="rt",ylab="mean of log(counts)")
  dev.off()
}
graphics.off()
# rm(rt_cv)
dev.off()


#### Internal standard

#4473

png("internalStandard_reps.png",height=800,width=800)
par(xpd=TRUE,mar=c(10,4,4,2))
plot(1:311,log(ms_data[4473,]),xaxt='n',las=2,xlab="replicates",ylab="log(counts)",pch=1,main=paste("Internal standard replicates-fmoc"," ",row.names(ms_data)[4473]),
     col=ifelse(metadata$RunDay==15,'blue',
                ifelse(metadata$RunDay==17,'green',
                       ifelse(metadata$RunDay==21,'orange','red'))))

legend(2.8,-1.5,legend = c("15","17","21","23"), ncol=4,text.col = c("blue","green","orange","red"),cex=1)
dev.off()

internalStandard_cv<-mz_grp_cv[4473,]
internalStandard_mean<-mz_grp_mean[4473,]
internalStandard_mean_non_zero<-mz_grp_mean_non_zero[4473,]
internalStandard_zeroes<-mz_grp_zero[4473,]

png("internalStandard.png",height=800,width=800)
par(xpd=TRUE,mar=c(10,4,4,2))
plot(1:52,internalStandard_mean_non_zero,xaxt='n',las=2,xlab="",ylab="values",pch=1,main=paste("Internal standard replicates-fmoc"," ",row.names(ms_data)[4473]),
     col=ifelse(metadata_strains$RunDay==15,'blue',
                ifelse(metadata_strains$RunDay==17,'green',
                       ifelse(metadata_strains$RunDay==21,'orange','red'))))

points(1:52,internalStandard_mean,pch=15,col=ifelse(metadata_strains$RunDay==15,'blue',
                                                    ifelse(metadata_strains$RunDay==17,'green',
                                                           ifelse(metadata_strains$RunDay==21,'orange','red'))))

points(1:52,internalStandard_cv,pch=16,col=ifelse(metadata_strains$RunDay==15,'blue',
                                                  ifelse(metadata_strains$RunDay==17,'green',
                                                         ifelse(metadata_strains$RunDay==21,'orange','red'))))

points(1:52,internalStandard_zeroes,pch=17,col=ifelse(metadata_strains$RunDay==15,'blue',
                                                      ifelse(metadata_strains$RunDay==17,'green',
                                                             ifelse(metadata_strains$RunDay==21,'orange','red'))))

axis(1, at=1:52, labels=colnames(mz_grp_mean_non_zero)[1:52],las=3)
legend(2.8,-2.5,legend = c("15","17","21","23"), ncol=2,text.col = c("blue","green","orange","red"),cex=1)
legend(9.8,-2.5, legend = c("Mean without zero","Mean","CV","Zeroes"), ncol=2,pch = c(1,15,16,17),cex=1)
dev.off()

png("internalStandard_reps_boxplot.png",width=2000)
boxplot(log(as.numeric(ms_data[4473,25:311]))~as.factor(SampleGroup[25:311]),las=2,cex=0.7)
dev.off()


groups<-unique(SampleGroup)

internalStandard_with_outliers<-0
internalStandard_outlier_replicates<-0
a<-1

for(i in 1:length(groups))
{
  ms_data_grp1<-ms_data[4473, grep(groups[i], names(ms_data))]
  ms_data_grp1<-log(ms_data_grp1)
  boxplot_cv<-boxplot(as.numeric(ms_data_grp1),plot=FALSE)
  #print(paste0(i," ",boxplot_cv$out))
  if(length(boxplot_cv$out)>0)
  {
    internalStandard_with_outliers[a]<-groups[i]
    internalStandard_outlier_replicates[a]<-list(names(ms_data_grp1)[which(ms_data_grp1%in%boxplot_cv$out)])
    a<-a+1
  }
}
internalStandard_outlier_replicates<-unlist(internalStandard_outlier_replicates)

# > internalStandard_with_outliers
# [1] "D12_001"  "D12_84"   "D12_94"   "D12_207"  "D12_245"  "D12_255"  "D4_14"    "D4_177"   "D4_187"   "D4_207"   "D4_245"   "D4_253"   "D4_254"   "D4_258"   "D4_325"  
# [16] "matri_02"
# 
# > internalStandard_outlier_replicates
# [1] "D12_001b1_r001" "D12_84b3_r002"  "D12_94b1_r001"  "D12_207b1_r001" "D12_207b2_r002" "D12_245b1_r002" "D12_255b2_r002" "D4_14b3_r001"   "D4_177b1_r001"  "D4_187_b3_r002"
# [11] "D4_207_b2_r001" "D4_245b1_r001"  "D4_253b1_r002"  "D4_254b2_r001"  "D4_258b3_r002"  "D4_325_b3_r002" "matri_02_r004" 

######### Converting values based on 50% zero count per sample group

new_ms_data<-ms_data_tst[, lapply(.SD, function(v) { 
  len <- length(v)
  if((sum(v< 1e-3)/len)>0.5) rep(0,len) else v
}), by=V2]

rownames(new_ms_data)<-colnames(ms_data)
new_ms_data<-as.data.frame(t(new_ms_data))


########################################
###### Multivariate data analysis ######
########################################

#
png("heatmap_fulldata_day12.png",height=800,width=800)
heatmap.2(as.matrix(mz_grp_zero_day12),dendrogram="col",scale="none",main="day12-full data-zerocount",col = brewer.pal(6,"Accent"),trace="none",mar=c(10,10),Rowv=FALSE,labRow= NA)
dev.off()
#heatmap.2(as.matrix(mz_grp_zero_day4),dendrogram="none",trace="none",mar=c(10,10),Rowv=FALSE,Colv=FALSE,labRow= NA, ylab=NULL,)

######## Multinomial test

# observed <- c(4,1)   # observed data: 5 items in category one, 2 items in category two, 1 item in category three
# prob <- c(0.8,0.2)# model: hypothetical probability that an item falls into category one, two, or three
# out <- multinomial.test(observed, prob)

multinomial_prob_value<-0
multinomial_prob_pval<-0
feature_probability<-list()

for (i in 1:5)#nrow(mz_grp_total))
{
  mz_observed<-as.numeric(mz_grp_observed[i,5:26]) #for day 12, it is 5:28, day 4 it is 27:28
  mz_probability<-as.numeric(c(mz_grp_total[i,5:26]/sum(mz_grp_total[i,5:26])))
  multinomial_prob_value<-multinomial.test(mz_observed,mz_probability,MonteCarlo=TRUE)
  feature_probability[i]<-list(c(row.names(mz_grp_total)[i],multinomial_prob_value$p.value))
  multinomial_prob_pval[i]<-multinomial_prob_value$p.value
}

########ttest between columns 

combos <- combn(ncol(ms_data),2)

adply(combos, 2, function(x) {
  test <- t.test(ms_data_test[, x[1]], ms_data_test[, x[2]],paired=F,var.equal=F,p.adjust="BH")
  
  out <- data.frame("var1" = colnames(ms_data)[x[1]]
                    , "var2" = colnames(ms_data[x[2]])
                    , "t.value" = sprintf("%.3f", test$statistic)
                    ,  "df"= test$parameter
                    ,  "p.value" = sprintf("%.3f", test$p.value)
  )
  return(out)
  
})

#t.test(log(ms_data[,1]),log(ms_data[,5]),paired=F,var.equal=F,p.adjust="BH")
#pairwise_ttest_strains<-pairwise.t.test(t(ms_data_test),StrainName_test,p.adj="BH", paired=F, var.equal=T)
#pairwise.wilcox.test(t(ms_data_test),p.adj="BH", paired=F, var.equal=T)

## ttest on rows
t.list <- apply(ms_data,1,function(x){t.test(x[1:6],x[61:66],paired=TRUE,var.equal=F,p.adjust=BH)$p.value}) 
a.list <- apply(ms_data,1,function(x){wilcox.test(x[1:42],x[43:84],paired=TRUE,var.equal=F,p.adj=BH,exact=F)$p.value})

id.sig <- which(a.list < 0.01 );
metab.sig<-names(id.sig)
ttest.table <- cbind(ms_data[metab.sig, ])

#fold change
fc_val<-apply(ttest.table,1,function(x){foldchange(mean(x[1:42]),mean(x[43:84]))})
fc_sig<- which(fc_val>10);

fc.metab.sig<-names(fc_sig)
significant.metab <- cbind(ms_data[fc.metab.sig,],fc_val[fc_sig])

write.table(significant.metab,"comparison_highVslow_DAY12.txt",sep='\t',col.names=NA,quote=FALSE)

# Writing to a file

write.table(pairwise_ttest_strains$p.value,"paiwise_day4.txt",sep='\t',col.names=NA,quote=FALSE)
day4_strains<-read.table("paiwise_day4.txt",sep='\t',header=TRUE,row.names=1)

#apply(d, 2, function(d) {t.test(x = d[,1], y = d[,2])$p.value})

########## Identifying differential features ######

#day4-complete cases

day4_sig_features<-apply(ms_data_day4_nonzero,1,function(x){kruskal.test(x ~ as.factor(SampleGroup_day4))$p.value})
#day4_sig_features<-apply(ms_data_day4_nonzero,1,function(x){kruskal(y=x,trt=as.factor(RunDay_day4),group=TRUE,p.adj="bon")$statistics$p.chisq})
id.sig_day4 <- which(day4_sig_features < 0.05 );
metab.sig_day4<-cbind(ms_data_day4_nonzero[id.sig_day4,],round(day4_sig_features[id.sig_day4],5))

#day12-complete cases

day12_sig_features<-apply(ms_data_day12_nonzero,1,function(x){kruskal.test(x ~ as.factor(SampleGroup_day12))$p.value})
#day12_sig_features<-apply(ms_data_day12_nonzero,1,function(x){kruskal(y=x,trt=as.factor(RunDay_day12),group=TRUE,p.adj="bon")$statistics$p.chisq})
id.sig_day12 <- which(day12_sig_features < 0.05 );
metab.sig_day12<-cbind(ms_data_day12_nonzero[id.sig_day12,],round(day12_sig_features[id.sig_day12],5))

#plotting results

#pval_bins<-c(0,0.001,0.01,0.05,1)

png('pval_sig_strain.png',width=800,height=800)
par(mfrow=c(2,1))
hist(day4_sig_features,breaks=20,main="day4",freq=TRUE)
hist(day12_sig_features,breaks=20,main="day12",freq=TRUE)
dev.off()

mz_rt_day4 <- strsplit(rownames(metab.sig_day4), "\\@")
mz_day4<-sapply(mz_rt_day4 , function (x) if(length(x) == 2) x[1] else as.character(NA))
rt_day4<-sapply(mz_rt_day4 , function (x) if(length(x) == 2) x[2] else as.character(NA))

mz_rt_day12 <- strsplit(rownames(metab.sig_day12), "\\@")
mz_day12<-sapply(mz_rt_day12 , function (x) if(length(x) == 2) x[1] else as.character(NA))
rt_day12<-sapply(mz_rt_day12 , function (x) if(length(x) == 2) x[2] else as.character(NA))

png('mzrt_sig_features_strain_kruskal.png',height=800,width=800)
plot(rt,mz,pch=1,main="sig featuresVsStrain(0.05)",ylab="m/z",xlab="rt",col="#00000033")
points(rt_day4,mz_day4,pch=1,col="green")
points(rt_day12,mz_day12,pch=1,col="red")
legend("topright", title="Colors", c("full data","day4","day12"), fill=c("#00000033","green","red"), cex=0.5)
dev.off()

png("sig_features_strain_kruskal.png",height=800,width=800)
par(mfrow=c(2,2))
hist(as.numeric(mz_day12),main="mz-day12",xlab=paste("sig n:",length(mz_day12),"total=",nrow(ms_data_day12_nonzero)))
hist(as.numeric(rt_day12),main="rt-day12")
hist(as.numeric(mz_day4),main="mz-day4",xlab=paste("sig n:",length(mz_day4),"total=",nrow(ms_data_day4_nonzero)))
hist(as.numeric(rt_day4),main="rt-day4")
#mtext("Kruskal-featuresVsStrains",line=30,outer=TRUE)
dev.off()


########## Clustering ######

d <- dist(t(mz_grp_mean_non_zero), method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
#pdf('hclust.pdf',width=800)
png('hclust_mz_grp_non_zero.png',width=800)
plot (fit)
dev.off()

clusDendro<-as.dendrogram(fit,hang=2)
labelColors<-c("black","red","blue")
clusMember<-rep(1,length(names(new_ms_data1)))
clusMember[grep("D12",names(new_ms_data1))]<-2
clusMember[grep("D4",names(new_ms_data1))]<-3
names(clusMember)<-names(new_ms_data1)

colLab <- function(n)
{
  if(is.leaf(n)) {
    a <- attributes(n)
    # clusMember - a vector designating leaf grouping
    # labelColors - a vector of colors for the above grouping
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}

clusDendro<-dendrapply(clusDendro, colLab)

png('hclust_reps_blanks_matrix_rem5rep.png',width=8000,height=1200)
plot(clusDendro)
dev.off()

#### PCOA

bray_grps.D <- vegdist(t(log(ms_data_day12_nonzero)), "bray") #calculating distance matrix using bray curtis
# res_grps <- pcoa(bray_grps.D)
# res_grps$values
#PCoA.res<-capscale(bray_grps.D~1,distance="bray") 

pc<-cmdscale(bray_grps.D, eig=TRUE, add=TRUE, x.ret =TRUE) 
scores_pcoa<-as.data.frame(pc$x)
eig<-eigenvals(pc) 
var_exp<-cumsum(eig)/sum(eig)
var_exp[1:3]

#day4 non zero  cumsum 0.1906465 0.3345288 0.4094569
#day4           cumsum 0.2041412 0.3663344 0.4326650
#day12 non zero cumsum 0.2075193 0.3309860 0.4153886
#day12          cumsum 0.2972361 0.3891566 0.4529161

#number of zeroes play a huge role in differentiating between day12 and day12non_zero

#eig/sum(eig)

plot(cumsum(pc$eig) / sum(pc$eig), type="h", lwd=5, las=1,xlab="Number of dimensions",ylab=expression(R^2))
plot(pc$eig, type="h", lwd=5, las=1, xlab="Number of dimensions", ylab="Eigenvalues")

kruskal_pcoa_scores<-sapply(scores_pcoa, function(x) kruskal.test(x~ as.factor(RunDay)) )
aov_pcoa_scores<-sapply(scores_pcoa, function(x) aov(x~ as.factor(RunDay)) )

lm_pcoa_scores<-sapply(scores_pcoa, function(x) 
{
  lm_val<-lm(x~ as.factor(RunDay[24:268])) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})

## for day 4
bray_grps.D_day4 <- vegdist(t(log(ms_data_day4)), "bray") #calculating distance matrix using bray curtis
pc_day4<-cmdscale(bray_grps.D_day4, eig=TRUE, add=TRUE, x.ret =TRUE) 
scores_pcoa_day4<-as.data.frame(pc_day4$x)
eig<-eigenvals(pc_day4)
cumsum(eig/sum(eig))
kruskal_pcoa_scores_day4<-sapply(scores_pcoa_day4, function(x) kruskal.test(x~ as.factor(RunDay_day4)) )
aov_pcoa_scores_day4<-sapply(scores_pcoa_day4, function(x) aov(x~ as.factor(RunDay_day4)) )

cor_pcoa_scores_day4<-sapply(scores_pcoa_day4, function(x) cor.test(x,RunDay_day4,method="kendall",use="all.obs"))

lm_pcoa_scores_day4<-sapply(scores_pcoa_day4, function(x) 
{
  lm_val<-lm(x~ as.factor(RunDay_day4)) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})

## for day 12
bray_grps.D_day12 <- vegdist(t(log(ms_data_day12)), "bray") #calculating distance matrix using bray curtis
pc_day12<-cmdscale(bray_grps.D_day12, eig=TRUE, add=TRUE, x.ret =TRUE) 
scores_pcoa_day12<-as.data.frame(pc_day12$x)
eig<-eigenvals(pc_day12)
cumsum(eig/sum(eig))
kruskal_pcoa_scores_day12<-sapply(scores_pcoa_day12, function(x) kruskal.test(x~ as.factor(RunDay_day12)) )
aov_pcoa_scores_day12<-sapply(scores_pcoa_day12, function(x) aov(x~ as.factor(RunDay_day12)) )
cor_pcoa_scores_day12<-sapply(scores_pcoa_day12, function(x) cor.test(x,RunDay_day12,method="kendall",use="all.obs"))

lm_pcoa_scores_day12<-sapply(scores_pcoa_day12, function(x) 
{
  lm_val<-lm(x~ as.factor(RunDay_day12)) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})

# pdf(file="PCOA_samples.pdf")
# biplot(res_grps)
# dev.off()

# Create ordination plot
# pdf(file="PCOA_ordplot_zeroCor_rem5rep.pdf",width=12, height=12, paper="a4r")
# fig<-ordiplot(scores(pc)[,c(1,2)], type="t", main="PCoA Samples",cex=0.5)
# dev.off()

x11()
ordiplot (scores(pc)[,c(1,2)], display = 'sp', type = 'n',main="PCoA Samples", cex=0.5)
points(scores(pc)[,c(1,2)], col = clusMember , pch = clusMember )
legend("bottomleft", legend = c("Day4","Day12"), pch = 1:2,col = c("blue","red"))

#### Testing various correlations ########

# Correlation between PCOA-axis1 and number of zeroes in each column(sample)
# Testing the influence of number of zeroes on separation of samples

pcoa_scores_axis1<-scores(pc)[,c(1)]
cor_zero_pcoa1<-cor.test(pcoa_scores_axis1,zeroes_in_column, use="all.obs", method="kendall")  
#(for complete dataset, including blanks)

cor_zero_pcoa1<-cor.test(pcoa_scores_axis1,zeroes_in_column[25:286], use="all.obs", method="kendall")  
#(without blanks,matrix)

cor_zero_pcoa1<-cor.test(pcoa_scores_axis1,zeroes_in_column_rem[25:281], use="all.obs", method="pearson") 
#(without blanks,matrix and 5 bioreps)

cor_zero_pcoa1<-cor.test(pcoa_scores_axis1,zeroes_in_column_day4, use="all.obs", method="kendall") 
#for day 4 strains
cor_zero_pcoa1<-cor.test(pcoa_scores_axis1,zeroes_in_column_day12, use="all.obs", method="kendall") 
#for day 12 strains


# Correlation PCOA-axis2 and RunDay of each column(sample)
# Testing the influence of runday on separation of samples

pcoa_scores_axis2<-scores(pc)[,c(2)]
cor_runday_pcoa2<-cor.test(pcoa_scores_axis2,as.numeric(RunDay_rem[25:281]), use="all.obs", method="kendall")                
#(without blanks,matrix and 5 bioreps)

cor_runday_pcoa2<-cor.test(pcoa_scores_axis2,as.numeric(RunDay_day4), use="all.obs", method="kendall")
#for day4 strains
cor_runday_pcoa2<-cor.test(pcoa_scores_axis2,as.numeric(RunDay_day12), use="all.obs", method="kendall")                
#for day12 strains

# Correlation PCOA-axis2 and GrowthStage of each column(sample)
# Testing the influence of GrowthStage on separation of samples

pcoa_scores_axis2<-scores(pc)[,c(2)]
cor_runday_pcoa3<-cor.test(pcoa_scores_axis2,as.numeric(GrowthStage_rem[25:281]), use="all.obs", method="kendall") 
#(without blanks,matrix and 5 bioreps)
# Pearson's product-moment correlation

#Correlation between zeroes in column and runday for strains
cor.test(zeroes_in_column[25:286],RunDay_strains,method="kendall",use="all.obs")
# Kendall's rank correlation tau

#Correlation between zeroes in column and runday for strains
cor.test(zeroes_in_column[25:286],as.numeric(GrowthStage_strains),method="kendall",use="all.obs")
# Kendall's rank correlation tau

#Correlation between Growth stage and RunDay for strains
cor.test(as.numeric(GrowthStage_strains),RunDay_strains,method="kendall",use="all.obs")
# Kendall's rank correlation tau

#Correlation between zeroes in column and RunDay for strain at groth stage-day4
cor.test(zeroes_in_column_day4,RunDay_day4,method="kendall",use="all.obs")
# Kendall's rank correlation tau

#Correlation between zeroes in column and RunDay for strain at groth stage-day12
cor.test(zeroes_in_column_day12,RunDay_day12,method="kendall",use="all.obs")
# Kendall's rank correlation tau

###############################################
## PCA using prcomp on samp by mz ###########
###############################################

fit <- princomp(log(ms_data_rem), scale=TRUE)
summary(fit) # print variance accounted for 
d<-loadings(fit) # pc loadings 
screeplot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

fit <- prcomp(log(t(ms_data_rem[24:268])),scale=TRUE) #uses SVD
#mycolors <- c("red", "green", "blue", "magenta", "black")
#plot(fit_day12$x, type="n"); text(fit_day12$x, rownames(fit_day12$x), cex=0.5)
#summary(fit) # print variance accounted for 
#d<-loadings(fit) # pc loadings 
#screeplot(fit,type="lines") # scree plot 
scores_pca_prcomp<-as.data.frame(fit$x)
kruskal_pca_scores<-sapply(scores_pca_prcomp, function(x) kruskal.test(x~ as.factor(RunDay)) )
aov_pca_scores<-sapply(scores_pca_prcomp, function(x) 
{
  aov_model<-aov(x~ as.factor(RunDay)) 
  #aic_value<-extractAIC(aov_model)
  #return(aic_value)
}
)

lm_pca_scores<-sapply(scores_pca_prcomp, function(x) 
{
  lm_val<-lm(x~ as.factor(RunDay_rem[24:268])) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})

# for day 4

d<-as.data.frame(t(ms_data_day4))
d<-d[,apply(d, 2, var, na.rm=TRUE) != 0] #removing mz features which have constant variance

fit_day4 <- prcomp(log(d),scale=TRUE) #uses SVD
scores_pca_prcomp_day4<-as.data.frame(fit_day4$x)
kruskal_pca_scores_day4<-sapply(scores_pca_prcomp_day4, function(x) kruskal.test(x~ as.factor(RunDay_day4)))
aov_pca_scores_day4<-sapply(scores_pca_prcomp_day4, function(x) aov_model<-aov(x~ as.factor(RunDay_day4)))
cor_pca_scores_day4<-sapply(scores_pca_prcomp_day4, function(x) cor.test(x,RunDay_day4,method="kendall",use="all.obs"))

lm_pca_scores_runday4<-apply(scores_pca_prcomp_day4,2, function(x) 
{
  lm_val<-lm(x~ as.factor(RunDay_day4)) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})

lm_pca_scores_strain_day4<-apply(scores_pca_prcomp_day4,2, function(x) 
{
  lm_val<-lm(x~ as.factor(SampleGroup_day4)) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})
# head(sort(lm_pca_scores_day4,decreasing=TRUE))
# PC8      PC11       PC6       PC7       PC9      PC21 
# 0.2152218 0.1933180 0.1894001 0.1865741 0.1360567 0.1352894

#for non-zero
d<-as.data.frame(t(ms_data_day4_nonzero))
d<-d[,apply(d, 2, var, na.rm=TRUE) != 0] #removing mz features which have constant variance
fit_day4 <- prcomp(log(d),scale=TRUE) #uses SVD
scores_pca_prcomp_day4_nonzero<-as.data.frame(fit_day4$x)
lm_pca_scores_runday4_nonzero<-apply(scores_pca_prcomp_day4_nonzero,2, function(x) 
{
  lm_val<-lm(x~ as.factor(RunDay_day4)) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})
lm_pca_scores_strain_day4_nonzero<-apply(scores_pca_prcomp_day4_nonzero,2, function(x) 
{
  lm_val<-lm(x~ as.factor(SampleGroup_day4)) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})
# > head(sort(lm_pca_scores_runday12_nonzero,decreasing=TRUE))
# PC12      PC14       PC2       PC6       PC8      PC20 
# 0.2687165 0.1748026 0.1737839 0.1549156 0.1445072 0.1373533


#for day 12

d<-as.data.frame(t(ms_data_day12))
d<-d[,apply(d, 2, var, na.rm=TRUE) != 0] #removing mz features which have constant variance
fit_day12 <- prcomp(log(d),scale=TRUE) #uses SVD
scores_pca_prcomp_day12<-as.data.frame(fit_day12$x)
kruskal_pca_scores_day12<-sapply(scores_pca_prcomp_day12, function(x) kruskal.test(x~ as.factor(RunDay_day12)))
aov_pca_scores_day12<-sapply(scores_pca_prcomp_day12, function(x) aov_model<-aov(x~ as.factor(RunDay_day12)))
cor_pca_scores_day12<-sapply(scores_pca_prcomp_day12, function(x) cor.test(x,RunDay_day12,method="kendall",use="all.obs"))

lm_pca_scores_runday12<-apply(scores_pca_prcomp_day12,2, function(x) 
{
  lm_val<-lm(x~ as.factor(RunDay_day12)) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})
lm_pca_scores_strain_day12<-apply(scores_pca_prcomp_day12,2, function(x) 
{
  lm_val<-lm(x~ as.factor(SampleGroup_day12)) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})

# head(sort(lm_pca_scores_day12,decreasing=TRUE))
# PC2       PC8      PC11       PC6      PC12      PC13 
# 0.2133751 0.1975413 0.1872261 0.1438924 0.1260683 0.1243714 

#for non-zero
d<-as.data.frame(t(ms_data_day12_nonzero))
#d<-as.data.frame(ms_data_day12_nonzero)
d<-d[,apply(d, 2, var, na.rm=TRUE) != 0] #removing mz features which have constant variance
fit_day12 <- prcomp(log(d),scale=TRUE) #uses SVD
#matplot(cbind(fit_day12_mzbysam$sdev,fit_day12_sambymz$sdev))
scores_pca_prcomp_day12_nonzero<-as.data.frame(fit_day12$x)
#per_var_exp<-fit_day12$sdev^2/sum(fit_day12$sdev^2)

summary(fit_day12)

lm_pca_scores_runday12_nonzero_loadings<-apply(fit_day12_mzbysam_princomp$loadings,2, function(x) 
{
  lm_val<-lm(x~ as.factor(RunDay_day12)) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})

lm_pca_scores_strain_day12_nonzero_loadings<-apply(fit_day12_mzbysam_princomp$loadings,2, function(x) 
{
  lm_val<-lm(x~ as.factor(SampleGroup_day12)) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})

plot(1:length(lm_pca_scores_runday12_nonzero_loadings),lm_pca_scores_runday12_nonzero_loadings,col="green",ylab="Multiple R2",xlab="PC's",main="Day 12",
     ylim=c(0,max(c(max(lm_pca_scores_runday12_nonzero_loadings),max(lm_pca_scores_strain_day12_nonzero_loadings)))))
points(1:length(lm_pca_scores_strain_day12_nonzero_loadings),lm_pca_scores_strain_day12_nonzero_loadings,col="red")
rug(c(1:length(lm_pca_scores_runday12_nonzero_loadings)), side = 1, col="#00000033")
legend("topright", inset=.05,c("Strain","RunDay"),fill=c("red","green"),cex=0.5)


######

pdf("pcscores.pdf",height=16,width=16)
par(mfrow=c(2,2))
plot(1:length(lm_pca_scores_runday12),lm_pca_scores_runday12,col="green",ylab="Multiple R2",xlab="PC's",main="Day 12",
     ylim=c(0,max(c(max(lm_pca_scores_runday12),max(lm_pca_scores_strain_day12)))))
points(1:length(lm_pca_scores_strain_day12),lm_pca_scores_strain_day12,col="red")
rug(c(1:length(lm_pca_scores_runday12)), side = 1, col="#00000033")
legend("topright", inset=.05,c("Strain","RunDay"),fill=c("red","green"),cex=1)

plot(1:length(lm_pca_scores_runday12_nonzero),lm_pca_scores_runday12_nonzero,col="green",ylab="Multiple R2",xlab="PC's",main="Day 12(non zero)",
     ylim=c(0,max(c(max(lm_pca_scores_runday12_nonzero),max(lm_pca_scores_strain_day12_nonzero)))))
points(1:length(lm_pca_scores_strain_day12_nonzero),lm_pca_scores_strain_day12_nonzero,col="red")
rug(c(1:length(lm_pca_scores_runday12_nonzero)), side = 1, col="#00000033")
legend("topright", inset=.05,c("Strain","RunDay"),fill=c("red","green"),cex=1)

plot(1:length(lm_pca_scores_runday4),lm_pca_scores_runday4,col="green",ylab="Multiple R2",xlab="PC's",main="Day 4",
     ylim=c(0,max(c(max(lm_pca_scores_runday4),max(lm_pca_scores_strain_day4)))))
points(1:length(lm_pca_scores_strain_day4),lm_pca_scores_strain_day4,col="red")
rug(c(1:length(lm_pca_scores_runday4)), side = 1, col="#00000033")
legend("topright", inset=.05,c("Strain","RunDay"),fill=c("red","green"),cex=1)

plot(1:length(lm_pca_scores_runday4_nonzero),lm_pca_scores_runday4_nonzero,col="green",ylab="Multiple R2",xlab="PC's",main="Day 4(non zero)",
     ylim=c(0,max(c(max(lm_pca_scores_runday4_nonzero),max(lm_pca_scores_strain_day4_nonzero)))))
points(1:length(lm_pca_scores_strain_day4_nonzero),lm_pca_scores_strain_day4_nonzero,col="red")
rug(c(1:length(lm_pca_scores_runday4_nonzero)), side = 1, col="#00000033")
legend("topright", inset=.05,c("Strain","RunDay"),fill=c("red","green"),cex=1)

dev.off()

# > head(sort(lm_pca_scores_runday12_nonzero,decreasing=TRUE))
# PC12      PC14       PC2       PC6       PC8      PC20 
# 0.2687165 0.1748026 0.1737839 0.1549156 0.1445072 0.1373533

#for non-zero with biochemical parameters

d<-as.data.frame(t(ms_data_day12_nonzero_bc))
d<-d[,apply(d, 2, var, na.rm=TRUE) != 0] #remoing mz features which have constant variance
fit_day12 <- prcomp(log(d),scale=TRUE) #uses SVD
scores_pca_prcomp_day12_nonzero_bc<-as.data.frame(fit_day12$x)

lm_pca_scores_runday12_nonzero_bc<-apply(scores_pca_prcomp_day12_nonzero_bc,2, function(x) 
{
  lm_val<-lm(x~ as.factor(biochem_para$max.biomass.density)) 
  lm_cor<-summary(lm_val)
  return(lm_cor$r.squared)
})

#Strain > head(sort(lm_pca_scores_runday12_nonzero_bc,decreasing=TRUE))
# PC2       PC6       PC5       PC3       PC1       PC7 
# 0.9730438 0.9624849 0.9447193 0.9439467 0.8268231 0.8241385 
#Runday > head(sort(lm_pca_scores_runday12_nonzero_bc,decreasing=TRUE))
# PC3        PC1        PC2        PC4        PC5        PC9 
# 0.83266228 0.48948676 0.40686647 0.10063913 0.02586452 0.02507212 


###############################################
## R functions- batch effect removal ##########
###############################################

#function to compute PCA and perform linear models against runday and strain
#PCA using princomp on mz by samp

compute_pca<-function(dataset,preprocess_method) {
  dataset<-as.data.frame(log(dataset)) #log transform the data using natural log
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
      colnames(dataset_rmbatch1)<-colnames(dataset)
      dataset_rmbatch[[i]]<-dataset_rmbatch1
    }       
  } else{
    dataset_rmbatch<-list()
    svd_dataset<-svd(processed_data)
    svd_dataset$d1<-svd_dataset$d
    svd_dataset$d1[start_pc_comp]<-0 
    end_pc_comp#not used
    dataset_rmbatch1<-svd_dataset$u %*% diag(svd_dataset$d1) %*% t(svd_dataset$v)
    rownames(dataset_rmbatch)<-rownames(dataset)
    colnames(dataset_rmbatch)<-colnames(dataset)
    dataset_rmbatch[[1]]<-dataset_rmbatch1
  }
  return(dataset_rmbatch)
}

# compute permutative f test statistics to difentify signigicant features
compute_perm_ftest<-function(dataset,classlabel) {
  classlabel_factor<-as.numeric(as.factor(classlabel))-1
  if(class(dataset) == "list") {
    p.values<-vector("list", length(dataset)) ### Change to p.values<-vector("list", length(dataset)) # make it faster
    sig_metab_dataset<-vector("list", length(dataset)) 
    for(i in 1:length(dataset))
    { data_matrix<-as.data.frame(dataset[i])
      dataset_sig_features<-mt.maxT(data_matrix,classlabel_factor,test="f",side="abs",fixed.seed.sampling="y",B=100000,nonpara="n")
      p.values[[i]]<-dataset_sig_features
      id.sig_dataset <- which(dataset_sig_features$adjp < 0.05 );
      metab_sig<-cbind(data_matrix[id.sig_dataset,],round(dataset_sig_features$adjp[id.sig_dataset],5))
      sig_metab_dataset[[i]]<-metab_sig
    }
  } else { # only a single dataset
    cat ("not a list")
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
lm_pca_scores_strain_day4_nonzero_loadings<-compute_linearModel(fit_day4_mzbysam_princomp,SampleGroup_day4)
lm_pca_scores_strain_day4_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_strain_day4_nonzero_loadings,"r2")
lm_pca_scores_strain_day4_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_strain_day4_nonzero_loadings,"pval")

lm_pca_scores_runday4_nonzero_loadings<-compute_linearModel(fit_day4_mzbysam_princomp,RunDay_day4)
lm_pca_scores_runday4_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_runday4_nonzero_loadings,"r2")
lm_pca_scores_runday4_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_runday4_nonzero_loadings,"pval")

svd_day4_nonzero <- compute_svd(ms_data_day4_nonzero,"scale",1,ncol(ms_data_day4_nonzero),TRUE)
save(svd_day4_nonzero,file='svd_day4_nonzero.rda')
#str(svd_day4_nonzero)
day4_nonzero_sigfeat_r<-compute_perm_ftest(svd_day4_nonzero,RunDay_day4)
day4_nonzero_sigfeat_r_pval<-day4_nonzero_sigfeat_r[[1]]
day4_nonzero_sigfeat_r_matrix<-day4_nonzero_sigfeat_r[[2]]
day4_nonzero_sigfeat_r_pvaldf<-compute_pval_list(day4_nonzero_sigfeat_r_pval)
day4_no_nonzero_sigfeat_r<-compute_no_features(day4_nonzero_sigfeat_r_matrix)
save(day4_nonzero_sigfeat_r,file='day4_nonzero_sigfeat_r.rda')

day4_nonzero_sigfeat_s<-compute_perm_ftest(svd_day4_nonzero,SampleGroup_day4)
day4_nonzero_sigfeat_s_pval<-day4_nonzero_sigfeat_s[[1]]
day4_nonzero_sigfeat_s_matrix<-day4_nonzero_sigfeat_s[[2]]
day4_nonzero_sigfeat_s_pvaldf<-compute_pval_list(day4_nonzero_sigfeat_s_pval)
day4_no_nonzero_sigfeat_s<-compute_no_features(day4_nonzero_sigfeat_s_matrix)
save(day4_nonzero_sigfeat_s,file='day4_nonzero_sigfeat_s.rda')

###################
#### day12
###################

fit_day12_mzbysam_princomp<-compute_pca(ms_data_day12_nonzero,"scale")

lm_pca_scores_strain_day12_nonzero_loadings<-compute_linearModel(fit_day12_mzbysam_princomp,SampleGroup_day12)
lm_pca_scores_strain_day12_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_strain_day12_nonzero_loadings,"r2")
lm_pca_scores_strain_day12_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_strain_day12_nonzero_loadings,"pval")

lm_pca_scores_runday12_nonzero_loadings<-compute_linearModel(fit_day12_mzbysam_princomp,RunDay_day12)
lm_pca_scores_runday12_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_runday12_nonzero_loadings,"r2")
lm_pca_scores_runday12_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_runday12_nonzero_loadings,"pval")

svd_day12_nonzero <- compute_svd(ms_data_day12_nonzero,"scale",1,ncol(ms_data_day12_nonzero),TRUE)
save(svd_day12_nonzero,file='svd_day12_nonzero.rda')
#str(svd_day12_nonzero)
day12_nonzero_sigfeat_r<-compute_perm_ftest(svd_day12_nonzero,RunDay_day12)
day12_nonzero_sigfeat_r_pval<-day12_nonzero_sigfeat_r[[1]]
day12_nonzero_sigfeat_r_matrix<-day12_nonzero_sigfeat_r[[2]]
day12_nonzero_sigfeat_r_pvaldf<-compute_pval_list(day12_nonzero_sigfeat_r_pval)
day12_no_nonzero_sigfeat_r<-compute_no_features(day12_nonzero_sigfeat_r_matrix)
save(day12_nonzero_sigfeat_r,file='day12_nonzero_sigfeat_r.rda')

day12_nonzero_sigfeat_s<-compute_perm_ftest(svd_day12_nonzero,SampleGroup_day12)
day12_nonzero_sigfeat_s_pval<-day12_nonzero_sigfeat_s[[1]]
day12_nonzero_sigfeat_s_matrix<-day12_nonzero_sigfeat_s[[2]]
day12_nonzero_sigfeat_s_pvaldf<-compute_pval_list(day12_nonzero_sigfeat_s_pval)
day12_no_nonzero_sigfeat_s<-compute_no_features(day12_nonzero_sigfeat_s_matrix)
save(day12_nonzero_sigfeat_s,file='day12_nonzero_sigfeat_s.rda')

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


day12_nonzero_sigfeat_r_pval<-day12_nonzero_sigfeat_r[[1]]
day12_nonzero_sigfeat_r_matrix<-day12_nonzero_sigfeat_r[[2]]
day12_nonzero_sigfeat_r_pvaldf<-compute_pval_list(day12_nonzero_sigfeat_r_pval)
day12_no_nonzero_sigfeat_r<-compute_no_features(day12_nonzero_sigfeat_r_matrix)

png("pval_factor_new.png",height=800,width=800)
par(mfrow=c(2,1))
matplot(day4_x138_nonzero_sigfeat_r_pvaldf,type='l',ylab="pval",xlab="day4 Runday")
matplot(day4_x138_nonzero_sigfeat_s_pvaldf,type='l',ylab="pval",xlab="day4 Strain")
dev.off()

#matplot(cbind(fit_day12_mzbysam_princomp$sdev,svd_day12_nonzero$d))

#################### Plots for batch effect correction ####################

pdf("LinearModel_PCStrainRunday.pdf",width=16)
m<- matrix(c(1,2,3,4,5,5),ncol = 2,byrow = TRUE,heights = c(0.5,0.5,0.1))
layout(mat = m)

plot(lm_pca_scores_strain_day4_nonzero_loadings_r.sq,type="l", xaxt="n",xlab="PC's",ylab="mulriple r2",main="day4")
points(lm_pca_scores_runday4_nonzero_loadings_r.sq,type="l",lty=2,col="red")
axis(1, at = seq(0, 125, by = 5))

plot(lm_pca_scores_strain_day12_nonzero_loadings_r.sq,type="l",xaxt="n",xlab="PC's",ylab="mulriple r2",main="day12")
points(lm_pca_scores_runday12_nonzero_loadings_r.sq,type="l",lty=2,col="red")
axis(1, at = seq(0, 125, by = 5))

plot(lm_pca_scores_strain_day4_nonzero_loadings_pval,type="l",xaxt="n",xlab="PC's",ylab="pval",main="day4")
points(lm_pca_scores_runday4_nonzero_loadings_pval,type="l",lty=2,col="red")
axis(1, at = seq(0, 125, by = 5))

plot(lm_pca_scores_strain_day12_nonzero_loadings_pval,type="l",xaxt="n",xlab="PC's",ylab="pval",main="day12")
points(lm_pca_scores_runday12_nonzero_loadings_pval,type="l",lty=2,col="red")
axis(1, at = seq(0, 125, by = 5))

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("Strains", "RunDay"),  col=c("black","red"), lty=c(1,2),horiz = TRUE)

dev.off()

pdf("BatchCorrectedPermuteTest_FeatureVsStrainRunday.pdf",width=12,height=12)
m<- matrix(c(1,2,3,4,5,5),ncol = 2,byrow = TRUE)
layout(mat = m)

plot(1:length(day4_sig_features_nc$adjp),day4_sig_features_nc$adjp,main="day4 strains",xlab="features",ylab="adj. pval")
points(1:length(day4_sig_features$adjp),day4_sig_features$adjp,col="red")

plot(1:length(day4_sig_features_nc_r$adjp),day4_sig_features_nc_r$adjp,main="day4 runday",xlab="features",ylab="adj. pval")
points(1:length(day4_sig_features_r$adjp),day4_sig_features_r$adjp,col="red")

plot(1:length(day12_sig_features_nc$adjp),day12_sig_features_nc$adjp,main="day12 strains",xlab="features",ylab="adj. pval")
points(1:length(day12_sig_features$adjp),day12_sig_features$adjp,col="red")

plot(1:length(day12_sig_features_nc_r$adjp),day12_sig_features_nc_r$adjp,main="day12 runday",xlab="features",ylab="adj. pval")
points(1:length(day12_sig_features_r$adjp),day12_sig_features_r$adjp,col="red")

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("Non-corrected", "Runday effect Corrected"),   text.col=c("black","red"), horiz = TRUE)
dev.off()


#################### Correlation with environmental parameters ############

#Percent of vairation explained by biochemical parameters

a<-adonis(t(log(ms_data_day12_nonzero_bc))~ biochem_para$Protein.percen+biochem_para$Carb.percent+
            biochem_para$lipid.produc+biochem_para$max.biomass.density+biochem_para$avg.OD+biochem_para$RunDay+ biochem_para$Strain, 
          method = "bray", perm=999)

# Call:
#   adonis(formula = t(log(ms_data_day12_nonzero_bc)) ~ biochem_para$Protein.percen +biochem_para$Carb.percent + biochem_para$lipid.produc + biochem_para$max.biomass.density +      biochem_para$avg.OD + biochem_para$RunDay + biochem_para$Strain,      permutations = 999, method = "bray") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
# biochem_para$Protein.percen       1  0.003704 0.0037038  9.4517 0.04512  0.001 ***
#   biochem_para$Carb.percent         1  0.005983 0.0059830 15.2679 0.07289  0.001 ***
#   biochem_para$lipid.produc         1  0.009956 0.0099565 25.4077 0.12130  0.001 ***
#   biochem_para$max.biomass.density  1  0.008528 0.0085277 21.7617 0.10389  0.001 ***
#   biochem_para$avg.OD               1  0.009422 0.0094218 24.0433 0.11478  0.001 ***
#   biochem_para$RunDay               1  0.008998 0.0089983 22.9625 0.10962  0.001 ***
#   biochem_para$Strain               4  0.015115 0.0037788  9.6430 0.18414  0.001 ***
#   Residuals                        52  0.020377 0.0003919         0.24825           
# Total                            62  0.082083                   1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


####nn-graphs

ms_graph_data_strains<-t(log1p(ms_data))

rownames(ms_graph_data_strains)<-rownames(ms_graph_data_strains)
ms.d_str<-as.matrix(dist(ms_graph_data_strains,"euclidean"))
ms.1nn_str<-define.1nn.graph(ms.d_str)

pdf(file="NN_graph_strains.pdf",paper="a4r")
plot(define.1nn.graph(ms.d_str),layout=layout.fruchterman.reingold,vertex.label=V(ms.1nn_str)$name,vertex.label.cex=0.8,main="Strains")
dev.off()


################PERMANOVA

#ms_data_total$DayStrain<-paste0(ms_data_total$Day,ms_data_total$Strain)
a<-adonis(t(log(ms_data+1)) ~ GrowthStage*StrainId+RunDay, strata=StrainId  , data=ms_data,method = "bray", perm=999)
b1<-metaMDS(log(ms_data+1),k=2)

plot(b)
ordiplot(b,type="n")
orditorp(b,display="species",col="red",air=0.01)
orditorp(b,display="sites",cex=0.01,air=0.01)

distance<-vegdist(log(t(ms_data+1)), method="bray")
model1<-betadisper(distance, SampleGroup)
permute_data<-permutest(model1, pairwise = TRUE,permutations=999)
posthoc_test<-TukeyHSD(model1)

write.table(permute_data$pairwise,"permute_data_day12",quote=F)
write.table(posthoc_test$group,"posthoc_test_day12",quote=F)

###############Plotting TIC's

#combined

getTICs(files=mzxmlfiles,pdfname=paste("xcms_agilent_TIC_algae_15th.pdf", sep = ""),rt="raw")

#ws

mzxmlfiles1<-mzxmlfiles[1:6] 
getTICs(files=mzxmlfiles1,pdfname=paste("15th_ACN_Blank.pdf", sep = ""),rt="raw")

#tt8
tt8_mzxmlfiles<-mzxmlfiles[-grep('WS',mzxmlfiles)] 
getTICs(files=tt8_mzxmlfiles,pdfname=paste("xcms_agilent_TIC_TT8_140513.pdf", sep = ""),rt="raw")

######### ggplot ########

ggplot(ms_data_nd_new, aes_string(x=rt,y=n)) + geom_line(aes(colour = series)) + facet_grid(series ~ .)+theme(axis.ticks = element_blank(), axis.text.y = element_blank())
lapply(list("00111","12922"), function(i) ggplot(i,aes(x=rt,value))+geom_point())
log_blanks<-log1p(blanks)
log_blanks$rt<-ms_data_total[,5]
log_blanks<- melt(log_blanks,  id = 'rt', variable_name = 'series')
ggplot(log_blanks, aes(x=rt,value)) + geom_line(aes(colour = series)) +geom_rug(sides="b")

######### TIC overlays #############################################

getTIC <- function(file,rtcor=NULL) {
  object <- xcmsRaw(file)
  cbind(if (is.null(rtcor)) object@scantime else rtcor, rawEIC(object,mzrange=range(object@env$mz))$intensity) 
}

##
##  overlay TIC from all files in current folder or from xcmsSet, create pdf
##
getTICs <- function(xcmsSet=NULL,files=NULL, pdfname="TICs.pdf",rt=c("raw","corrected")) {
  if (is.null(xcmsSet)) {
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
    if (is.null(files))
      files <- getwd()
    info <- file.info(files)
    listed <- list.files(files[info$isdir], pattern = filepattern,
                         recursive = TRUE, full.names = TRUE)
    files <- c(files[!info$isdir], listed)
  } else {
    files <- filepaths(xcmsSet)
  }
  
  N <- length(files)
  TIC <- vector("list",N)
  
  for (i in 1:N) {
    cat(files[i],"\n")
    if (!is.null(xcmsSet) && rt == "corrected")
      rtcor <- xcmsSet@rt$corrected[[i]] else 
        rtcor <- NULL
    TIC[[i]] <- getTIC(files[i],rtcor=rtcor)
  }
  
  pdf(pdfname,w=16,h=10)
  cols <- rainbow(N)
  lty = 1:N
  pch = 1:N
  xlim = range(sapply(TIC, function(x) range(x[,1])))
  ylim = range(sapply(TIC, function(x) range(x[,2])))
  plot(0, 0, type="n", xlim = xlim, ylim = ylim, main = "Total Ion Chromatograms", xlab = "Retention Time", ylab = "TIC")
  for (i in 1:N) {
    tic <- TIC[[i]]
    points(tic[,1], tic[,2], col = cols[i], pch = pch[i], type="l")
  }
  legend("topright",paste(basename(files)), col = cols, lty = lty, pch = pch)
  dev.off()
  
  invisible(TIC)
}

######### Averaging samples # ref http://stackoverflow.com/questions/10704344/averaging-every-16-columns-in-r

# Average function 
byapply <- function(x, by, fun, ...)
{
  # Create index list
  if (length(by) == 1)
  {
    nc <- ncol(x)
    split.index <- rep(1:ceiling(nc / by), each = by, length.out = nc)
  } else # 'by' is a vector of groups
  {
    nc <- length(by)
    split.index <- by
  }
  index.list <- split(seq(from = 1, to = nc), split.index)
  
  # Pass index list to fun using sapply() and return object
  sapply(index.list, function(i)
  {
    do.call(fun, list(x[, i], ...))
  })
}


################################## Nearest neighbours algorithm #########################
# Author(s): Shiv
# Version :19092013

#--code for constructing and analysing nearest neighbour graphs
#--started by RW on 9 Decvember 2012--
#--in this instance, we're restricting attention to 1-nn graphs

#--take a (named) distance matrix, for each element, find it's nearest neighbour, construct a graph based on non-redundant edges
define.1nn.graph<-function(Dmat)
{
  #--take a (named) distance matrix, Dmat, for each element, find it's nearest neighbour, construct a graph based on non-redundant edges
  #--Dmat should be class 'matrix' not 'dist'
  #--returns a graph
  
  res=NULL
  
  #--define nodes and targets... 
  nodes<-colnames(Dmat)
  nnodes<-rep("",nrow(Dmat))
  val<-rep("",nrow(Dmat))
  
  #--remove diagonal...
  ndDmat<-Dmat
  diag(ndDmat)<-NA
  
  #--go through set of nodes and find nearest neighbour...
  for(i in (1:length(nodes)))
  {
    curcol<-ndDmat[,nodes[i]]
    nnodes[i]<-nodes[which.min(curcol)]
    val[i]<-sort(curcol, FALSE)[1]
  }
  
  #--make edgelist...and define a directed graph..
  el<-cbind(nodes,nnodes)
  el1<-cbind(nodes,nnodes,val)
  res<-graph.edgelist(el)
  
  return(res)
  
}

#--example
#a<-rnorm(100,0,1)
#b<-matrix(a,20,5)
#rownames(b)<-letters[1:nrow(b)]
#--make a distance matrix...
#bd<-as.matrix(dist(b))
#--and compute 1-nn graph
#b1.1nn<-define.1nn.graph(bd)
#--and plot...
#plot(define.1nn.graph(bd),layout=layout.fruchterman.reingold)

#--Rohan finished editing here on 9 December 2012--



######### Metabolite Identification ################################

aracyc<-read.table("D:/Shiv_data/Shiv_NUS/PhD_Project/Integrating_Omics/Analysis_R/MetabolomicsSoftwareComparison/Aracyc/aracyc_20120827_uniqueMass.txt",sep="\t",header=FALSE)
str(aracyc)

# mzlist<-read.table("D:/Shiv_data/Shiv_NUS/PhD_Project/Integrating_Omics/Analysis_R/MetabolomicsSoftwareComparison/Data/TT8_RawData_Metabolomics/Agilent/MZXML_MS1/RemovedNoisyBiorep/mzlist_fc2p05.txt",header=FALSE)
# 
# mzlist<-mzlist[order(mzlist$V1),]
# ### Finding the closest value
# 
# #Modified by Shiv --from The R Book Second edition pg 47
# 
# aracyc_cpds<-aracyc[,1]
# mz_match <- list()
# 
# ## Function 
# 
# closest<-function(xv,sv)
# {
#   #xv[which(abs(xv-sv)==min(abs(xv-sv)))]
#   ifelse(min(abs(xv-sv))<0.5,xv[which(abs(xv-sv)==min(abs(xv-sv)))],0)
# }
# 
# for(i in 1:length(mzlist))
# {
#   cpd<-closest(aracyc_cpds,mzlist[i])
#   #print(mzlist[i,1])
#   print(paste(mzlist[i],cpd,sep=" "))
#   if(cpd >0)
#   {
#     mz_match[[i]]<-cbind(aracyc[which(aracyc$V1==cpd),],mzlist[i])#aracyc[cpd,]
#   }
# }
# 
# mz_id <- do.call("rbind",mz_match) 
# colnames(mz_id)<-c('AracycMZ_Mass','Metabolite_name','MZ')
# mz_id<-mz_id[,c('MZ','AracycMZ_Mass','Metabolite_name')]
# 
# write.table(mz_id,"mz_id.txt",quote=FALSE,sep="\t",row.name=FALSE)


