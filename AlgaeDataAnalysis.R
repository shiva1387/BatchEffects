###########################################
# Data Analysis in R-Malaysian algae data #
###########################################
# Author(s): Shiv
# Version: 06102013 
# Input: ".tsv" file from XCMS 
# Software: XCMS
# Modified By :Shivshankar Umashankar 

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
# 
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

#ms_data_day4<-ms_data[, grep('D4', names(ms_data))] 
#ms_data_day12<-ms_data[, grep('D12', names(ms_data))] 

#name_list <- strsplit(SampleGroup, "_")
#GrowthStage<-sapply(name_list , function (x) if(length(x) == 2) x[1] else as.character(NA))
#StrainId<-sapply(name_list , function (x) if(length(x) == 2) x[2] else as.character(NA))
#StrainName<-as.factor(names(ms_data))
#StrainName0<-names(ms_data)


###### Data: inclusion and exclusion, removing replicates(samples with large variance) from the dataset 

rem <- c('D4_094_b1_r002','D4_104_b3_r002','D4_245b1_r001','D4_254b1_r001','D4_325_b1_r002')

ms_data_rem<-ms_data[, !names(ms_data) %in% rem]
names(ms_data_rem[25:281])

metadata_rem<-metadata[, !colnames(metadata) %in% rem]
RunDay_rem<-metadata_rem[3,]
GrowthStage_rem<-metadata_rem[1,]
RunDay_rem<-as.vector(RunDay_rem)
GrowthStage_rem<-as.vector(GrowthStage_rem)

GrowthStage_rem<-gsub('D4', 1, GrowthStage_rem)
GrowthStage_rem<-gsub('D12', 2, GrowthStage_rem)

GrowthStage_strains<-GrowthStage[25:286]
GrowthStage_strains<-gsub('D4', 1, GrowthStage_strains)
GrowthStage_strains<-gsub('D12', 2, GrowthStage_strains)

RunDay_strains<-as.numeric(RunDay[25:286])

z1<-t(zeroes_in_column)
colnames(z1)<-colnames(ms_data)
zeroes_in_column_rem<-z1[, !colnames(z1) %in% rem]

GrowthStage_day4<-GrowthStage[grep('D4', names(GrowthStage))]
GrowthStage_day12<-GrowthStage[grep('D12', names(GrowthStage))] 

RunDay_day4<-as.numeric(RunDay[grep('D4', names(RunDay))])
RunDay_day12<-as.numeric(RunDay[grep('D12', names(RunDay))])

zeroes_in_column_day4<-z1[grep('D4', colnames(z1))] 
zeroes_in_column_day12<-z1[grep('D12', colnames(z1))] 

#Measure of variation between replicates of each strain #based on boxplot

groups<-unique(SampleGroup)

groups_with_outliers<-0
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
    a<-a+1
  }
}
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

## Counting zeros in each column

zero_count_column<-function(x){ d=as.vector(x); sum(d < 1+1e-3) }# here zero is given as 1
zeroes_in_column<-lapply(ms_data,zero_count_column) 
zeroes_in_column<-as.numeric(zeroes_in_column)

## Generating boxplot and cv/mean for each sample group

# #scatterplot matrix for mean,cv,no of zeroes
# d<-data.frame(c(1,2,3,4,5,6),
#               c(round(mean(ms_data_grp1[,1]), 2),round(mean(ms_data_grp1[,2]), 2),round(mean(ms_data_grp1[,3]), 2),round(mean(ms_data_grp1[,4]), 2),round(mean(ms_data_grp1[,5]), 2),round(mean(ms_data_grp1[,6]), 2)),
#               c(round(cv(ms_data_grp1[,1]), 2),round(cv(ms_data_grp1[,2]), 2),round(cv(ms_data_grp1[,3]), 2),round(cv(ms_data_grp1[,4]), 2),round(cv(ms_data_grp1[,5]), 2),round(cv(ms_data_grp1[,6]), 2)),
#               c(zeroes_in_column[no+1],zeroes_in_column[no+2],zeroes_in_column[no+3],zeroes_in_column[no+4],zeroes_in_column[no+5],zeroes_in_column[no+6]))
# colnames(d)<-c('replicates','mean','cv','zeroes')
# 
# plot(d,pch=16,col=rainbow(9),cex=2,main=paste0(groups[i]))

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

# converting to data.table format

ms_data_tst<-data.table(t(log(ms_data)))
ms_data_tst<-cbind(ms_data_tst,as.factor(SampleGroup))
lapply(ms_data_tst,class) #checking col classes

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


### Counting number of zeroes for each m/z feature in each group

zero_count<-function(x) sum(x < 1e-3) 
zero_mz_feature2<-t(ms_data_tst[,lapply(.SD,zero_count),by=V2])
zero_mz_feature<-as.data.frame(zero_mz_feature2[2:nrow(zero_mz_feature2),])
colnames(zero_mz_feature)<-unique(SampleGroup)

write.table(zero_mz_feature,"mz_features_zeroes.txt",quote=FALSE,sep="\t")
mz_grp_zero<-read.table("mz_features_zeroes.txt",sep='\t',header=TRUE,row.name=1)

mz_grp_zero<-data.frame(cbind(mz,rt,mz_grp_zero))

for (i in seq(3,ncol(mz_grp_cv),2))
{
  mz<-as.vector(mz_grp_cv$mz)
  rt<-as.vector(mz_grp_cv$rt)
 ## now plot scatter plot
  png(paste("noofZeros/", i, "_mz.png", sep = ""),height=800,width=1200)
  par(mfrow=c(2,2))
  plot(mz,mz_grp_zero[,i],main=names(mz_grp_zero)[i],xlab="mz",ylab="No of zeros")
  hist(mz_grp_zero[,i],main=names(mz_grp_zero)[i],xlab="no of zeroes",ylab="Features", breaks=6)
  plot(mz,mz_grp_zero[,i+1],main=names(mz_grp_zero)[i+1],xlab="mz",ylab="No of zeros")
  hist(mz_grp_zero[,i+1],main=names(mz_grp_zero)[i+1],xlab="no of zeroes",ylab="Features", breaks=6)
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

#heatmap.2(as.matrix(mz_grp_cv),dendrogram="column",trace="none",mar=c(10,10),Rowv=FALSE,labRow= NULL, ylab=NULL)

######## Multinomial test

observed <- c(4,2)   # observed data: 5 items in category one, 2 items in category two, 1 item in category three
prob <- c(0.6, 0.4) # model: hypothetical probability that an item falls into category one, two, or three
out <- multinomial.test(observed, prob)

# Exact Multinomial Test, distance measure: p
# 
# Events    pObs    p.value
# 7   0.311          1

########ttest

## ttest between columns 

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

########## Clustering ######

d <- dist(t(new_ms_data1), method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
#pdf('hclust.pdf',width=800)
png('hclust_reps_zeroCor.png',width=8000)
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

bray_grps.D <- vegdist(t(log(ms_data_day12)), "bray") #calculating distance matrix using bray curtis
# res_grps <- pcoa(bray_grps.D)
# res_grps$values
pc<-cmdscale(bray_grps.D, k=10, eig=TRUE, add=TRUE, x.ret =TRUE) 
#PCoA.res<-capscale(bray_grps.D~1,distance="bray") 


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
# Pearson's product-moment correlation
# data:  pcoa_scores_axis1 and zeroes_in_column
# t = -123.828, df = 309, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.9920528 -0.9876049
# sample estimates:
# cor 
# -0.9900737 
# Kendall's rank correlation tau 
# -0.8620787

cor_zero_pcoa1<-cor.test(pcoa_scores_axis1,zeroes_in_column[25:286], use="all.obs", method="kendall")  
#(without blanks,matrix)
# Pearson's product-moment correlation
# data:  pcoa_scores_axis1 and zeroes_in_column[25:286]
# t = 40.629, df = 260, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.9108919 0.9442969
# sample estimates:
#   cor 
# 0.9294758 
# Kendall's rank correlation tau 
# 0.7178032

cor_zero_pcoa1<-cor.test(pcoa_scores_axis1,zeroes_in_column_rem[25:281], use="all.obs", method="pearson") 
#(without blanks,matrix and 5 bioreps)
# Pearson's product-moment correlation
# data:  pcoa_scores_axis1 and as.numeric(as.character(zeroes_in_column_v1[25:281]))
# t = 27.7962, df = 255, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.8331291 0.8945476
# sample estimates:
# cor 
# 0.8670966 
# Kendall's rank correlation tau 
# 0.592765


cor_zero_pcoa1<-cor.test(pcoa_scores_axis1,zeroes_in_column_day4, use="all.obs", method="kendall") 
#for day 4 strains
# Kendall's rank correlation tau
# 
# data:  pcoa_scores_axis1 and zeroes_in_column_day4
# z = -12.7119, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
# tau 
# -0.7477437 

cor_zero_pcoa1<-cor.test(pcoa_scores_axis1,zeroes_in_column_day12, use="all.obs", method="kendall") 
#for day 12 strains
# Kendall's rank correlation tau
# 
# data:  pcoa_scores_axis1 and zeroes_in_column_day12
# z = -12.0514, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
# tau 
# -0.7145159 


# Correlation PCOA-axis2 and RunDay of each column(sample)
# Testing the influence of runday on separation of samples

pcoa_scores_axis2<-scores(pc)[,c(2)]
cor_runday_pcoa2<-cor.test(pcoa_scores_axis2,as.numeric(RunDay_rem[25:281]), use="all.obs", method="kendall")                
#(without blanks,matrix and 5 bioreps)
# Pearson's product-moment correlation
#data:  pcoa_scores_axis2 and as.numeric(RunDay_rem[25:281])
# t = 23.5807, df = 255, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.7851964 0.8629379
# sample estimates:
#   cor 
# 0.8280054
# Kendall's rank correlation tau 
# 0.6448284

cor_runday_pcoa2<-cor.test(pcoa_scores_axis2,as.numeric(RunDay_day4), use="all.obs", method="kendall")
#for day4 strains
# Kendall's rank correlation tau
# 
# data:  pcoa_scores_axis2 and as.numeric(RunDay_day4)
# z = 10.3382, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
# tau 
# 0.7000252 

cor_runday_pcoa2<-cor.test(pcoa_scores_axis2,as.numeric(RunDay_day12), use="all.obs", method="kendall")                
#for day12 strains
# Kendall's rank correlation tau
# 
# data:  pcoa_scores_axis2 and as.numeric(RunDay_day12)
# z = -7.9067, p-value = 2.643e-15
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
# tau 
# -0.5436908 


# Correlation PCOA-axis2 and GrowthStage of each column(sample)
# Testing the influence of GrowthStage on separation of samples

pcoa_scores_axis2<-scores(pc)[,c(2)]
cor_runday_pcoa3<-cor.test(pcoa_scores_axis2,as.numeric(GrowthStage_rem[25:281]), use="all.obs", method="kendall") 
#(without blanks,matrix and 5 bioreps)
# Pearson's product-moment correlation
# # data:  pcoa_scores_axis2 and as.numeric(GrowthStage_rem[25:281])
# t = -15.5856, df = 255, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.7562032 -0.6299474
# sample estimates:
# cor 
# -0.6984706 
# Kendall's rank correlation tau 
# -0.570355 

#Correlation between zeroes in column and runday for strains
cor.test(zeroes_in_column[25:286],RunDay_strains,method="kendall",use="all.obs")

# Kendall's rank correlation tau
# 
# data:  zeroes_in_column[25:286] and RunDay_strains
# z = 0.3597, p-value = 0.7191
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#        tau 
# 0.01676437 

#Correlation between zeroes in column and runday for strains
cor.test(zeroes_in_column[25:286],as.numeric(GrowthStage_strains),method="kendall",use="all.obs")
# Kendall's rank correlation tau
# 
# data:  zeroes_in_column[25:286] and as.numeric(GrowthStage_strains)
# z = -2.2307, p-value = 0.0257
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
# tau 
# -0.1130019 

#Correlation between Growth stage and RunDay for strains
cor.test(as.numeric(GrowthStage_strains),RunDay_strains,method="kendall",use="all.obs")
# Kendall's rank correlation tau
# 
# data:  as.numeric(GrowthStage_strains) and RunDay_strains
# z = -9.8696, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
# tau 
# -0.5618739 

#Correlation between zeroes in column and RunDay for strain at groth stage-day4
cor.test(zeroes_in_column_day4,RunDay_day4,method="kendall",use="all.obs")

# Kendall's rank correlation tau
# 
# data:  zeroes_in_column_day4 and RunDay_day4
# z = 1.9832, p-value = 0.04735
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
# tau 
# 0.1343476 
# 

#Correlation between zeroes in column and RunDay for strain at groth stage-day12
cor.test(zeroes_in_column_day12,RunDay_day12,method="kendall",use="all.obs")
# 
# Kendall's rank correlation tau
# 
# data:  zeroes_in_column_day12 and RunDay_day12
# z = -4.6396, p-value = 3.491e-06
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# -0.3192004 

#### PCA

fit <- princomp(new_ms_data, scale=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
screeplot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

fit <- prcomp(new_ms_data,scale=TRUE)
mycolors <- c("red", "green", "blue", "magenta", "black")
plot(fit$x, type="n"); text(fit$x, rownames(fit$x), cex=0.5)

library(psych)
fit <- principal(log1p(ms_data), nfactors=5, rotate="varimax")
fit


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




