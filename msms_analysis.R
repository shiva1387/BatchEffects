################################################
# MSMS Data Analysis in R-Malaysian algae data #
################################################
# Author(s): Shiv
# Version: 10112014
# The different modes represent different regions of mz mass coverage.
# Mode1 is from 50 to 250, mode 2 from 250 to 550 and so on.

#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 

library(ggplot2)
library(data.table)
library(stringr)
library(multtest)
library(Hmisc)
library(gtools)
library(RColorBrewer)

######## reading in MS1 data


setwd('../../AlgaeData/results.from.scelse.cluster.211213/')

load("svd_day4_x138_nonzero.rda")
batch_corrected_mat_d4<-svd_day4_nonzero[[7]]

load("svd_day12_x138_nonzero.rda")
batch_corrected_mat_d12<-svd_day12_nonzero[[4]]

load("../algae-data-objects.RData")

#Ensure all strain ids are in the same format
colnames(batch_corrected_mat_d12)<-as.character(gsub('D12_14','D12_014',colnames(batch_corrected_mat_d12)))
colnames(batch_corrected_mat_d12)<-as.character(gsub('D12_84','D12_084',colnames(batch_corrected_mat_d12)))
colnames(batch_corrected_mat_d12)<-as.character(gsub('D12_87','D12_087',colnames(batch_corrected_mat_d12)))
colnames(batch_corrected_mat_d12)<-as.character(gsub('D12_94','D12_094',colnames(batch_corrected_mat_d12)))
colnames(batch_corrected_mat_d12)<-as.character(gsub('D12_94','D12_094',colnames(batch_corrected_mat_d12)))
colnames(batch_corrected_mat_d4)<-as.character(gsub('D4_14','D4_014',colnames(batch_corrected_mat_d4)))

colnames(ms_data_day12_nonzero)<-as.character(gsub('D12_14','D12_014',colnames(ms_data_day12_nonzero)))
colnames(ms_data_day12_nonzero)<-as.character(gsub('D12_84','D12_084',colnames(ms_data_day12_nonzero)))
colnames(ms_data_day12_nonzero)<-as.character(gsub('D12_87','D12_087',colnames(ms_data_day12_nonzero)))
colnames(ms_data_day12_nonzero)<-as.character(gsub('D12_94','D12_094',colnames(ms_data_day12_nonzero)))
colnames(ms_data_day12_nonzero)<-as.character(gsub('D12_94','D12_094',colnames(ms_data_day12_nonzero)))
colnames(ms_data_day4_nonzero)<-as.character(gsub('D4_14','D4_014',colnames(ms_data_day4_nonzero)))

names(SampleGroup_day12)<-as.character(gsub('D12_14','D12_014',names(SampleGroup_day12)))
names(SampleGroup_day12)<-as.character(gsub('D12_84','D12_084',names(SampleGroup_day12)))
names(SampleGroup_day12)<-as.character(gsub('D12_87','D12_087',names(SampleGroup_day12)))
names(SampleGroup_day12)<-as.character(gsub('D12_94','D12_094',names(SampleGroup_day12)))
names(SampleGroup_day12)<-as.character(gsub('D12_94','D12_094',names(SampleGroup_day12)))
names(SampleGroup_day4)<-as.character(gsub('D4_14','D4_014',names(SampleGroup_day4)))

names(RunDay_day12)<-as.character(gsub('D12_14','D12_014',names(RunDay_day12)))
names(RunDay_day12)<-as.character(gsub('D12_84','D12_084',names(RunDay_day12)))
names(RunDay_day12)<-as.character(gsub('D12_87','D12_087',names(RunDay_day12)))
names(RunDay_day12)<-as.character(gsub('D12_94','D12_094',names(RunDay_day12)))
names(RunDay_day12)<-as.character(gsub('D12_94','D12_094',names(RunDay_day12)))
names(RunDay_day4)<-as.character(gsub('D4_14','D4_014',names(RunDay_day4)))


## The data scaled data is stored here
#This is done as ta subset of this data is used to compare the relationship between strains before and after batch correction
ms_data_day12_nonzero_scale<-ScaleData(ms_data_day12_nonzero)
colnames(ms_data_day12_nonzero_scale)<-colnames(ms_data_day12_nonzero)
rownames(ms_data_day12_nonzero_scale)<-rownames(ms_data_day12_nonzero)

ms_data_day4_nonzero_scale<-ScaleData(ms_data_day4_nonzero)
colnames(ms_data_day4_nonzero_scale)<-colnames(ms_data_day4_nonzero)
rownames(ms_data_day4_nonzero_scale)<-rownames(ms_data_day4_nonzero)

# Adding colnames to names of Sample group
names(SampleGroup_day12)<-colnames(ms_data_day12_nonzero)
names(SampleGroup_day4)<-colnames(ms_data_day4_nonzero)

## Raw data
ms_data_day4_msms_names<-c("D4_001","D4_051","D4_187","D4_253","D4_254","D4_322")
SampleGroup_day4_msms<-SampleGroup_day4[SampleGroup_day4 %in% ms_data_day4_msms_names]#Getting the Sample group and strain names
# unique(SampleGroup_day4_msms)
ms_data_day4_msms<-ms_data_day4_nonzero_scale[,colnames(ms_data_day4_nonzero_scale) %in% names(SampleGroup_day4_msms)]
batch_corrected_mat_d4_msms<-batch_corrected_mat_d4[,colnames(batch_corrected_mat_d4) %in% names(SampleGroup_day4_msms)]

ms_data_day12_msms_names<-c("D12_001","D12_051","D12_187","D12_253","D12_254","D12_322")
SampleGroup_day12_msms<-SampleGroup_day12[SampleGroup_day12 %in% ms_data_day12_msms_names]#Getting the Sample group and strain names
# unique(SampleGroup_day12_msms)
ms_data_day12_msms<-ms_data_day12_nonzero_scale[,colnames(ms_data_day12_nonzero_scale) %in% names(SampleGroup_day12_msms)]
batch_corrected_mat_d12_msms<-batch_corrected_mat_d12[,colnames(batch_corrected_mat_d12) %in% names(SampleGroup_day12_msms)]

####################################### Analysis of distance
## Estimating the relations between strains used in MSMS

############### DAY12

#raw
ms_data_day12_msms_aod<-adonis(t(ms_data_day12_msms)~ SampleGroup_day12_msms, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(ms_data_day12_msms) ~ SampleGroup_day12_msms,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleGroup_day12_msms  5     39261  7852.2  4.4525 0.47104  0.001 ***
#   Residuals              25     44088  1763.5         0.52896           
# Total                  30     83349                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#batch corrected
batch_corrected_mat_d12_msms_aod<-adonis(t(batch_corrected_mat_d12_msms)~ SampleGroup_day12_msms, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(batch_corrected_mat_d12_msms) ~ SampleGroup_day12_msms,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleGroup_day12_msms  5     38042  7608.5  6.5802 0.56823  0.001 ***
#   Residuals              25     28907  1156.3         0.43177           
# Total                  30     66949                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

############### DAY4

#raw
ms_data_day4_msms_aod<-adonis(t(ms_data_day4_msms)~ SampleGroup_day4_msms, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(ms_data_day4_msms) ~ SampleGroup_day4_msms,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleGroup_day4_msms  5     33367  6673.5  3.8113 0.39654  0.001 ***
#   Residuals             29     50779  1751.0         0.60346           
# Total                 34     84146                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#batch corrected
batch_corrected_mat_d4_msms_aod<-adonis(t(batch_corrected_mat_d4_msms)~ SampleGroup_day4_msms, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(batch_corrected_mat_d4_msms) ~ SampleGroup_day4_msms,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model   R2 Pr(>F)    
# SampleGroup_day4_msms  5     21817  4363.4  3.2625 0.36  0.001 ***
#   Residuals             29     38785  1337.4         0.64           
# Total                 34     60602                 1.00           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###########################################################################################################################
######## reading in the MS1 data from (MSMS data)
###########################################################################################################################

setwd("F:/Vinay's Algae Data/Aug2013/data/MSMS")
msms_data_mode1<-read.table("VJ-MSMS_101114_mode1_xcms1_38_R.txt",sep="\t",header=T,row.names=1)
str(msms_data_mode1)

metainfo_mode1<-read.table("metainfo_MSMS_101114_mode1_xcms1_38_R.txt",sep="\t",header=T,row.names=1)
str(metainfo_mode1)

msms_data_mode2<-read.table("VJ-MSMS_101114_mode2_xcms1_38_R.txt",sep="\t",header=T,row.names=1)
str(msms_data_mode2)

metainfo_mode2<-read.table("metainfo_MSMS_101114_mode2_xcms1_38_R.txt",sep="\t",header=T,row.names=1)
str(metainfo_mode2)

msms_data_mode3<-read.table("VJ-MSMS_101114_mode3_xcms1_38_R.txt",sep="\t",header=T,row.names=1)
str(msms_data_mode3)

metainfo_mode3<-read.table("metainfo_MSMS_101114_mode3_xcms1_38_R.txt",sep="\t",header=T,row.names=1)
str(metainfo_mode3)

msms_data_mode4<-read.table("VJ-MSMS_101114_mode4_xcms1_38_R.txt",sep="\t",header=T,row.names=1)
str(msms_data_mode4)

metainfo_mode4<-read.table("metainfo_MSMS_101114_mode4_xcms1_38_R.txt",sep="\t",header=T,row.names=1)
str(metainfo_mode4)

### Samples for which MS/MS data is available
#D4_001 D4_187 D4_253 D4_254 D4_322 D4_051
#D12_001 D12_187 D12_253 D12_254 D12_322 D12_051

#### Formatting dataset

## Analysis

msms_data_mode1_data<-msms_data_mode1[, -grep('pool', colnames(msms_data_mode1))] 
msms_data_mode1_data<-msms_data_mode1_data[, -grep('blank', colnames(msms_data_mode1_data))] 

msms_data_mode2_data<-msms_data_mode2[, -grep('pool', colnames(msms_data_mode2))] 
msms_data_mode2_data<-msms_data_mode2_data[, -grep('blank', colnames(msms_data_mode2_data))] 

msms_data_mode3_data<-msms_data_mode3[, -grep('pool', colnames(msms_data_mode3))] 
msms_data_mode3_data<-msms_data_mode3_data[, -grep('blank', colnames(msms_data_mode3_data))] 

msms_data_mode4_data<-msms_data_mode4[, -grep('pool', colnames(msms_data_mode4))] 
msms_data_mode4_data<-msms_data_mode4_data[, -grep('blank', colnames(msms_data_mode4_data))] 


metainfo_mode1_data<-metainfo_mode1[-grep('pool', rownames(metainfo_mode1)),] 
metainfo_mode1_data<-metainfo_mode1_data[-grep('blank', rownames(metainfo_mode1_data)),] 

metainfo_mode2_data<-metainfo_mode2[-grep('pool', rownames(metainfo_mode2)),] 
metainfo_mode2_data<-metainfo_mode2_data[-grep('blank', rownames(metainfo_mode2_data)),] 

metainfo_mode3_data<-metainfo_mode3[-grep('pool', rownames(metainfo_mode3)),] 
metainfo_mode3_data<-metainfo_mode3_data[-grep('blank', rownames(metainfo_mode3_data)),] 

metainfo_mode4_data<-metainfo_mode4[-grep('pool', rownames(metainfo_mode4)),] 
metainfo_mode4_data<-metainfo_mode4_data[-grep('blank', rownames(metainfo_mode4_data)),] 

### Combined data matrix

## lazy workaround to avoid colnames issue
a<-msms_data_mode1_data;
b<-msms_data_mode2_data;names(b)<-colnames(a);
c<-msms_data_mode3_data;names(c)<-colnames(a);
d<-msms_data_mode4_data;names(d)<-colnames(a);

msms_data<-rbind(a,b,c,d)#only the first mode column names will be used
metainfo_msms_data<-metainfo_mode1_data#as only the first mode's column names are used

##NOTE: ensure that the order of the columns are the same in each dataset before merging

##Day4 samples
msms_data_d4<-msms_data[, grep('tp1', colnames(msms_data))] 
msms_data_d4_scale<-ScaleData(msms_data_d4)
metainfo_msms_data_d4<-metainfo_msms_data[grep('tp1', rownames(metainfo_msms_data)),]

##Day12 samples

msms_data_d12<-msms_data[, grep('tp2', colnames(msms_data))] 
msms_data_d12_scale<-ScaleData(msms_data_d12)
metainfo_msms_data_d12<-metainfo_msms_data[grep('tp2', rownames(metainfo_msms_data)),]

### Functions

ScaleData<-function(data_matrix){
  processed_data<-scale(data_matrix,center=T,scale=T)
  colnames(processed_data)<-colnames(data_matrix)
  rownames(processed_data)<-rownames(data_matrix)
  return(processed_data)
}

### PCA ggplot
pcaPlot<-function(axisA,axisB,Strains,residualVariance,plotTitle) {
  set.seed(1)
  forPlot<-data.frame(PCaxisA = axisA,PCaxisB = axisB, Strains=Strains)
  #Choosing color
  colourCount = length(unique(forPlot$Strains))
  getPalette = colorRampPalette(brewer.pal(length(unique(forPlot$Strains)),"Paired"))
  set.seed(1) #important to set seed so that we obtain the same shapes for strains all the time
  pch_types<-c(15, 16, 17, 18, 25, 8)
  pch_values<-sample(pch_types, length(unique(forPlot$Strains)), replace = TRUE)
  
  plot1<- ggplot(data=forPlot, aes(x=PCaxisA, y=PCaxisB,colour= factor(Strains), shape = factor(Strains))) + geom_point(size=4) #for samples
  plot2<- plot1 +  scale_colour_manual('Strains', values=getPalette(colourCount)) + scale_shape_manual('Strains',values=pch_values)
  plot3<- plot2+ theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1,size=12),axis.text.y=element_text(size=12),
                                    panel.grid.major.x = element_blank(), # to x remove gridlines
                                    panel.grid.major.y = element_blank(), # to y remove gridlines
                                    panel.border = element_blank(),  # remove top and right border
                                    panel.background = element_blank(),
                                    axis.line = element_line(color = 'black'))+ 
    xlab(paste0("PC 1 loadings","\n","Variation exp= ",round(residualVariance[1]*100,2),"%")) + 
    ylab(paste0("PC 2 loadings","\n","Variation exp= ",round(residualVariance[2]*100,2),"%")) +
    ggtitle(plotTitle)
  return(plot3)
}

##Multiple ggplots

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
##Averaging replicates
avg.strains<-function(datamatrix0,strains){
  datamatrix<-as.data.frame(t(datamatrix0))
  datamatrix$strain<-strains
  datamatrix <- data.table(datamatrix)
  datamatrix<-datamatrix[,lapply(.SD, mean),by=strain]
  datamatrix<-as.data.frame(datamatrix)
  rownames(datamatrix)<-datamatrix[,1]
  datamatrix<-datamatrix[,2:ncol(datamatrix)]
  return(datamatrix)
}

####################################### Analysis of distance
## Estimating the relations between strains before and after batch correction procedure

############### DAY4

msms_data_d4_aod<-adonis(t(msms_data_d4_scale)~ metainfo_msms_data_d4$SampleName, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(msms_data_d4_scale) ~ metainfo_msms_data_d4$SampleName,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# metainfo_msms_data_d4$SampleName  5    5531.9 1106.38  3.9194 0.62022  0.001 ***
#   Residuals                        12    3387.4  282.28         0.37978           
# Total                            17    8919.3                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

################ DAY12
msms_data_d12_aod<-adonis(t(msms_data_d12_scale)~ metainfo_msms_data_d12$SampleName, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(msms_data_d12_scale) ~ metainfo_msms_data_d12$SampleName,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# metainfo_msms_data_d12$SampleName  5    6821.8 1364.36  5.6202 0.70076  0.001 ***
#   Residuals                         12    2913.1  242.76         0.29924           
# Total                             17    9734.9                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pdf("aod_plots_MSMS.pdf",height=8,width=12)
par(mfrow=c(2,3))
plot(density(ms_data_day4_msms_aod), xlim=c(0,5), main="DAY4 raw")
plot(density(ms_data_day12_msms_aod), xlim=c(0,5), main="DAY12 raw")
plot(density(msms_data_d4_aod), xlim=c(0,5), main="DAY4 MSMS")
plot(density(batch_corrected_mat_d4_msms_aod), xlim=c(0,5), main="DAY4 corrected")
plot(density(batch_corrected_mat_d12_msms_aod), xlim=c(0,5), main="DAY12 corrected")
plot(density(msms_data_d12_aod), xlim=c(0,5), main="DAY12 MSMS")
dev.off()

####################################### PCA

##Averaging replicates 
ms_data_day4_msms_reps<-t(avg.strains(ms_data_day4_msms,SampleGroup_day4_msms))
ms_data_day12_msms_reps<-t(avg.strains(ms_data_day12_msms,SampleGroup_day12_msms))
batch_corrected_mat_d4_msms_reps<-t(avg.strains(batch_corrected_mat_d4_msms,SampleGroup_day4_msms))
batch_corrected_mat_d12_msms_reps<-t(avg.strains(batch_corrected_mat_d12_msms,SampleGroup_day12_msms))
msms_data_d4_scale_reps<-t(avg.strains(msms_data_d4_scale,metainfo_msms_data_d4$SampleName))
msms_data_d12_scale_reps<-t(avg.strains(msms_data_d12_scale,metainfo_msms_data_d12$SampleName))

ms_data_day4_msms_pca<-princomp(ms_data_day4_msms_reps,cor=F,scores=T)
residual_variance<-ms_data_day4_msms_pca$sdev^2/sum(ms_data_day4_msms_pca$sdev^2)
ms_data_day4_msms_pcaPlot<-pcaPlot(ms_data_day4_msms_pca$loadings[,1],ms_data_day4_msms_pca$loadings[,2],unique(SampleGroup_day4_msms),residual_variance,"Day 4 raw")

ms_data_day12_msms_pca<-princomp(ms_data_day12_msms_reps,cor=F,scores=T)
residual_variance<-ms_data_day12_msms_pca$sdev^2/sum(ms_data_day12_msms_pca$sdev^2)
ms_data_day12_msms_pcaPlot<-pcaPlot(ms_data_day12_msms_pca$loadings[,1],ms_data_day12_msms_pca$loadings[,2],unique(SampleGroup_day12_msms),residual_variance,"Day 12 raw")

batch_corrected_mat_d4_msms_pca<-princomp(batch_corrected_mat_d4_msms_reps,cor=F,scores=T)
residual_variance<-batch_corrected_mat_d4_msms_pca$sdev^2/sum(batch_corrected_mat_d4_msms_pca$sdev^2)
batch_corrected_mat_d4_msms_pcaPlot<-pcaPlot(batch_corrected_mat_d4_msms_pca$loadings[,1],batch_corrected_mat_d4_msms_pca$loadings[,2],unique(SampleGroup_day4_msms),residual_variance,"Day 4 corrected")

batch_corrected_mat_d12_msms_pca<-princomp(batch_corrected_mat_d12_msms_reps,cor=F,scores=T)
residual_variance<-batch_corrected_mat_d12_msms_pca$sdev^2/sum(batch_corrected_mat_d12_msms_pca$sdev^2)
batch_corrected_mat_d12_msms_pcaPlot<-pcaPlot(batch_corrected_mat_d12_msms_pca$loadings[,1],batch_corrected_mat_d12_msms_pca$loadings[,2],unique(SampleGroup_day12_msms),residual_variance,"Day 12 corrected")

msms_data_d4_scale_pca<-princomp(msms_data_d4_scale_reps,cor=F,scores=T)
residual_variance<-msms_data_d4_scale_pca$sdev^2/sum(msms_data_d4_scale_pca$sdev^2)
msms_data_d4_pcaPlot<-pcaPlot(msms_data_d4_scale_pca$loadings[,1],msms_data_d4_scale_pca$loadings[,2],unique(as.factor(metainfo_msms_data_d4$SampleName)),residual_variance,"Day 4 MSMS")

msms_data_d12_scale_pca<-princomp(msms_data_d12_scale_reps,cor=F,scores=T)
residual_variance<-msms_data_d12_scale_pca$sdev^2/sum(msms_data_d12_scale_pca$sdev^2)
msms_data_d12_pcaPlot<-pcaPlot(msms_data_d12_scale_pca$loadings[,1],msms_data_d12_scale_pca$loadings[,2],unique(as.factor(metainfo_msms_data_d12$SampleName)),residual_variance,"Day 12 MSMS")

pdf("PCA_D4-D12_MSMSComparison.pdf",height=12,width=16)
multiplot(ms_data_day4_msms_pcaPlot, batch_corrected_mat_d4_msms_pcaPlot, ms_data_day12_msms_pcaPlot, 
          batch_corrected_mat_d12_msms_pcaPlot,msms_data_d4_pcaPlot, msms_data_d12_pcaPlot, cols=3)
dev.off()