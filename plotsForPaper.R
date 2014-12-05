##Plots for Shiv thesis
##Author: Shiv
##Version :05-12-14
#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 

library(ggplot2)
library(preprocessCore)
library(data.table)
library(raster)
library(vegan)
library(xcms)
library(RColorBrewer)
library(gtools)
library(pheatmap)
library(GGally)
library(grid)
library(vegan)

#############
# User      #
# Specific  #
# Variables #
#############

### Functions
rm(list=ls())
#### Figure 1

PrimerData<-read.table("F:/Projects/NUS/Vinay's Algae Data/Aug2013/dataAnalysis/Figures/primer/wo_outliers/xcms_138/batchEffectsPaper/pcoa_all.txt",sep="\t",header=T,check.names=FALSE,row.names=1)

day4<-PrimerData[PrimerData$DataType=="Day4",]
day12<-PrimerData[PrimerData$DataType=="Day12",]
blanks<-PrimerData[PrimerData$DataType=="Blank",]
matrix<-PrimerData[PrimerData$DataType=="Matrix",]

test<-rbind(day4,day12)

set.seed(1)
# ggplot 
plot1<- ggplot(data=blanks, aes(x=blanks$Axis1, y=blanks$Axis2, colour= factor(RunDay), shape = factor(RunDay))) + geom_point(size=1)
#plot1<- ggplot(data=forPlot, aes(x=PCaxisA, y=PCaxisB, colour= factor(Strain), shape = factor(Growth))) + geom_point(size=4)
blanks_plot<- plot1+ theme_bw() + theme(axis.text.x=element_text(hjust = 1,size=8),axis.text.y=element_text(size=8),
                                  panel.grid.major.x = element_blank(), # to x remove gridlines
                                  panel.grid.major.y = element_blank(), # to y remove gridlines
                                  panel.border = element_blank(),  # remove top and right border
                                  panel.background = element_blank(),
                                  axis.line = element_line(color = 'black'))+ 
xlab(paste0("PCO 1","\n","44.8% of total variation")) + 
ylab(paste0("PCO 2","\n","10% of total variation")) + #changed for pcoa plots
ggtitle("blanks")

pdf("PCOA_ggplot.pdf",height=8,width=10)
multiplot(blank_plot,day4_plot, matrix_plot, day12_plot,cols=2)
dev.off()

#### Figure 2
r2.pval<-read.table("F:/Projects/NUS/Vinay's Algae Data/Aug2013/data/x138_lm_model_loadings_scale.txt",sep="\t",header=T,check.names=FALSE,row.names=1)
r2.pval.melt<-melt(r2.pval,id.vars = c("DataType", "Day", "PrincipalComponents"))
colnames(r2.pval.melt)[4]<-"param"
colnames(r2.pval.melt)[5]<-"param.values"
set.seed(1)
# ggplot 
plot1<- ggplot(data=r2.pval.melt, aes(x=PrincipalComponents, y=param.values, colour= factor(DataType), shape = factor(DataType))) + geom_point(size=2) +facet_grid(param ~ Day) +
  theme_bw() + theme(axis.text.x=element_text(hjust = 1,size=8),axis.text.y=element_text(size=8),
                                        panel.grid.major.x = element_blank(), # to x remove gridlines
                                        panel.grid.major.y = element_blank(), # to y remove gridlines
                                        #panel.border = element_blank(),  # remove top and right border
                                        panel.background = element_blank(),
                                        axis.line = element_line(color = 'black'))
ggsave("PC_associations_ggplot.pdf",plot1,height=8,width=10)

#### Supplementary Figure
svd_filters<-read.table("F:/Projects/NUS/Vinay's Algae Data/Aug2013/data/SVDFilters_features.txt",sep="\t",header=T,check.names=FALSE,row.names=1)
svd_filters.melt<-melt(svd_filters,id.vars = c("Day", "PrincipalComponents"))
colnames(svd_filters.melt)[3]<-"param"
colnames(svd_filters.melt)[4]<-"param.values"
set.seed(1)
dummy2 <- data.frame(Day = c("Exponential", "Stationary"), Z = c(7, 4)) # a dummy dataset created to plot indicate PC filters used
# ggplot 
plot1<- ggplot(data=svd_filters.melt, aes(x=PrincipalComponents, y=param.values, colour= factor(param), shape = factor(param))) + geom_point(size=2) +facet_grid(Day~.) +
  geom_vline(data = dummy2, aes(xintercept = Z),linetype = "dotted") +
  theme_bw() + theme(axis.text.x=element_text(hjust = 1,size=8),axis.text.y=element_text(size=8),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     #panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black'))
ggsave("SVDFilters_ggplot.pdf",plot1,height=8,width=10)

