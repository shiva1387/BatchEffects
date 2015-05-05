##Plots for Shiv thesis
##Author: Shiv
##Version :22-04-15
#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

### Check for MAY17- D12_104b3_r001.d
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
library(scales)

#############
# User      #
# Specific  #
# Variables #
#############

### Functions
rm(list=ls())
#### Figure 1

PrimerData<-read.table("../../dataAnalysis/Figures/primer/wo_outliers/xcms_138/batchEffectsPaper/pcoa_all.txt",sep="\t",header=T,check.names=FALSE,row.names=1)

day4<-PrimerData[PrimerData$DataType=="Day4",]
day12<-PrimerData[PrimerData$DataType=="Day12",]
blanks<-PrimerData[PrimerData$DataType=="Blank",]
matrix<-PrimerData[PrimerData$DataType=="Matrix",]

test<-rbind(day4,day12)

set.seed(1)
# ggplot 
plot1<- ggplot(data=day4, aes(x=day4$Axis1, y=day4$Axis2, colour= factor(RunDay), shape = factor(RunDay))) + geom_point(size=1)
#plot1<- ggplot(data=forPlot, aes(x=PCaxisA, y=PCaxisB, colour= factor(Strain), shape = factor(Growth))) + geom_point(size=4)
day4_plot<- plot1+ theme_bw() + theme(axis.text.x=element_text(hjust = 1,size=8),axis.text.y=element_text(size=8),
                                  panel.grid.major.x = element_blank(), # to x remove gridlines
                                  panel.grid.major.y = element_blank(), # to y remove gridlines
                                  panel.border = element_blank(),  # remove top and right border
                                  panel.background = element_blank(),
                                  axis.line = element_line(color = 'black'))+ 
xlab(paste0("PCO 1","\n","44.8% of total variation")) + 
ylab(paste0("PCO 2","\n","10% of total variation")) + #changed for pcoa plots
ggtitle("day4")

pdf("PCOA_ggplot.pdf",height=8,width=10)
multiplot(blanks_plot,day4_plot, matrix_plot, day12_plot,cols=2)
dev.off()

## Plotting Exponential phase with colors according to Strain for thesis presentation

colourCount = length(unique(day4$Strain))
getPalette = colorRampPalette(brewer.pal(10,"Paired"))
set.seed(1) #important to set seed so that we obtain the same shapes for strains all the time
pch_types<-c(15, 16, 17, 18, 25, 8)
pch_values<-sample(pch_types, 22, replace = TRUE)

##plotting
plot1<- ggplot(data=day4, aes(x=day4$Axis1, y=day4$Axis2, colour= factor(Strain), shape = factor(Strain))) + geom_point(size=1)
plot2 <- plot1 + scale_colour_manual('Strain', values=getPalette(colourCount)) + scale_shape_manual('Strain',values=pch_values)
day4_plot1<- plot2+ theme_bw() + theme(axis.text.x=element_text(hjust = 1,size=8),axis.text.y=element_text(size=8),
                                      panel.grid.major.x = element_blank(), # to x remove gridlines
                                      panel.grid.major.y = element_blank(), # to y remove gridlines
                                      panel.border = element_blank(),  # remove top and right border
                                      panel.background = element_blank(),
                                      axis.line = element_line(color = 'black'))+ 
  xlab(paste0("PCO 1","\n","19.1% of total variation")) + 
  ylab(paste0("PCO 2","\n","12.3% of total variation")) + #changed for pcoa plots
  ggtitle("day4")

pdf("PCOA_ggplot_thesisPresentation.pdf",height=8,width=10)
multiplot(day4_plot1,day4_plot1, day4_plot1, day4_plot1,cols=2)
dev.off()


#### PCA


ms_data_princomp<-compute_pca_batcheffect(batch_corrected_mat_d4,"norm")
residual_variance<-ms_data_princomp$sdev^2/sum(ms_data_princomp$sdev^2)

plot(ms_data_princomp$loadings[,2]~ms_data_princomp$loadings[,1])
text(ms_data_princomp$loadings[,2]~ms_data_princomp$loadings[,1], labels = colnames(batch_corrected_mat_d4), cex=0.6, pos=4)

### PCA ggplot
forPlot<-data.frame(PCaxisA = ms_data_princomp$loadings[,1],PCaxisB = ms_data_princomp$loadings[,2], 
                    RunDay=day4$RunDay) #Subset SampleGroups as the last 9 are blanks
##plotting
plot1<- ggplot(data=forPlot, aes(x=PCaxisA, y=PCaxisB, colour= factor(RunDay), shape = factor(RunDay))) + geom_point(size=1)
day4_plot1<- plot1+  theme_bw() + theme(axis.text.x=element_text(hjust = 1,size=8),axis.text.y=element_text(size=8),
                                        panel.grid.major.x = element_blank(), # to x remove gridlines
                                        panel.grid.major.y = element_blank(), # to y remove gridlines
                                        panel.border = element_blank(),  # remove top and right border
                                        panel.background = element_blank(),
                                        axis.line = element_line(color = 'black'))+ 
  xlab(paste0("PCO 1","\n",round(residual_variance[1]*100,2),"% of total variation")) + 
  ylab(paste0("PCO 2","\n",round(residual_variance[2]*100,2),"% of total variation")) + #changed for pcoa plots
  ggtitle("batch corrected day4")

pdf("PCOA_batchCorrectedData.pdf",height=10,width=8)
multiplot(day4_plot1,day12_plot1,cols=1)
dev.off()


########### Calculation of AOD statistics for Fig 1

ScaleDataMin<-function(data_matrix){
  processed_data<-scale(data_matrix,center=T,scale=T)
  processed_data<-processed_data-min(processed_data)
  colnames(processed_data)<-colnames(data_matrix)
  rownames(processed_data)<-rownames(data_matrix)
  return(processed_data)
}

### Against RunDay

RunDay_blanks<-as.vector(sapply(names(ms_data_total_blanks[3:26]), function(x) strsplit(x,"[.]")[[1]][1]))
RunDay_blanks<-gsub('X','',RunDay_blanks)

RunDay_matrix<-as.vector(sapply(names(ms_data_total_matrix[3:27]), function(x) strsplit(x,"[_]")[[1]][2]))
RunDay_matrix<-as.vector(sapply(RunDay_matrix, function(x) strsplit(x,"[.]")[[1]][1]))

aod_blanks_nc<-adonis(t(ScaleDataMin(ms_data_total_blanks[,3:26]))~ RunDay_blanks, method = "bray", perm=999)#first 2 columns are mz and rt
# Call:
#   adonis(formula = t(ScaleDataMin(ms_data_total_blanks[, 3:26])) ~      RunDay_blanks, permutations = 999, method = "bray") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# RunDay_blanks  3   0.75686 0.252285  4.0319 0.37686  0.001 ***
#   Residuals     20   1.25145 0.062573         0.62314           
# Total         23   2.00831                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

aod_blanks_nc<-adonis(t(ScaleDataMin(ms_data_total_blanks[,3:26]))~ RunDay_blanks, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(ScaleDataMin(ms_data_total_blanks[, 3:26])) ~      RunDay_blanks, permutations = 999, method = "euclidean") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_blanks  3    124585   41528  3.4989 0.34419  0.001 ***
#   Residuals     20    237383   11869         0.65581           
# Total         23    361968                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

aod_matrix_nc<-adonis(t(ScaleDataMin(ms_data_total_matrix[,3:27]))~ RunDay_matrix, method = "bray", perm=999)
# Call:
#   adonis(formula = t(ScaleDataMin(ms_data_total_matrix[, 3:27])) ~      RunDay_matrix, permutations = 999, method = "bray") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_matrix  3   0.64885 0.21628  3.2711 0.31847  0.001 ***
#   Residuals     21   1.38851 0.06612         0.68153           
# Total         24   2.03736                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

aod_matrix_nc<-adonis(t(ScaleDataMin(ms_data_total_matrix[,3:27]))~ RunDay_matrix, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(ScaleDataMin(ms_data_total_matrix[, 3:27])) ~      RunDay_matrix, permutations = 999, method = "euclidean") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_matrix  3    112424   37475  2.5391 0.26618  0.001 ***
#   Residuals     21    309937   14759         0.73382           
# Total         24    422361                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


aod_d12_nc<-adonis(t(ScaleDataMin(ms_data_day12_nonzero))~ RunDay_day12, method = "bray", perm=999)
# Call:
#   adonis(formula = t(ScaleDataMin(ms_data_day12_nonzero)) ~ RunDay_day12,      permutations = 999, method = "bray") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_day12   2    1.1987 0.59935  14.844 0.20519  0.001 ***
#   Residuals    115    4.6433 0.04038         0.79481           
# Total        117    5.8420                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
aod_d12_nc<-adonis(t(ScaleDataMin(ms_data_day12_nonzero))~ RunDay_day12, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(ScaleDataMin(ms_data_day12_nonzero)) ~ RunDay_day12,      permutations = 999, method = "euclidean") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# RunDay_day12   2     48187   24093  10.051 0.1488  0.001 ***
#   Residuals    115    275657    2397         0.8512           
# Total        117    323844                 1.0000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

aod_d4_nc<-adonis(t(ScaleDataMin(ms_data_day4_nonzero))~ RunDay_day4, method = "bray", perm=999)
# Call:
#   adonis(formula = t(ScaleDataMin(ms_data_day4_nonzero)) ~ RunDay_day4,      permutations = 999, method = "bray") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_day4   3    1.0299 0.34329  9.4705 0.19016  0.001 ***
#   Residuals   121    4.3861 0.03625         0.80984           
# Total       124    5.4159                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
aod_d4_nc<-adonis(t(ScaleDataMin(ms_data_day4_nonzero))~ RunDay_day4, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(ScaleDataMin(ms_data_day4_nonzero)) ~ RunDay_day4,      permutations = 999, method = "euclidean") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_day4   3     36079 12026.4  4.8135 0.10662  0.001 ***
#   Residuals   121    302318  2498.5         0.89338           
# Total       124    338397                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

### Against Strain

aod_d12_nc<-adonis(t(ScaleDataMin(ms_data_day12_nonzero))~ SampleGroup_day12, method = "bray", perm=999)
# Call:
#   adonis(formula = t(ScaleDataMin(ms_data_day12_nonzero)) ~ SampleGroup_day12,      permutations = 999, method = "bray") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs  MeanSqs F.Model    R2 Pr(>F)    
# SampleGroup_day12  21    3.7038 0.176371  7.9187 0.634  0.001 ***
#   Residuals          96    2.1382 0.022273         0.366           
# Total             117    5.8420                  1.000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
aod_d12_nc<-adonis(t(ScaleDataMin(ms_data_day12_nonzero))~ SampleGroup_day12, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(ScaleDataMin(ms_data_day12_nonzero)) ~ SampleGroup_day12,      permutations = 999, method = "euclidean") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleGroup_day12  21    169158  8055.2  4.9991 0.52235  0.001 ***
#   Residuals          96    154685  1611.3         0.47765           
# Total             117    323844                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

aod_d4_nc<-adonis(t(ScaleDataMin(ms_data_day4_nonzero))~ SampleGroup_day4, method = "bray", perm=999)
# Call:
#   adonis(formula = t(ScaleDataMin(ms_data_day4_nonzero)) ~ SampleGroup_day4,      permutations = 999, method = "bray") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# SampleGroup_day4  21    3.0563 0.145537  6.3528 0.56431  0.001 ***
#   Residuals        103    2.3597 0.022909         0.43569           
# Total            124    5.4159                  1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
aod_d4_nc<-adonis(t(ScaleDataMin(ms_data_day4_nonzero))~ SampleGroup_day4, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(ScaleDataMin(ms_data_day4_nonzero)) ~ SampleGroup_day4,      permutations = 999, method = "euclidean") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleGroup_day4  21    143357  6826.5   3.605 0.42363  0.001 ***
#   Residuals        103    195041  1893.6         0.57637           
# Total            124    338397                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#### Figure 2
r2.pval<-read.table("F:/Projects/NUS/Vinay's Algae Data/Aug2013/data/x138_lm_model_loadings_scale.txt",sep="\t",header=T,check.names=FALSE,row.names=1,,stringsAsFactors=FALSE)
r2.pval$PVAL<-p.adjust(as.numeric(r2.pval$PVAL),method="BH") ### Added on 090115 to perform FDR on pvalue associations as nexted model had FDR corrections for pbalues
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
ggsave("PC_associations_ggplot_BHCorrected.pdf",plot1,height=8,width=10)

### modified on jan 14
r2.pval.V1<-r2.pval
r2.pval.V1[r2.pval.V1[,3]>0.05,4]<-paste0("NS-",r2.pval.V1[r2.pval.V1[,3]>0.05,4])
#r2.pval.runday_day4_nonzero[r2.pval.runday_day4_nonzero[,3]> 0.05,5]<-"NS-R2(RunDay)"

r2.pval.V1$DataType = factor(r2.pval.V1$DataType,levels=c("Strain","NS-Strain","RunDay","NS-RunDay"))
r2.pval.V1<-r2.pval.V1[-c(3)]
r2.pval.V1$dataGroup<-r2.pval$DataType #for finding groups for lines
## ggplot
set.seed(1)
# ggplot 
plot1<- ggplot(data=r2.pval.V1, aes(x=PrincipalComponents, y=R2, group=dataGroup,
                                                          colour= factor(DataType), 
                                                          shape = factor(DataType))) +
  geom_point(size=2) +facet_grid(Day~.) +
  scale_colour_manual(values=c("#7CAE00", "#7CAE00","#F8766D", "#F8766D")) + 
  scale_shape_manual(values=c(16,1,17,2)) + #geom_line(size=0) +
  theme_bw() + theme(axis.text.x=element_text(hjust = 1,size=8),axis.text.y=element_text(size=8),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     #panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black'))
# plot2<-plot1+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                             labels = trans_format("log10", math_format(10^.x)))

ggsave("PC_associations_ggplot_linearmodel_V1.pdf",plot1,height=8,width=10)

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


############ Fig 5


for(i in 1:15)#ncol(ms_data_day4_nonzero))
{
  mz_rt_day4 <- strsplit(rownames(ms_data_day4_nonzero), "\\@")
  mz_day4<-sapply(mz_rt_day4 , function (x) if(length(x) == 2) x[1] else as.character(NA))
  rt_day4<-sapply(mz_rt_day4 , function (x) if(length(x) == 2) x[2] else as.character(NA))
  
  sigfeat_day4_bc_r<-strsplit(rownames(day4_nonzero_sigfeat_r_matrix[[i]]), "\\@")
  sigfeat_day4_bc_mz_r<-sapply(sigfeat_day4_bc_r , function (x) if(length(x) == 2) x[1] else as.character(NA))
  sigfeat_day4_bc_rt_r<-sapply(sigfeat_day4_bc_r , function (x) if(length(x) == 2) x[2] else as.character(NA))
  tot_mz_r<-length(sigfeat_day4_bc_mz_r)
  
  sigfeat_day4_bc_s<-strsplit(rownames(day4_nonzero_sigfeat_s_matrix[[i]]), "\\@")
  sigfeat_day4_bc_mz_s<-sapply(sigfeat_day4_bc_s , function (x) if(length(x) == 2) x[1] else as.character(NA))
  sigfeat_day4_bc_rt_s<-sapply(sigfeat_day4_bc_s , function (x) if(length(x) == 2) x[2] else as.character(NA))
  tot_mz_s<-length(sigfeat_day4_bc_mz_s)
  
  resid_var<-round(100-sum(residual_variance_day4_mzbysam_princomp[1:i])*100,2)
  day4_numdenum_r<-day4_nonzero_denum_values_r[,i]/day4_nonzero_num_values_r[,i]
  day4_numdenum_r_ind<-match(sigfeat_day4_bc_r,mz_rt_day4); day4_numdenum_r_ind<-day4_numdenum_r_ind[!is.na(day4_numdenum_r_ind)] #obtaining indices only for the significant features
  day4_numdenum_s<-day4_nonzero_denum_values_s[,i]/day4_nonzero_num_values_s[,i]
  day4_numdenum_s_ind<-match(sigfeat_day4_bc_s,mz_rt_day4); day4_numdenum_s_ind<-day4_numdenum_s_ind[!is.na(day4_numdenum_s_ind)]
  tmp1<-day4_numdenum_r[day4_numdenum_r_ind]
  tmp1<-tmp1[tmp1<0.1]
  tmp2<-day4_numdenum_s[day4_numdenum_s_ind]
  tmp2<-tmp2[tmp2<0.1]
  
  
  #plotting test statistic vs no-of-sig-features runday
  sigfeat_day4_r<-rep('nosig',length(rt_day4))
  sigfeat_day4_r[day4_numdenum_r_ind]<-'sig'
  
  df1_r<-data.frame(as.numeric(day4_nonzero_denum_values_r[,i]),sigfeat_day4_r)
  colnames(df1_r)<-c('denom','sigfeat')
  df1_r$bins <- cut2(df1_r$denom, c(0.025,0.05,0.1,0.5,1,1.5,2,3))
  df2_r<-as.data.frame.matrix(table(df1_r$bins,df1_r$sigfeat))
  #colnames(df2_r)<-c('bin','nosig','sig')
  df2_r$bin<-gsub('\\[|\\)|\\]', '', rownames(df2_r))
  splits<-strsplit(df2_r$bin, ",")
  df2_r$newBin<-as.numeric(sapply(splits , function (x) if(length(x) == 2) x[2] else as.character(NA)))
  
  #plotting test statistic vs no-of-sig-features strain
  sigfeat_day4_s<-rep('nosig',length(rt_day4))
  sigfeat_day4_s[day4_numdenum_s_ind]<-'sig'
  
  df1_s<-data.frame(as.numeric(day4_nonzero_denum_values_s[,i]),sigfeat_day4_s)
  colnames(df1_s)<-c('denom','sigfeat')
  df1_s$bins <- cut2(df1_s$denom, c(0.025,0.05,0.1,0.5,1,1.5,2,3))
  df2_s<-as.data.frame.matrix(table(df1_s$bins,df1_s$sigfeat))
  #colnames(df2_s)<-c('bin','nosig','sig')
  df2_s$bin<-gsub('\\[|\\)|\\]', '', rownames(df2_s))
  splits<-strsplit(df2_s$bin, ",")
  df2_s$newBin<-as.numeric(sapply(splits , function (x) if(length(x) == 2) x[2] else as.character(NA)))
  
  
  #plot
  tiff(paste('batchAnalysis_thesis_V1/day4/',i,'.tiff',sep=""), width = 180, height = 150,units = 'mm', compression="lzw", res=600)
  #pdf(paste('batchAnalysis_thesis_V1/day4/',i,'.pdf',sep=""),height=4,width=7)
  #m<- matrix(c(1,2,3,4,5,6,7,8,9,9),ncol = 2,byrow = TRUE)
#   m<- matrix(c(1,2,3,4),ncol = 2,byrow = TRUE)
#   layout(mat = m)
# 
  tiff(paste('batchAnalysis_thesis_V1/day4/',i,'a.tiff',sep=""), width = 80, height = 75,units = 'mm', compression="lzw", res=600)
   plot(as.numeric(day4_nonzero_fvalues_r[,i]),as.numeric(day4_nonzero_r2F_r[,i]),pch=1,ylab="r2",xlab="fvalue",col=ifelse(sigfeat_day4_r=="sig","red","#00000033"))
  #legend(x="topright",inset=0, legend=c("sig","non-sig"),col=c("red","black"),pch=1)  
  dev.off()  
  tiff(paste('batchAnalysis_thesis_V1/day4/',i,'b.tiff',sep=""), width = 80, height = 75,units = 'mm', compression="lzw", res=600)
  plot(as.numeric(day4_nonzero_fvalues_s[,i]),as.numeric(day4_nonzero_r2F_s[,i]),pch=1,ylab="r2",xlab="fvalue",col=ifelse(sigfeat_day4_s=="sig","red","#00000033"))
  #legend(x="topright",inset=0, legend=c("sig","non-sig"),col=c("red","black"),pch=1)  
  dev.off()  
  tiff(paste('batchAnalysis_thesis_V1/day4/',i,'c.tiff',sep=""), width = 80, height = 75,units = 'mm', compression="lzw", res=600)
  plot(log(df2_r$newBin),df2_r$nosig,type='l', ylim=c(0,max(df2_r$nosig,df2_r$sig)))
  lines(log(df2_r$newBin),df2_r$sig,type='l',col="red")
  #legend(x="topright",inset=0, legend=c("sig","non-sig"),col=c("red","black"),lty=1)  
  dev.off()  
  tiff(paste('batchAnalysis_thesis_V1/day4/',i,'d.tiff',sep=""), width = 80, height = 75,units = 'mm', compression="lzw", res=600)
  plot(log(df2_s$newBin),df2_s$nosig,type='l', ylim=c(0,max(df2_s$nosig,df2_s$sig)))
  lines(log(df2_s$newBin),df2_s$sig,col="red")
  #legend(x="topright",inset=0, legend=c("sig","non-sig"),col=c("red","black"),lty=1)              
  dev.off()
   
} # end of for loop

#### Figure 2 - Eignvector vs. factors using ANOVA, but using the nested model (Ref Rohan email 17-Dec-2014)

compute_nested_linearModel_associations<-function(results.from.pca,StrainId,RunDayId) { #dependent.factor1 is Strain id(sample groups) and dependent.factor2 is RunDay 
  lm_pca_scores<-apply(results.from.pca$loadings,2, function(x) {
    lm_val<-lm(x~ as.factor(RunDayId) + as.factor(RunDayId)/as.factor(StrainId))
    lm_cor<-summary(lm_val)
    p.val_runday_strain<-anova(lm_val)$'Pr(>F)'[1:2]
    fvalue_runday_strain<-anova(lm_val)$'F value'[1:2] # modified by shiv to return r2 for each model term on 09-Jan 2015
    r2_runday_strain<-anova(lm_val)$'Sum Sq'/sum(anova(lm_val)$'Sum Sq') #includes the r2 term for residuals
    r2_runday_strain<-r2_runday_strain[1:2] #includes the r2 term only for runday and runday/strain
    return(list(lm_cor$r.squared,p.val_runday_strain,fvalue_runday_strain,r2_runday_strain))
   })
}

compute_pca_batcheffect<-function(dataset,preprocess_method) {
  pca_results <- princomp(dataset,cor=F,scores=T) ### IMP: choose quantile normalized or scaled data
  return(pca_results)
}

# #function to extract values from a list
# compute.r2.pval<-function(linearmodel_list,r2.pval) {
#   if(r2.pval=="r2") { #WARNING:code implictly assumes r2 is in the first column and p.val in the second
#     return (sapply(linearmodel_list, function(x){as.numeric(x[1])}))
#   } else{
#     return (sapply(linearmodel_list, function(x){x[2]}))
#   }
# }

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



#Day4
fit_day4_mzbysam_princomp<-compute_pca_batcheffect(batch_corrected_mat_d4,"scale") # From AlgaeDataAnalysis_modular.R
residual_variance_day4_mzbysam_princomp<-fit_day4_mzbysam_princomp$sdev^2/sum(fit_day4_mzbysam_princomp$sdev^2)
lm_pca_strain_runday_day4_nonzero_loadings<-compute_nested_linearModel_associations(fit_day4_mzbysam_princomp,SampleGroup_day4,RunDay_day4)
lm_pca_strain_runday_day4_nonzero_loadings_r.sq.model<-extract.variables.pc.associations(lm_pca_strain_runday_day4_nonzero_loadings,"r2model")
lm_pca_strain_runday_day4_nonzero_loadings_pval<-extract.variables.pc.associations(lm_pca_strain_runday_day4_nonzero_loadings,"pvalue")
lm_pca_strain_runday_day4_nonzero_loadings_pval<-do.call(rbind.data.frame, lm_pca_strain_runday_day4_nonzero_loadings_pval)
lm_pca_strain_runday_day4_nonzero_loadings_pval<-as.data.frame(sapply(lm_pca_strain_runday_day4_nonzero_loadings_pval,
                                                         function(x) p.adjust(as.numeric(x),method="BH"))) # FDR correction

lm_pca_strain_runday_day4_nonzero_loadings_fval<-extract.variables.pc.associations(lm_pca_strain_runday_day4_nonzero_loadings,"fvalue")
lm_pca_strain_runday_day4_nonzero_loadings_fval<-do.call(rbind.data.frame, lm_pca_strain_runday_day4_nonzero_loadings_fval)
lm_pca_strain_runday_day4_nonzero_loadings_r2<-extract.variables.pc.associations(lm_pca_strain_runday_day4_nonzero_loadings,"r2")
lm_pca_strain_runday_day4_nonzero_loadings_r2<-do.call(rbind.data.frame, lm_pca_strain_runday_day4_nonzero_loadings_r2)

r2.pval.strain_runday_day4_nonzero<-cbind(rep("Exponential",length(lm_pca_strain_runday_day4_nonzero_loadings_r.sq.model)),
                                          as.numeric(paste0(1:length(lm_pca_strain_runday_day4_nonzero_loadings_r.sq.model))),
                                          lm_pca_strain_runday_day4_nonzero_loadings_r.sq.model,lm_pca_strain_runday_day4_nonzero_loadings_pval,
                                          lm_pca_strain_runday_day4_nonzero_loadings_fval,lm_pca_strain_runday_day4_nonzero_loadings_r2)

colnames(r2.pval.strain_runday_day4_nonzero)<-c("Day","PrincipalComponents",
                                                "R2(model)","P-value(RunDay)","P-value(RunDay/Strain)",
                                                "F-value(RunDay)","F-value(RunDay/Strain)",
                                                "R2(RunDay)","R2(RunDay/Strain)")


#### Splitting runday and strain and then merging them again
#### This is done to add a column called Significant to mark significant R2 values

# RunDay
r2.pval.runday_day4_nonzero<-data.frame(rep("Exponential",length(lm_pca_strain_runday_day4_nonzero_loadings_r.sq.model)),
                                   as.numeric(paste0(1:length(lm_pca_strain_runday_day4_nonzero_loadings_r.sq.model))),
                                   lm_pca_strain_runday_day4_nonzero_loadings_pval[,1],lm_pca_strain_runday_day4_nonzero_loadings_r2[,1],
                                   rep("R2(RunDay)",length(lm_pca_strain_runday_day4_nonzero_loadings_r.sq.model)),stringsAsFactors=FALSE)
colnames(r2.pval.runday_day4_nonzero)<-c("Day","PrincipalComponents","P-value(RunDay)","R2(RunDay)","Significance")
r2.pval.runday_day4_nonzero[r2.pval.runday_day4_nonzero[,3]> 0.05,5]<-"NS-R2(RunDay)"

# Strain
r2.pval.strain_day4_nonzero<-data.frame(rep("Exponential",length(lm_pca_strain_runday_day4_nonzero_loadings_r.sq.model)),
                                        as.numeric(paste0(1:length(lm_pca_strain_runday_day4_nonzero_loadings_r.sq.model))),
                                        lm_pca_strain_runday_day4_nonzero_loadings_pval[,2],lm_pca_strain_runday_day4_nonzero_loadings_r2[,2],
                                        rep("R2(RunDay/Strain)",length(lm_pca_strain_runday_day4_nonzero_loadings_r.sq.model)),stringsAsFactors=FALSE)
colnames(r2.pval.strain_day4_nonzero)<-c("Day","PrincipalComponents","P-value(RunDay/Strain)","R2(RunDay/Strain)","Significance")
r2.pval.strain_day4_nonzero[r2.pval.strain_day4_nonzero[,3]> 0.05,5]<-"NS-R2(RunDay/Strain)"

#Combine runday and strain
r2.pval.runday_day4_nonzero.melt<-melt(r2.pval.runday_day4_nonzero[,c(1,2,4,5)], id.vars = c("Day","PrincipalComponents","Significance"))
r2.pval.strain_day4_nonzero.melt<-melt(r2.pval.strain_day4_nonzero[,c(1,2,4,5)], id.vars = c("Day","PrincipalComponents","Significance"))

r2.pval.strain_runday_day4_nonzero.V1<-rbind(r2.pval.strain_day4_nonzero.melt,r2.pval.runday_day4_nonzero.melt)

###################################
write.table(r2.pval.strain_runday_day4_nonzero,"PC_associations_nested_model_day4.txt",sep="\t",quote=FALSE,row.names=FALSE)

r2.pval.strain_runday_day4_nonzero.melt<-melt(r2.pval.strain_runday_day4_nonzero[,c(1,2,4,5,8,9)],
                                               id.vars = c("Day","PrincipalComponents"))


#Day12

fit_day12_mzbysam_princomp<-compute_pca_batcheffect(batch_corrected_mat_d12,"scale") # From AlgaeDataAnalysis_modular.R
residual_variance_day12_mzbysam_princomp<-fit_day12_mzbysam_princomp$sdev^2/sum(fit_day12_mzbysam_princomp$sdev^2)
lm_pca_strain_runday_day12_nonzero_loadings<-compute_nested_linearModel_associations(fit_day12_mzbysam_princomp,SampleGroup_day12,RunDay_day12)
lm_pca_strain_runday_day12_nonzero_loadings_r.sq.model<-extract.variables.pc.associations(lm_pca_strain_runday_day12_nonzero_loadings,"r2model")
lm_pca_strain_runday_day12_nonzero_loadings_pval<-extract.variables.pc.associations(lm_pca_strain_runday_day12_nonzero_loadings,"pvalue")
lm_pca_strain_runday_day12_nonzero_loadings_pval<-do.call(rbind.data.frame, lm_pca_strain_runday_day12_nonzero_loadings_pval)
lm_pca_strain_runday_day12_nonzero_loadings_pval<-as.data.frame(sapply(lm_pca_strain_runday_day12_nonzero_loadings_pval,
                                                                      function(x) p.adjust(as.numeric(x),method="BH"))) # FDR correction

lm_pca_strain_runday_day12_nonzero_loadings_fval<-extract.variables.pc.associations(lm_pca_strain_runday_day12_nonzero_loadings,"fvalue")
lm_pca_strain_runday_day12_nonzero_loadings_fval<-do.call(rbind.data.frame, lm_pca_strain_runday_day12_nonzero_loadings_fval)
lm_pca_strain_runday_day12_nonzero_loadings_r2<-extract.variables.pc.associations(lm_pca_strain_runday_day12_nonzero_loadings,"r2")
lm_pca_strain_runday_day12_nonzero_loadings_r2<-do.call(rbind.data.frame, lm_pca_strain_runday_day12_nonzero_loadings_r2)

r2.pval.strain_runday_day12_nonzero<-cbind(rep("Stationary",length(lm_pca_strain_runday_day12_nonzero_loadings_r.sq.model)),
                                          as.numeric(paste0(1:length(lm_pca_strain_runday_day12_nonzero_loadings_r.sq.model))),
                                          lm_pca_strain_runday_day12_nonzero_loadings_r.sq.model,lm_pca_strain_runday_day12_nonzero_loadings_pval,
                                          lm_pca_strain_runday_day12_nonzero_loadings_fval,lm_pca_strain_runday_day12_nonzero_loadings_r2)

colnames(r2.pval.strain_runday_day12_nonzero)<-c("Day","PrincipalComponents",
                                                "R2(model)","P-value(RunDay)","P-value(RunDay/Strain)",
                                                "F-value(RunDay)","F-value(RunDay/Strain)",
                                                "R2(RunDay)","R2(RunDay/Strain)")

write.table(r2.pval.strain_runday_day12_nonzero,"PC_associations_nested_model_day12.txt",sep="\t",quote=FALSE,row.names=FALSE)


r2.pval.strain_runday_day12_nonzero.melt<-melt(r2.pval.strain_runday_day12_nonzero[,c(1,2,4,5,8,9)],
                                               id.vars = c("Day","PrincipalComponents"))

#### Splitting runday and strain and then merging them again
#### This is done to add a column called Significant to mark significant R2 values

# RunDay
r2.pval.runday_day12_nonzero<-data.frame(rep("Stationary",length(lm_pca_strain_runday_day12_nonzero_loadings_r.sq.model)),
                                        as.numeric(paste0(1:length(lm_pca_strain_runday_day12_nonzero_loadings_r.sq.model))),
                                        lm_pca_strain_runday_day12_nonzero_loadings_pval[,1],lm_pca_strain_runday_day12_nonzero_loadings_r2[,1],
                                        rep("R2(RunDay)",length(lm_pca_strain_runday_day12_nonzero_loadings_r.sq.model)),stringsAsFactors=FALSE)
colnames(r2.pval.runday_day12_nonzero)<-c("Day","PrincipalComponents","P-value(RunDay)","R2(RunDay)","Significance")
r2.pval.runday_day12_nonzero[r2.pval.runday_day12_nonzero[,3]> 0.05,5]<-"NS-R2(RunDay)"

# Strain
r2.pval.strain_day12_nonzero<-data.frame(rep("Stationary",length(lm_pca_strain_runday_day12_nonzero_loadings_r.sq.model)),
                                        as.numeric(paste0(1:length(lm_pca_strain_runday_day12_nonzero_loadings_r.sq.model))),
                                        lm_pca_strain_runday_day12_nonzero_loadings_pval[,2],lm_pca_strain_runday_day12_nonzero_loadings_r2[,2],
                                        rep("R2(RunDay/Strain)",length(lm_pca_strain_runday_day12_nonzero_loadings_r.sq.model)),stringsAsFactors=FALSE)
colnames(r2.pval.strain_day12_nonzero)<-c("Day","PrincipalComponents","P-value(RunDay/Strain)","R2(RunDay/Strain)","Significance")
r2.pval.strain_day12_nonzero[r2.pval.strain_day12_nonzero[,3]> 0.05,5]<-"NS-R2(RunDay/Strain)"

#Combine runday and strain
r2.pval.runday_day12_nonzero.melt<-melt(r2.pval.runday_day12_nonzero[,c(1,2,4,5)], id.vars = c("Day","PrincipalComponents","Significance"))
r2.pval.strain_day12_nonzero.melt<-melt(r2.pval.strain_day12_nonzero[,c(1,2,4,5)], id.vars = c("Day","PrincipalComponents","Significance"))

r2.pval.strain_runday_day12_nonzero.V1<-rbind(r2.pval.strain_day12_nonzero.melt,r2.pval.runday_day12_nonzero.melt)


###################### Merge d4 and d12

r2.pval.strain_runday_nonzero.V1<- rbind(r2.pval.strain_runday_day4_nonzero.V1, r2.pval.strain_runday_day12_nonzero.V1)
colnames(r2.pval.strain_runday_nonzero.V1)[4]<-"VarType"
colnames(r2.pval.strain_runday_nonzero.V1)[5]<-"VarValues"

cols = gg_color_hue(4)
# [1] orang "#F8766D" green "#7CAE00" blue "#00BFC4" violet "#C77CFF"

r2.pval.strain_runday_nonzero.V1$Significance = factor(r2.pval.strain_runday_nonzero.V1$Significance, 
                                                       levels=c("R2(RunDay/Strain)","NS-R2(RunDay/Strain)",
                                                                "R2(RunDay)","NS-R2(RunDay)"))


## ggplot
set.seed(1)
# ggplot 
plot1<- ggplot(data=r2.pval.strain_runday_nonzero.V1, aes(x=PrincipalComponents, y=VarValues, group=VarType,
                                                          colour= factor(Significance), 
                                                               shape = factor(Significance))) +
                     geom_point(size=2) +facet_grid(Day~.) +
                     scale_colour_manual(values=c("#7CAE00", "#7CAE00","#F8766D", "#F8766D")) + 
                     scale_shape_manual(values=c(16,1,17,2)) + geom_line(size=0.1) +
                     theme_bw() + theme(axis.text.x=element_text(hjust = 1,size=8),axis.text.y=element_text(size=8),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     #panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black'))
# plot2<-plot1+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                             labels = trans_format("log10", math_format(10^.x)))

ggsave("PC_associations_ggplot_nestedmodel_V3(line)_batchCorr.pdf",plot1,height=8,width=10)


### Version used in Jan 2015 with R2 for individual model terms (plots both r2 and p value in 4 panels)

r2.pval.strain_runday_nonzero.melt<-rbind(r2.pval.strain_runday_day4_nonzero.melt,r2.pval.strain_runday_day12_nonzero.melt)
colnames(r2.pval.strain_runday_nonzero.melt)[3]<-"VarType"
colnames(r2.pval.strain_runday_nonzero.melt)[4]<-"VarValues"

r2.pval.strain_runday_nonzero.melt$VarType.r2.pval<-sapply(as.character(r2.pval.strain_runday_nonzero.melt$VarType), function(x) strsplit(as.character(x),"\\(")[[1]][1])
r2.pval.strain_runday_nonzero.melt$VarType.day<-sapply(as.character(r2.pval.strain_runday_nonzero.melt$VarType), function(x) gsub(")","",strsplit(as.character(x),"\\(")[[1]][2]))

r2.pval.strain_runday_nonzero.melt$VarType.r2.pval = factor(r2.pval.strain_runday_nonzero.melt$VarType.r2.pval, levels=c('R2','P-value'))

# > head(r2.pval.strain_runday_nonzero.melt.sig)
# Day PrincipalComponents         VarType    VarValues VarType.r2.pval VarType.day
# 1 Exponential                   1 P-value(RunDay) 6.405984e-25         P-value      RunDay
# 2 Exponential                   2 P-value(RunDay) 1.863475e-66         P-value      RunDay
# 3 Exponential                   3 P-value(RunDay) 4.639174e-01         P-value      RunDay
# 4 Exponential                   4 P-value(RunDay) 7.600919e-13         P-value      RunDay
# 5 Exponential                   5 P-value(RunDay) 8.972560e-44         P-value      RunDay
# 6 Exponential                   6 P-value(RunDay) 3.148959e-42         P-value      RunDay

set.seed(1)
# ggplot 
plot1<- ggplot(data=r2.pval.strain_runday_nonzero.melt, aes(x=PrincipalComponents, y=VarValues, colour= factor(VarType.day), shape = factor(VarType.day))) + geom_point(size=2) +facet_grid(VarType.r2.pval~Day) +
  theme_bw() + theme(axis.text.x=element_text(hjust = 1,size=8),axis.text.y=element_text(size=8),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     #panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black'))
plot2<-plot1+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = trans_format("log10", math_format(10^.x)))

ggsave("PC_associations_ggplot_nestedmodel_V2(log).pdf",plot2,height=8,width=10)

# 
# r2.pval<-read.table("F:/Projects/NUS/Vinay's Algae Data/Aug2013/data/x138_lm_model_loadings_scale.txt",sep="\t",header=T,check.names=FALSE,row.names=1)
# r2.pval.melt<-melt(r2.pval,id.vars = c("DataType", "Day", "PrincipalComponents"))
# colnames(r2.pval.melt)[4]<-"param"
# colnames(r2.pval.melt)[5]<-"param.values"
# set.seed(1)
# # ggplot 
# plot1<- ggplot(data=r2.pval.melt, aes(x=PrincipalComponents, y=param.values, colour= factor(DataType), shape = factor(DataType))) + geom_point(size=2) +facet_grid(param ~ Day) +
#   theme_bw() + theme(axis.text.x=element_text(hjust = 1,size=8),axis.text.y=element_text(size=8),
#                      panel.grid.major.x = element_blank(), # to x remove gridlines
#                      panel.grid.major.y = element_blank(), # to y remove gridlines
#                      #panel.border = element_blank(),  # remove top and right border
#                      panel.background = element_blank(),
#                      axis.line = element_line(color = 'black'))
# ggsave("PC_associations_ggplot.pdf",plot1,height=8,width=10)
# 


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