## Shiv thesis comments

##Author: Shiv
##Version :04-11-14

### Library

library('ggplot2')
library('RColorBrewer')
library('data.table')

### Reading Data-svd and differential metabolites

setwd('../../AlgaeData/results.from.scelse.cluster.211213/')

load("svd_day4_x138_nonzero.rda")
batch_corrected_mat_d4<-svd_day4_nonzero[[7]]

load("svd_day12_x138_nonzero.rda")
batch_corrected_mat_d12<-svd_day12_nonzero[[4]]

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


##Getting the full dataset
load("../algae-data-objects.RData")
## Complete.dataset (without zeroes)
ms_data_d4_d12<-ms_data
ms_data_d4_d12[ms_data_d4_d12<1+1e-3]<-NA
ms_data_d4_d12<- ms_data_d4_d12[complete.cases(ms_data_d4_d12),]
SampleGroup_d4_d12<-c(SampleGroup_day12,SampleGroup_day4)
RunDay_d4_d12<-c(RunDay_day12,RunDay_day4)

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


### Functions
ScaleData<-function(data_matrix){
  processed_data<-scale(data_matrix,center=T,scale=T)
  processed_data<-processed_data-min(processed_data)
  return(processed_data)
}
compute_linearModel_batchEffect<-function(data_matrix,StrainId,RunDayId) { #dependent.factor1 is Strain id(sample groups) and dependent.factor2 is RunDay 
  lm_pca_scores<-apply(data_matrix,2, function(x) {
    lm_val<-lm(x~ as.factor(RunDayId))
    p.val_runday<-anova(lm_val)$'Pr(>F)'[1]
    lm_val_strain<-lm(lm_val$residuals~ as.factor(StrainId))
    p.val_strain<-anova(lm_val_strain)$'Pr(>F)'[1]
    return(list(p.val_runday,p.val_strain))
  })
}
compute_nested_linearModel<-function(data_matrix,StrainId,RunDayId) { #dependent.factor1 is Strain id(sample groups) and dependent.factor2 is RunDay 
  lm_pca_scores<-apply(data_matrix,2, function(x) {
    lm_val<-lm(x~ as.factor(RunDayId) + as.factor(RunDayId)/as.factor(StrainId))
    p.val_runday_strain<-anova(lm_val)$'Pr(>F)'[1:2]
    residuals_runday_strain<-resid(lm_val)
    #return(list(residuals_runday_strain,p.val_runday_strain))
    return(p.val_runday_strain)
  })
}
compute_nestedModel_significance<-function(data_matrix,StrainId,RunDayId) { #testing independently for strain and runday
  lm_pca_scores<-apply(data_matrix,2, function(x) {
    lm_val<-lm(x~ as.factor(RunDayId))
    p.val_runday<-anova(lm_val)$'Pr(>F)'[1]
    lm_val_strain<-lm(x~ as.factor(StrainId))
    p.val_strain<-anova(lm_val_strain)$'Pr(>F)'[1]
    return(list(p.val_runday,p.val_strain))
  })
}
compute_linearModel_strain_batch<-function(data_matrix,StrainId,RunDayId) { #Linear model
  lm_pca_scores<-apply(data_matrix,2, function(x) {
    lm_val<-lm(x~  as.factor(StrainId) + as.factor(RunDayId))
    p.val_runday_strain<-anova(lm_val)$'Pr(>F)'[1]
    residuals_runday_strain<-resid(lm_val)
    return(list(residuals_runday_strain))
  })
}
#function to extract r2 value from list containing r2 and p.val returned from linear model
compute.r2.pval<-function(linearmodel_list,r2.pval) {
  if(r2.pval=="r2") { #WARNING:code implictly assumes r2 is in the first column and p.val in the second
    return (sapply(linearmodel_list, function(x){as.numeric(x[1])}))
  } else{
    return (sapply(linearmodel_list, function(x){as.numeric(x[2])}))
  }
}

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


####################################### Linear models
## Estimating Nested linear models

#lm_val<-lm(as.numeric(ms_data_day12_nonzero[1,])~ as.factor(RunDay_day12) + as.factor(RunDay_day12)/as.factor(SampleGroup_day12))
lm_val<-lm(as.numeric(ms_data_day12_nonzero[1,])~ as.factor(SampleGroup_day12) + as.factor(RunDay_day12))
p.val_runday_strain<-anova(lm_val)$'Pr(>F)'[1]
residuals_runday_strain<-sapply(linearmodel_list, function(x){as.numeric(x[1])})
plot(residuals_runday_strain~fitted(lm_val))
plot(lm_val)

lm_val<-lm(residuals_runday_strain_day12[,1]~ as.factor(RunDay_day12))
p.val_runday<-anova(lm_val)$'Pr(>F)'[1]
lm_val_strain<-lm(lm_val$residuals~ as.factor(SampleGroup_day12))
p.val_strain<-anova(lm_val_strain)$'Pr(>F)'[1]
return(list(p.val_runday,p.val_strain))


## store ANOVA from model
msum <- anova(lm_val)
## divide all Sums of Squares by the sum() of the Sums of Squares
msum[["Sum Sq"]]/sum(msum[["Sum Sq"]])


##### Nested Linear models

#day12
#RawData
lm_day12_nested_results<-compute_nested_linearModel(t(ScaleData(ms_data_day12_nonzero)),SampleGroup_day12,RunDay_day12)
#residuals_runday_strain_day12<-as.data.frame(lm_day12_nested_results,"r2"))#extracts the first item in the list
lm_day12_nested_results_pval<-as.data.frame(lm_day12_nested_results)
rownames(lm_day12_nested_results_pval)<-c("RunDay","Strain")
pdf("lm_day12.pdf",height=8,width=8)
par(mfrow=c(2,2))
hist(as.numeric(lm_day12_nested_results_pval[1,]),main="Runday",ylab="No. of features", xlab=paste0("Raw p value","\n","P<0.05= ",sum(as.numeric(lm_day12_nested_results_pval[1,])<0.05)))
plot(as.numeric(lm_day12_nested_results_pval[2,]),as.numeric(lm_day12_nested_results_pval[1,]),main="Nested linear model",xlab="pvalue Runday/Strain", ylab="p value Runday"); rug(as.numeric(lm_day12_nested_results_pval[2,]))
hist(as.numeric(lm_day12_nested_results_pval[2,]),main="Runday/Strain",ylab="No. of features", xlab=paste0("Raw p value","\n","P<0.05= ",sum(as.numeric(lm_day12_nested_results_pval[2,])<0.05)))
plot(log10(as.numeric(lm_day12_nested_results_pval[2,])),log10(as.numeric(lm_day12_nested_results_pval[1,])),main="Nested linear model(log)",xlab="pvalue Runday/Strain", ylab="p value Runday"); rug(as.numeric(lm_day12_nested_results_pval[2,]))
dev.off()



#day4
#RawData
lm_day4_nested_results<-compute_nested_linearModel(t(ScaleData(ms_data_day4_nonzero)),SampleGroup_day4,RunDay_day4)
#residuals_runday_strain_day4<-as.data.frame(lm_day4_nested_results,"r2"))#extracts the first item in the list
lm_day4_nested_results_pval<-as.data.frame(lm_day4_nested_results)
rownames(lm_day4_nested_results_pval)<-c("RunDay","Strain")
pdf("lm_day4.pdf",height=8,width=8)
par(mfrow=c(2,2))
hist(as.numeric(lm_day4_nested_results_pval[1,]),main="Runday",ylab="No. of features", xlab=paste0("Raw p value","\n","P<0.05= ",sum(as.numeric(lm_day4_nested_results_pval[1,])<0.05)))
plot(as.numeric(lm_day4_nested_results_pval[2,]),as.numeric(lm_day4_nested_results_pval[1,]),main="Nested linear model",xlab="pvalue Runday/Strain", ylab="p value Runday"); rug(as.numeric(lm_day4_nested_results_pval[2,]))
hist(as.numeric(lm_day4_nested_results_pval[2,]),main="Runday/Strain",ylab="No. of features", xlab=paste0("Raw p value","\n","P<0.05= ",sum(as.numeric(lm_day4_nested_results_pval[2,])<0.05)))
plot(log10(as.numeric(lm_day4_nested_results_pval[2,])),log10(as.numeric(lm_day4_nested_results_pval[1,])),main="Nested linear model(log)",xlab="pvalue Runday/Strain", ylab="p value Runday"); rug(as.numeric(lm_day4_nested_results_pval[2,]))
dev.off()

pdf("lm_day4_pval.pdf",height=8,width=8)
par(mfrow=c(2,1))
plot(as.numeric(lm_day4_nested_results_pval[2,]),as.numeric(lm_day4_nested_results_pval[1,]),main="Nested linear model",xlab="pvalue Runday/Strain", ylab="p value Runday"); rug(as.numeric(lm_day4_nested_results_pval[2,]))
dev.off()
# residuals_runday_strain_day12<-as.data.frame(lm_day12_nested_results)
# colnames(residuals_runday_strain_day12)<-rownames(ms_data_day12_nonzero)
# rownames(residuals_runday_strain_day12)<-colnames(ms_data_day12_nonzero)

# strain run day linear model
lm_day12_linear<-compute_linearModel_strain_batch(t(ScaleData(ms_data_day12_nonzero)),SampleGroup_day12,RunDay_day12)
residuals_linear_runday_strain_day12<-as.data.frame(lm_day12_linear)
colnames(residuals_linear_runday_strain_day12)<-rownames(ms_data_day12_nonzero)
rownames(residuals_linear_runday_strain_day12)<-colnames(ms_data_day12_nonzero)

##### Linear models
#day4
#RawData
lm_day12_raw<-compute_linearModel_batchEffect(t(ScaleData(ms_data_day12_nonzero)),SampleGroup_day12,RunDay_day12)
lm_runday_d12_pval_raw<-compute.r2.pval(lm_day12_raw,"r2")#extracts the first item in the list
lm_strain_d12_pval_raw<-compute.r2.pval(lm_day12_raw,"pval")

#BatchCorrected
lm_day12_corrected<-compute_linearModel_batchEffect(t(batch_corrected_mat_d12),SampleGroup_day12,RunDay_day12)
lm_runday_d12_pval_corrected<-compute.r2.pval(lm_day12_corrected,"r2")#extracts the first item in the list
lm_strain_d12_pval_corrected<-compute.r2.pval(lm_day12_corrected,"pval")

# Nested linear model
lm_day12_nested<-compute_nestedModel_significance(residuals_runday_strain_day12,SampleGroup_day12,RunDay_day12)
lm_runday_d12_pval_nested<-compute.r2.pval(lm_day12_nested,"r2")#extracts the first item in the list
lm_strain_d12_pval_nested<-compute.r2.pval(lm_day12_nested,"pval")

# Runday Strain linear model
lm_day12_linear<-compute_nestedModel_significance(t(ScaleData(batch_corrected_mat_d12)),SampleGroup_day12,RunDay_day12)
lm_runday_d12_pval_linear<-compute.r2.pval(lm_day12_linear,"r2")#extracts the first item in the list
lm_strain_d12_pval_linear<-compute.r2.pval(lm_day12_linear,"pval")

#pdf("lm_day12_comparison.pdf",height=8,width=8)
par(mfrow=c(2,2))
hist(lm_strain_d12_pval_raw,main="Raw data Strain")
hist(lm_runday_d12_pval_raw,main="Raw data Runday")
hist(lm_strain_d12_pval_corrected,main="After correction Strain")
hist(lm_runday_d12_pval_corrected,main="After correction Runday")
# hist(lm_strain_d12_pval_nested,main="Nested model correction Strain")
# hist(lm_runday_d12_pval_nested,main="Nested model correction Runday")

#plot(lm_strain_d4_nonzero_loadings_pval_corrected,lm_runday_d4_nonzero_loadings_pval_corrected)

######################################## DAY 12 and DAY4 combined
##### Nested Linear models

#RawData
lm_d4_d12_nested_results<-compute_nested_linearModel(t(ScaleData(ms_data_d4_d12)),SampleGroup_d4_d12,RunDay_d4_d12)
residuals_runday_strain_d4_d12<-as.data.frame(lm_d4_d12_nested_results)
colnames(residuals_runday_strain_d4_d12)<-rownames(ms_data_d4_d12)
rownames(residuals_runday_strain_d4_d12)<-colnames(ms_data_d4_d12)

lm_d4_d12_nested<-compute_nestedModel_significance(residuals_runday_strain_d4_d12,SampleGroup_d4_d12,RunDay_d4_d12)
lm_runday_d4_d12_pval_nested<-compute.r2.pval(lm_d4_d12_nested,"r2")#extracts the first item in the list
lm_strain_d4_d12_pval_nested<-compute.r2.pval(lm_d4_d12_nested,"pval")

####################################### Analysis of distance
## Estimating the relations between strains before and after batch correction procedure

############### DAY12

a<-adonis(t(batch_corrected_mat_d12)~ SampleGroup_day12, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(ms_data_day12_nonzero_scale) ~ SampleGroup_day12,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleGroup_day12  21    169158  8055.2  4.9991 0.52235  0.001 ***
#   Residuals          96    154685  1611.3         0.47765           
# Total             117    323844                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Subset data- from each batch
## A subset is used from the scaled data instead of subset the original data set and then scaline
## This is done to preserve the original variation in each datasets as scaling with a subset might have different distributions

## Raw data
d12_batch17_r<-RunDay_day12[RunDay_day12=="17"]
d12_batch17_s<-SampleGroup_day12[names(SampleGroup_day12) %in% names(d12_batch17_r)]#Getting the Sample group and strain names
d12_batch17_data<-ms_data_day12_nonzero_scale[,colnames(ms_data_day12_nonzero_scale) %in% names(d12_batch17_s)]

d12_batch17_data_aod<-adonis(t(d12_batch17_data)~ d12_batch17_s, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(d12_batch17_data) ~ d12_batch17_s, permutations = 999,      method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# d12_batch17_s 11     77775  7070.4   4.719 0.49013  0.001 ***
#   Residuals     54     80907  1498.3         0.50987           
# Total         65    158682                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Batch corrected data
d12_cor_batch17_data<-batch_corrected_mat_d12[,colnames(batch_corrected_mat_d12) %in% names(d12_batch17_s)]

d12_cor_batch17_data_aod<-adonis(t(d12_cor_batch17_data)~ d12_batch17_s, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(d12_cor_batch17_data) ~ d12_batch17_s, permutations = 999,      method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# d12_batch17_s 11     60677  5516.1  4.2697 0.46517  0.001 ***
#   Residuals     54     69762  1291.9         0.53483           
# Total         65    130439                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

############### DAY4

a<-adonis(t(ms_data_day4_nonzero_scale)~ SampleGroup_day4, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(ms_data_day4_nonzero_scale) ~ SampleGroup_day4,      permutations = 999, method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleGroup_day4  21    143357  6826.5   3.605 0.42363  0.001 ***
#   Residuals        103    195041  1893.6         0.57637           
# Total            124    338397                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Subset data- from each batch
## A subset is used from the scaled data instead of subset the original data set and then scaline
## This is done to preserve the original variation in each datasets as scaling with a subset might have different distributions

## Raw data
d4_batch23_r<-RunDay_day4[RunDay_day4=="23"]
d4_batch23_s<-SampleGroup_day4[names(SampleGroup_day4) %in% names(d4_batch23_r)]#Getting the Sample group and strain names
d4_batch23_data<-ms_data_day4_nonzero_scale[,colnames(ms_data_day4_nonzero_scale) %in% names(d4_batch23_s)]

d4_batch23_data_aod<-adonis(t(d4_batch23_data)~ d4_batch23_s, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(d4_batch23_data) ~ d4_batch23_s, permutations = 999,      method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# d4_batch23_s 13     61300  4715.4  2.4771 0.32792  0.001 ***
#   Residuals    66    125636  1903.6         0.67208           
# Total        79    186937                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

d4_cor_batch23_data_nmds <- metaMDS (d4_cor_batch23_data, distance= "euclidean")
d4_batch23_data_nmds <- metaMDS (d4_batch23_data, distance= "euclidean")

## Batch corrected data
d4_cor_batch23_data<-batch_corrected_mat_d4[,colnames(batch_corrected_mat_d4) %in% names(d4_batch23_s)]

d4_cor_batch23_data_aod<-adonis(t(d4_cor_batch23_data)~ d4_batch23_s, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(d4_cor_batch23_data) ~ d4_batch23_s, permutations = 999,      method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# d4_batch23_s 13     56657  4358.2  3.2337 0.3891  0.001 ***
#   Residuals    66     88951  1347.7         0.6109           
# Total        79    145608                 1.0000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pdf("aod_plots.pdf",height=8,width=8)
par(mfrow=c(2,2))
plot(density(d4_batch23_data_aod), xlim=c(0,5), main="DAY4 raw")
plot(density(d12_batch17_data_aod), xlim=c(0,5), main="DAY12 raw")
plot(density(d4_cor_batch23_data_aod), xlim=c(0,5), main="DAY4 corrected")
plot(density(d12_cor_batch17_data_aod), xlim=c(0,5), main="DAY12 corrected")
dev.off()

####################################### PCA

##Averaging replicates 
d4_batch23_data_reps<-t(avg.strains(d4_batch23_data,d4_batch23_s))
d4_cor_batch23_data_reps<-t(avg.strains(d4_cor_batch23_data,d4_batch23_s))
d12_batch17_data_reps<-t(avg.strains(d12_batch17_data,d12_batch17_s))
d12_cor_batch17_data_reps<-t(avg.strains(d12_cor_batch17_data,d12_batch17_s))

d12_batch17_data_reps_pca<-princomp(d12_batch17_data_reps,cor=F,scores=T)
residual_variance<-d12_batch17_data_reps_pca$sdev^2/sum(d12_batch17_data_reps_pca$sdev^2)

d4_cor_batch23_data_reps_pca<-princomp(d4_cor_batch23_data_reps,cor=F,scores=T)
residual_variance<-d4_cor_batch23_data_reps_pca$sdev^2/sum(d4_cor_batch23_data_reps_pca$sdev^2)

### PCA ggplot
set.seed(1)
# ggplot 

forPlot<-data.frame(PCaxisA = d12_batch17_data_reps_pca$loadings[,1],PCaxisB = d12_batch17_data_reps_pca$loadings[,2], Strains=colnames(d12_batch17_data_reps))
#Choosing color
colourCount = length(unique(forPlot$Strains))
getPalette = colorRampPalette(brewer.pal(length(unique(forPlot$Strains)),"Paired"))
set.seed(1) #important to set seed so that we obtain the same shapes for strains all the time
pch_types<-c(15, 16, 17, 18, 25, 8)
pch_values<-sample(pch_types, length(unique(forPlot$Strains)), replace = TRUE)

plot1<- ggplot(data=forPlot, aes(x=PCaxisA, y=PCaxisB,colour= factor(Strains), shape = factor(Strains))) + geom_point(size=4) #for samples
plot2<- plot1 +  scale_colour_manual('Strains', values=getPalette(colourCount)) + scale_shape_manual('Strains',values=pch_values)
plot3<- plot2+ theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1,size=4),axis.text.y=element_text(size=4),
                                  panel.grid.major.x = element_blank(), # to x remove gridlines
                                  panel.grid.major.y = element_blank(), # to y remove gridlines
                                  panel.border = element_blank(),  # remove top and right border
                                  panel.background = element_blank(),
                                  axis.line = element_line(color = 'black'))+ 
  xlab(paste0("PC 1 loadings","\n","Variation exp= ",round(residual_variance[1]*100,2),"%")) + 
  ylab(paste0("PC 2 loadings","\n","Variation exp= ",round(residual_variance[2]*100,2),"%")) +
  ggtitle("Strains-DAY12-raw") #

day4_raw<-plot3
day4_corrected<-plot3  
day12_raw<-plot3
day12_corrected<-plot3

pdf("PCA_D4-D12_batchCorrectionEffect.pdf",height=12,width=16)
multiplot(day4_raw, day4_corrected, day12_raw,day12_corrected,cols=2)
dev.off()

######################################## PCOA

bray_grps <- vegdist(t(d12_batch17_data_reps), "bray") #calculating distance matrix using bray curtis
d12_batch17_data_reps_pcoa<-cmdscale(bray_grps, eig=TRUE, add=TRUE, x.ret =TRUE) 
d12_batch17_data_reps_pcoa_scores<-as.data.frame(d12_batch17_data_reps_pcoa$x)
eig<-eigenvals(d12_batch17_data_reps_pcoa)
residual_variance<-cumsum(eig/sum(eig))

plot(d12_batch17_data_reps_pcoa_scores[,2]~d12_batch17_data_reps_pcoa_scores[,1])
text(d12_batch17_data_reps_pcoa_scores[,2]~d12_batch17_data_reps_pcoa_scores[,1], labels = colnames(d12_batch17_data_reps), cex=0.6, pos=4)
# 

### PCOA ggplot
set.seed(1)
# ggplot 

forPlot<-data.frame(PCaxisA = d12_batch17_data_reps_pcoa_scores[,1],PCaxisB = d12_batch17_data_reps_pcoa_scores[,2], Strains=colnames(d12_batch17_data_reps))
#Choosing color
colourCount = length(unique(forPlot$Strains))
getPalette = colorRampPalette(brewer.pal(length(unique(forPlot$Strains)),"Paired"))
set.seed(1) #important to set seed so that we obtain the same shapes for strains all the time
pch_types<-c(15, 16, 17, 18, 25, 8)
pch_values<-sample(pch_types, length(unique(forPlot$Strains)), replace = TRUE)

plot1<- ggplot(data=forPlot, aes(x=PCaxisA, y=PCaxisB,colour= factor(Strains), shape = factor(Strains))) + geom_point(size=4) #for samples
plot2<- plot1 +  scale_colour_manual('Strains', values=getPalette(colourCount)) + scale_shape_manual('Strains',values=pch_values)
plot3<- plot2+ theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1,size=4),axis.text.y=element_text(size=4),
                                  panel.grid.major.x = element_blank(), # to x remove gridlines
                                  panel.grid.major.y = element_blank(), # to y remove gridlines
                                  panel.border = element_blank(),  # remove top and right border
                                  panel.background = element_blank(),
                                  axis.line = element_line(color = 'black'))+ 
  xlab(paste0("PCOA 1 scores","\n","Variation exp= ",round(residual_variance[1]*100,2),"%")) + 
  ylab(paste0("PCOA 2 scores","\n","Variation exp= ",round(residual_variance[2]*100,2),"%")) +
  ggtitle("Strains-DAY12-raw") #

day4_raw<-plot3
day4_corrected<-plot3  
day12_raw<-plot3
day12_corrected<-plot3

pdf("PCA_D4-D12_batchCorrectionEffect.pdf",height=12,width=16)
multiplot(day4_raw, day4_corrected, day12_raw,day12_corrected,cols=2)
dev.off()

#### NN- graphs

ms_graph_data_strains<-avg.strains(d4_cor_batch23_data,d4_batch23_s)
rownames(ms_graph_data_strains)<-rownames(ms_graph_data_strains)
ms.d_str<-as.matrix(dist(ms_graph_data_strains,"euclidean"))
ms.1nn_str<-define.1nn.graph(ms.d_str)

pdf("d4_cor_batch17_nngraphs.pdf", width=8, height=8)
plot(define.1nn.graph(ms.d_str),layout=layout.fruchterman.reingold,vertex.label=V(ms.1nn_str)$name,vertex.label.cex=0.8,main="Strains-DAY4-batch corrected")
dev.off()