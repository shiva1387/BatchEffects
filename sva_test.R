##Analyzing batch effects in algal metabolomics data using SVA
##Author: Shiv
## Modified on 12 Sept 2014

## Packages

library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)

# > sessionInfo()
# R version 3.1.0 (2014-04-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] splines   parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] limma_3.20.9        pamr_1.55           survival_2.37-7     cluster_1.15.2      bladderbatch_1.2.0  Biobase_2.24.0      BiocGenerics_0.10.0 sva_3.10.0         
# [9] mgcv_1.8-3          nlme_3.1-117        corpcor_1.6.6      
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_1.2-4 digest_0.6.4     ggplot2_1.0.0    grid_3.1.0       gtable_0.1.2     lattice_0.20-29  MASS_7.3-31      Matrix_1.1-4     munsell_0.4.2    plyr_1.8.1      
# [11] proto_0.3-10     Rcpp_0.11.2      reshape2_1.4     scales_0.2.4     stringr_0.6.2    tools_3.1.0    

### Test example

## SVA example (Vignette SVA package)
pheno = pData(bladderEset)
edata = exprs(bladderEset)
mod = model.matrix(~as.factor(cancer), data=pheno)
mod0 = model.matrix(~1,data=pheno)
n.sv = num.sv(edata,mod,method="leek")
n.sv
svobj = sva(edata,mod,mod0,n.sv=n.sv)
#combat
combat_edata = ComBat(dat=edata, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)

########################################### Analyzing algal dataset

#Defining parameters
batch1=as.numeric(RunDay_day12)
batch1=gsub("15","1",batch1)
batch1=gsub("17","2",batch1)
batch1=gsub("21","3",batch1)
batch1=gsub("23","4",batch1)
batch1<-as.numeric(batch1)

## Sample description matrix 

batch_sample_desc_matrix<-data.frame(SampleNumber=seq(1:ncol(ms_data_day12_nonzero)),Batch=batch1,Strain=SampleGroup_day12)
rownames(batch_sample_desc_matrix)<-colnames(ms_data_day12_nonzero)

## Each strain (with its replicates) was analyzed only once.That is there are no continuous profiles of metabolite levels of a strain in all batches.
## Each batch analyzed a unique set of strains. For example strain A was profiles only in batch 1 and not in any other batches. 
# head(batch_sample_desc_matrix)
# SampleNumber Batch  Strain
# D12_001b1_r001            1     1 D12_001
# D12_001b1_r002            2     1 D12_001
# D12_001b2_r001            3     1 D12_001
# D12_001b2_r002            4     1 D12_001
# D12_001b3_r002            5     1 D12_001
# D12_006b1_r001            6     2 D12_006

## Parameters for SVA and Combat analysis

mod0_new = model.matrix(~1,,data=batch_sample_desc_matrix)
#mod0_new = model.matrix(~as.factor(Batch),,data=batch_sample_desc_matrix)
mod_new = model.matrix(~as.factor(Strain),data=batch_sample_desc_matrix)
edata_new = as.matrix(ms_data_day12_nonzero)


### FOR combat
combat_edata1 = ComBat(dat=edata_new, batch=batch_sample_desc_matrix$Batch, mod=mod_new, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
# Found 3 batches
# Found 21  categorical covariate(s)
# Standardizing Data across genes
# Error in solve.default(t(design) %*% design) : 
#   Lapack routine dgesv: system is exactly singular: U[24,24] = 0
# Error during wrapup: cannot open the connection

##### Possible causes
#https://stat.ethz.ch/pipermail/bioconductor/2013-August/054456.html
#http://permalink.gmane.org/gmane.science.biology.informatics.conductor/49933
#https://stat.ethz.ch/pipermail/bioconductor/2014-May/059597.html

### For SVA
# n.sv_new = num.sv(edata_new,mod_new,vfilter=2000,method="leek")
# svobj1 = sva(edata_new,mod_new,mod0_new,vfilter=2000, method="two-step")
#[1] "No significant surrogate variables"
pValues=f.pvalue(edata_new,mod_new,mod0_new)
qValues=p.adjust(pValues,method="BH")
modSv = cbind(mod_new,svobj1$sv)
mod0Sv = cbind(mod0,svobj1$sv)
pValuesSv = f.pvalue(edata_new,modSv,mod0Sv)
# Error in solve.default(t(mod) %*% mod) : 
# Lapack routine dgesv: system is exactly singular: U[23,23] = 0
# Error during wrapup: cannot open the connection
qValuesSv = p.adjust(pValuesSv,method="BH")

### For removing known batch effects with a linear model

modBatch = model.matrix(~as.factor(Strain) + as.factor(Batch), data=batch_sample_desc_matrix)
mod0Batch = model.matrix(~as.factor(Batch), data=batch_sample_desc_matrix)
pValuesBatch = f.pvalue(edata_new,modBatch,mod0Batch)
# Error in solve.default(t(mod) %*% mod) : 
# system is computationally singular: reciprocal condition number = 7.19153e-19
# Error during wrapup: cannot open the connection
qValuesBatch = p.adjust(pValuesBatch,method="BH")


