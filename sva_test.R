# Modified on 16 May 2014

library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)

mod1 = model.matrix(~as.factor(SampleGroup_day12))
edata1 = as.matrix(ms_data_day12_nonzero)
batch1=as.numeric(RunDay_day12)
batch1=gsub("15","1",batch1)
batch1=gsub("17","2",batch1)
batch1=gsub("21","3",batch1)
batch1=gsub("23","4",batch1)
batch1<-as.numeric(batch1)
#combat_edata = ComBat(dat=edata1, batch=batch1, mod=mod1, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)

# Found 3 batches
# Found 21  categorical covariate(s)
# Standardizing Data across genes
# Error in solve.default(t(design) %*% design) : 
#   Lapack routine dgesv: system is exactly singular: U[24,24] = 0
# Error during wrapup: cannot open the connection

#https://stat.ethz.ch/pipermail/bioconductor/2013-August/054456.html
#http://permalink.gmane.org/gmane.science.biology.informatics.conductor/49933

modBatch = model.matrix(~as.factor(SampleGroup_day12) + as.factor(batch1))
mod0Batch = model.matrix(~as.factor(batch1))
pValuesBatch = f.pvalue(edata1,modBatch,mod0Batch)
qValuesBatch = p.adjust(pValuesBatch,method="BH")

#https://stat.ethz.ch/pipermail/bioconductor/2014-May/059597.html
