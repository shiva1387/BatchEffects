### R code from vignette source 'vignettes/sva/inst/doc/sva.Rnw'

###################################################
### code chunk number 1: sva.Rnw:5-6
###################################################
options(width=65)


###################################################
### code chunk number 2: input
###################################################
library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)


###################################################
### code chunk number 3: input
###################################################
pheno = pData(bladderEset)


###################################################
### code chunk number 4: input
###################################################
edata = exprs(bladderEset)


###################################################
### code chunk number 5: input
###################################################
mod = model.matrix(~as.factor(cancer), data=pheno)


###################################################
### code chunk number 6: input
###################################################
mod0 = model.matrix(~1,data=pheno)


###################################################
### code chunk number 7: input
###################################################
n.sv = num.sv(edata,mod,method="leek")
n.sv


###################################################
### code chunk number 8: input
###################################################
svobj = sva(edata,mod,mod0,n.sv=n.sv)


###################################################
### code chunk number 9: input
###################################################
pValues = f.pvalue(edata,mod,mod0)
qValues = p.adjust(pValues,method="BH")


###################################################
### code chunk number 10: input
###################################################
modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)

pValuesSv = f.pvalue(edata,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")


###################################################
### code chunk number 11: input
###################################################
fit = lmFit(edata,modSv)


###################################################
### code chunk number 12: input
###################################################
contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
fitContrasts = contrasts.fit(fit,contrast.matrix)


###################################################
### code chunk number 13: input
###################################################
eb = eBayes(fitContrasts)
topTableF(eb, adjust="BH")


###################################################
### code chunk number 14: input
###################################################
batch = pheno$batch


###################################################
### code chunk number 15: input
###################################################
mod = model.matrix(~as.factor(cancer), data=pheno)


###################################################
### code chunk number 16: input
###################################################
combat_edata = ComBat(dat=edata, batch=batch, mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)


###################################################
### code chunk number 17: input
###################################################
pValuesComBat = f.pvalue(combat_edata,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")


###################################################
### code chunk number 18: input
###################################################
modBatch = model.matrix(~as.factor(cancer) + as.factor(batch),data=pheno)
mod0Batch = model.matrix(~as.factor(batch),data=pheno)
pValuesBatch = f.pvalue(edata,modBatch,mod0Batch)
qValuesBatch = p.adjust(pValuesBatch,method="BH")


###################################################
### code chunk number 19: input
###################################################
n.sv = num.sv(edata,mod,vfilter=2000,method="leek")
svobj = sva(edata,mod,mod0,n.sv=n.sv,vfilter=2000)


###################################################
### code chunk number 20: input
###################################################
set.seed(12354)
trainIndicator = sample(1:57,size=30,replace=F)
testIndicator = (1:57)[-trainIndicator]

trainData = edata[,trainIndicator]
testData = edata[,testIndicator]

trainPheno = pheno[trainIndicator,]
testPheno = pheno[testIndicator,]


###################################################
### code chunk number 21: input
###################################################
mydata = list(x=trainData,y=trainPheno$cancer)
mytrain = pamr.train(mydata)
table(pamr.predict(mytrain,testData,threshold=2),testPheno$cancer)


###################################################
### code chunk number 22: input
###################################################
trainMod = model.matrix(~cancer,data=trainPheno)
trainMod0 = model.matrix(~1,data=trainPheno)
trainSv = sva(trainData,trainMod,trainMod0)


###################################################
### code chunk number 23: input
###################################################
fsvaobj = fsva(trainData,trainMod,trainSv,testData)
mydataSv = list(x=fsvaobj$db,y=trainPheno$cancer)
mytrainSv = pamr.train(mydataSv)
table(pamr.predict(mytrainSv,fsvaobj$new,threshold=1),testPheno$cancer)
