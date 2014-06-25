##Identifying strain-specific metabolite associations
##Author: Shiv
##Version :20-06-14

# Rohan 13/06/14
# norm.to.total.expression<-function(x){x^2/sum(x^2)}
# energy.test<-t(apply(exp.ave.within.tissuesAverage,1,norm.to.total.expression))


# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2902135/?report=printable
# 
# Hi Shiv,
# This is the paper, the method is described in the Sup Mat (attached). Seems simple enough, you just calculate a 2-group test (e.g. t-statistic or whatever) 
# between each strain and each of the other strains, for each mass feature. 
# Strains that have a "specific" mass feature intensity will have high levels of the test statistic for all strain comparisons.

### Strain-specific association values for each metabolite-Algorithm by Rohan(2008)

norm.to.total.expression<-function(x){x^2/sum(x^2)}

association.values.d4<-apply(mappedMetabolites_mean_d4,2,norm.to.total.expression)
association.values.d12<-apply(mappedMetabolites_mean_d12,2,norm.to.total.expression)

plot(association.values.d4[,180], pch=20, xaxt="n", xlab="", ylab="Haygood measure",main=colnames(mappedMetabolites_mean_d4)[180])
# Plot the axis separately
axis(1, at=1:nrow(association.values.d4), labels=rownames(association.values.d4),las=2)
which(association.values.d4==max(association.values.d4),arr.ind=TRUE)
# row col
# D4_094   7 180
# > colnames(mappedMetabolites_mean_d4)[180]
# [1] "C09235"
#plot(mappedMetabolites_mean_d4[,180],association.values.d4[,180])

plot(1:ncol(association.values.d12),rowSums(association.values.d12), pch=20, xaxt="n", xlab="", ylab="Haygood measure-Sum",main="rowSums")
# Plot the axis separately
axis(1, at=1:nrow(association.values.d4), labels=rownames(association.values.d4),las=2)



plot(association.values.d12[,547], pch=20, xaxt="n", xlab="", ylab="Haygood measure",main=colnames(mappedMetabolites_mean_d12)[547])
# Plot the axis separately
axis(1, at=1:nrow(association.values.d12), labels=rownames(association.values.d12),las=2)
which(association.values.d12==max(association.values.d12),arr.ind=TRUE)
# row col
# D12_187  11 547
#plot(mappedMetabolites_mean_d12[,547],association.values.d12[,547])

### Plots

for( i in 1:nrow(association.values.d4)) {
png(paste0('day4/',rownames(association.values.d4)[i],".png"))
par(mfrow=c(2,1))
plot(association.values.d4[i,],as.numeric(mappedMetabolites_mean_d4[i,]),main=rownames(association.values.d4)[i])
hist(association.values.d4[i,],main=rownames(association.values.d4)[i])
dev.off()
}

## test statistic

#dataset_fvalue<-mt.teststat(data_matrix,classlabel_factor,test="f",nonpara="n")

