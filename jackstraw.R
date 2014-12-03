library(jackstraw)

set.seed(1234)
## simulate data from a latent variable model: Y = BL + E
B = c(rep(1,50),rep(-1,50), rep(0,900))
L = rnorm(20)
E = matrix(rnorm(1000*20), nrow=1000)
dat = B %*% t(L) + E
dat = t(scale(t(dat), center=TRUE, scale=TRUE))

## apply the jackstraw
out = jackstraw(dat, PC=1, r=1)
a<-permutationPA(dat, B = 100, threshold = 0.05,
                 verbose = TRUE, seed = NULL)

##### For algae dataset

#functions

####

##day4
ms_data_day4_nonzero_scale_perm<-permutationPA(ms_data_day4_nonzero_scale, B = 100, threshold = 0.05,  verbose = TRUE, seed = NULL)
# > str(ms_data_day4_nonzero_scale_perm)
# List of 2
# $ r: int 0
# $ p: num [1:125] 1 1 1 1 1 1 1 1 1 1 ...
ms_data_day4_nonzero_scale_jack<- jackstraw(ms_data_day4_nonzero_scale, PC=1, r=1)

##day12
ms_data_day12_nonzero_scale_perm<-permutationPA(ms_data_day12_nonzero_scale, B = 100, threshold = 0.05,  verbose = TRUE, seed = NULL)
# > str(ms_data_day12_nonzero_scale_perm)
# List of 2
# $ r: int 0
# $ p: num [1:118] 1 1 1 1 1 1 1 1 1 1 ...
# ms_data_day12_nonzero_scale_jack<- jackstraw(ms_data_day12_nonzero_scale, PC=1, r=1)

### Batch corrected data
#day4
batch_corrected_mat_d4_perm<-permutationPA(batch_corrected_mat_d4, B = 100, threshold = 0.05, verbose = TRUE, seed = NULL)
# str(batch_corrected_mat_d4_perm)
# List of 2
# $ r: int 26
# $ p: num [1:125] 0 0 0 0 0 0 0 0 0 0 ...

batch_corrected_mat_d4_jack<-jackstraw(batch_corrected_mat_d4, PC=1, r=batch_corrected_mat_d4_perm$r)

#day12
batch_corrected_mat_d12_perm<-permutationPA(batch_corrected_mat_d12, B = 100, threshold = 0.05,
                 verbose = TRUE, seed = NULL)
# > str(batch_corrected_mat_d12_perm)
# List of 2
# $ r: int 22
# $ p: num [1:118] 0 0 0 0 0 0 0 0 0 0 ...

batch_corrected_mat_d12_jack<-jackstraw(batch_corrected_mat_d12, PC=1, r=batch_corrected_mat_d12_perm$r)

#### Plots

pdf("Jackstraw_batchCorrectedData.pdf",width=8,height=8)
par(mfrow=c(2,1))
pca_d4<-princomp(batch_corrected_mat_d4,cor=F,scores=T)
residual_variance_d4<-pca_d4$sdev^2/sum(pca_d4$sdev^2)
plot(1:length(residual_variance_d4),residual_variance_d4*100,pch=1,ylab="% variation", main="Batch corrected matrix- day4",
     xlab=paste0("\nSig PCs (Jackstraw)=",batch_corrected_mat_d4_perm$r,
                 ", Variance explained in Sig PCs=",round(sum(residual_variance_d4[1:batch_corrected_mat_d4_perm$r])*100),"%"))
abline(v=batch_corrected_mat_d4_perm$r,lty=2,col="red")

pca_d12<-princomp(batch_corrected_mat_d12,cor=F,scores=T)
residual_variance_d12<-pca_d12$sdev^2/sum(pca_d12$sdev^2)
plot(1:length(residual_variance_d12),residual_variance_d12*100,pch=1,ylab="% variation", main="Batch corrected matrix- day12",
     xlab=paste0("\nSig PCs (Jackstraw)=",batch_corrected_mat_d12_perm$r,
                 ", Variance explained in Sig PCs=",round(sum(residual_variance_d12[1:batch_corrected_mat_d12_perm$r])*100),"%"))
abline(v=batch_corrected_mat_d12_perm$r,lty=2,col="red")

dev.off()

library(ggplot2)

#Sample data
dat <- data.frame(dens = c(rnorm(100), rnorm(100, 10, 5))
                  , lines = rep(c("a", "b"), each = 100))
#Plot.
ggplot(dat, aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)