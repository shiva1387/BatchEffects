
#d4.norm
# > head(sort(lm_pca_scores_runday4_nonzero_loadings_pval))
# Comp.2       Comp.1       Comp.8       Comp.6      Comp.11       Comp.5 
# 3.007729e-43 3.378846e-22 4.523162e-09 2.494628e-08 2.296277e-07 7.847854e-06 
# > head(sort(lm_pca_scores_runday4_nonzero_loadings_r.sq,decreasing=TRUE))
# Comp.2    Comp.1    Comp.8    Comp.6   Comp.11    Comp.5 
# 0.8110657 0.5750203 0.2928821 0.2720464 0.2439658 0.1967930 

#d4.scale
# > head(sort(lm_pca_scores_runday4_nonzero_loadings_pval))
# Comp.2       Comp.1       Comp.8       Comp.6      Comp.11       Comp.5 
# 2.484598e-43 4.553474e-21 3.490789e-09 1.489811e-08 3.768702e-07 6.002051e-06 
# > head(sort(lm_pca_scores_runday4_nonzero_loadings_r.sq,decreasing=TRUE))
# Comp.2    Comp.1    Comp.8    Comp.6   Comp.11    Comp.5 
# 0.8116675 0.5560723 0.2959875 0.2784050 0.2375382 0.2004892 


######## SVD
#for day 4

svd_day4_nonzero <- svd(d4.scale)
str(svd_day4_nonzero)

svd_day4_nonzero$d1<-svd_day4_nonzero$d
svd_day4_nonzero$d1[c(2,1,8,6,11)]<-0 #Removing the variation caused by runday(id using pca) where the multiple r2 correlation is above 0.5
d4_scale_rmbatch_v1<-svd_day4_nonzero$u %*% diag(svd_day4_nonzero$d1) %*% t(svd_day4_nonzero$v)

rownames(d4_scale_rmbatch_v1)<-rownames(ms_data_day4_nonzero)
colnames(d4_scale_rmbatch_v1)<-names(ms_data_day4_nonzero)

#u1 %*% d1 %*% t(v1)   #ms.svd$u%*%diag(ms.svd$d) #http://www.ats.ucla.edu/stat/r/pages/svd_demos.htm
###################
#### Against Strain
###################

##corrected
day4_sig_features_v1<-mt.maxT(d4_scale_rmbatch_v1,classlabel_samplegrp_d4,test="f",side="abs",fixed.seed.sampling="y",B=100000,nonpara="n") #using multtest package
id.sig_day4_v1 <- which(day4_sig_features_v1$adjp < 0.05 );
metab.sig_day4_v1<-cbind(d4.scale[id.sig_day4_v1,],round(day4_sig_features_v1$adjp[id.sig_day4_v1],5))

###################
#### Against runday
###################

##corrected
day4_sig_features_r_v1<-mt.maxT(d4_scale_rmbatch_v1,classlabel_runday_d4,test="f",side="abs",fixed.seed.sampling="y",B=100000,nonpara="n") #using multtest package
id.sig_day4_r_v1 <- which(day4_sig_features_r_v1$adjp < 0.05 );
metab.sig_day4_r_v1<-cbind(d4.scale[id.sig_day4_r_v1,],round(day4_sig_features_r_v1$adjp[id.sig_day4_r_v1],5))

#####################################################################################################

#d12.norm
# > head(sort(lm_pca_scores_runday12_nonzero_loadings_pval))
# Comp.3       Comp.1       Comp.4       Comp.2      Comp.12      Comp.10 
# 1.406748e-30 1.098893e-19 2.908428e-18 1.277628e-14 1.242495e-02 2.728178e-02
# > head(sort(lm_pca_scores_runday12_nonzero_loadings_r.sq,decreasing=TRUE))
# Comp.3     Comp.1     Comp.4     Comp.2    Comp.12    Comp.10 
# 0.68808406 0.52284414 0.49560140 0.41854605 0.07167529 0.05921720 
#d12.scale
# > head(sort(lm_pca_scores_runday12_nonzero_loadings_pval))
# Comp.3       Comp.4       Comp.1       Comp.2      Comp.12      Comp.10 
# 4.187723e-30 1.694542e-18 1.891703e-18 9.571292e-15 7.223229e-03 3.225789e-02 
# > head(sort(lm_pca_scores_runday12_nonzero_loadings_r.sq,decreasing=TRUE))
# Comp.3     Comp.4     Comp.1     Comp.2    Comp.12    Comp.10 
# 0.68226327 0.50019857 0.49926531 0.42138548 0.08017055 0.05654185 

######## SVD
#for day 12

svd_day12_nonzero <- svd(d12.scale)
str(svd_day12_nonzero)

svd_day12_nonzero$d1<-svd_day12_nonzero$d
svd_day12_nonzero$d1[c(3,4,1,2,12)]<-0 #Removing the variation caused by runday(id using pca) where the multiple r2 correlation is above 0.5
d12_scale_rmbatch_v1<-svd_day12_nonzero$u %*% diag(svd_day12_nonzero$d1) %*% t(svd_day12_nonzero$v)

rownames(d12_scale_rmbatch_v1)<-rownames(ms_data_day12_nonzero)
colnames(d12_scale_rmbatch_v1)<-names(ms_data_day12_nonzero)
#u1 %*% d1 %*% t(v1)   #ms.svd$u%*%diag(ms.svd$d) #http://www.ats.ucla.edu/stat/r/pages/svd_demos.htm

###################
#### Against Strain
###################

##corrected
day12_sig_features_v1<-mt.maxT(d12_scale_rmbatch_v1,classlabel_samplegrp_d12,test="f",side="abs",fixed.seed.sampling="y",B=100000,nonpara="n") #using multtest package
id.sig_day12_v1 <- which(day12_sig_features_v1$adjp < 0.05 );
metab.sig_day12_v1<-cbind(d12.scale[id.sig_day12_v1,],round(day12_sig_features_v1$adjp[id.sig_day12_v1],5))

###################
#### Against runday
###################

##corrected
day12_sig_features_r_v1<-mt.maxT(d12_scale_rmbatch_v1,classlabel_runday_d12,test="f",side="abs",fixed.seed.sampling="y",B=100000,nonpara="n") #using multtest package
id.sig_day12_r_v1 <- which(day12_sig_features_r_v1$adjp < 0.05 );
metab.sig_day12_r_v1<-cbind(d12.scale[id.sig_day12_r_v1,],round(day12_sig_features_r_v1$adjp[id.sig_day12_r_v1],5))

#plotting results

#pval_bins<-c(0,0.001,0.01,0.05,1)

png('pval_day12.png',width=1200,height=1200)
par(mfrow=c(2,2))
hist(day12_sig_features_nc,breaks=20,main="Strain",freq=TRUE,xlab=paste("no of. sig feat(0.05)=",length(id.sig_day12_nc)))
hist(day12_sig_features_nc_r,breaks=20,main="Runday",freq=TRUE,xlab=paste("no of. sig feat(0.05)=",length(id.sig_day12_nc_r)))
hist(day12_sig_features,breaks=20,main="Strain-corrected",freq=TRUE,xlab=paste("no of. sig feat(0.05)=",length(id.sig_day12)))
hist(day12_sig_features_r,breaks=20,main="Runday-corrected",freq=TRUE,xlab=paste("no of. sig feat(0.05)=",length(id.sig_day12_r)))
dev.off()

mz_rt_day12_nc <- strsplit(rownames(metab.sig_day12_nc), "\\@")
mz_day12_nc<-sapply(mz_rt_day12_nc , function (x) if(length(x) == 2) x[1] else as.character(NA))
rt_day12_nc<-sapply(mz_rt_day12_nc , function (x) if(length(x) == 2) x[2] else as.character(NA))

mz_rt_day12 <- strsplit(rownames(metab.sig_day12), "\\@")
mz_day12<-sapply(mz_rt_day12 , function (x) if(length(x) == 2) x[1] else as.character(NA))
rt_day12<-sapply(mz_rt_day12 , function (x) if(length(x) == 2) x[2] else as.character(NA))

png('d12_mzrt_sig_features_strain_kruskal.png',height=800,width=800)
plot(rt,mz,pch=1,main="sig featuresVsStrain(0.05)",ylab="m/z",xlab="rt",col="#00000033")
points(rt_day12_nc,mz_day12_nc,pch=1,col="green")
points(rt_day12,mz_day12,pch=1,col="red")
legend("topright", title="Colors", c("full data","non-cor","corr"), fill=c("#00000033","green","red"), cex=0.5)
dev.off()

mz_rt_day12_nc_r <- strsplit(rownames(metab.sig_day12_nc_r), "\\@")
mz_day12_nc_r<-sapply(mz_rt_day12_nc_r , function (x) if(length(x) == 2) x[1] else as.character(NA))
rt_day12_nc_r<-sapply(mz_rt_day12_nc_r , function (x) if(length(x) == 2) x[2] else as.character(NA))

mz_rt_day12_r <- strsplit(rownames(metab.sig_day12_r), "\\@")
mz_day12_r<-sapply(mz_rt_day12_r , function (x) if(length(x) == 2) x[1] else as.character(NA))
rt_day12_r<-sapply(mz_rt_day12_r , function (x) if(length(x) == 2) x[2] else as.character(NA))

png('d12_mzrt_sig_features_runday_kruskal.png',height=600,width=600)
par(mfrow=c(2,1))
plot(rt,mz,pch=1,main="d12 sig featuresVsRunday-nc",ylab="m/z",xlab="rt",col="#00000033")
points(rt_day12_nc_r,mz_day12_nc_r,pch=1,col="green")
plot(rt,mz,pch=1,main="d12 sig featuresVsRunday-Corrected",ylab="m/z",xlab="rt",col="#00000033")
points(rt_day12_r,mz_day12_r,pch=1,col="green")
dev.off()

#matplot(cbind(fit_day12_mzbysam_princomp$sdev,svd_day12_nonzero$d))

#################### Plots for batch effect correction ####################

pdf("LinearModel_PCStrainRunday.pdf",width=16)
m<- matrix(c(1,2,3,4,5,5),ncol = 2,byrow = TRUE,heights = c(0.5,0.5,0.1))
layout(mat = m)

plot(lm_pca_scores_strain_day4_nonzero_loadings_r.sq,type="l", xaxt="n",xlab="PC's",ylab="mulriple r2",main="day4")
points(lm_pca_scores_runday4_nonzero_loadings_r.sq,type="l",lty=2,col="red")
axis(1, at = seq(0, 125, by = 5))

plot(lm_pca_scores_strain_day12_nonzero_loadings_r.sq,type="l",xaxt="n",xlab="PC's",ylab="mulriple r2",main="day12")
points(lm_pca_scores_runday12_nonzero_loadings_r.sq,type="l",lty=2,col="red")
axis(1, at = seq(0, 125, by = 5))

plot(lm_pca_scores_strain_day4_nonzero_loadings_pval,type="l",xaxt="n",xlab="PC's",ylab="pval",main="day4")
points(lm_pca_scores_runday4_nonzero_loadings_pval,type="l",lty=2,col="red")
axis(1, at = seq(0, 125, by = 5))

plot(lm_pca_scores_strain_day12_nonzero_loadings_pval,type="l",xaxt="n",xlab="PC's",ylab="pval",main="day12")
points(lm_pca_scores_runday12_nonzero_loadings_pval,type="l",lty=2,col="red")
axis(1, at = seq(0, 125, by = 5))

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("Strains", "RunDay"),  col=c("black","red"), lty=c(1,2),horiz = TRUE)

dev.off()

pdf("BatchCorrectedPermuteTest_FeatureVsStrainRunday.pdf",width=12,height=12)
m<- matrix(c(1,2,3,4,5,5),ncol = 2,byrow = TRUE)
layout(mat = m)

plot(1:length(day4_sig_features_nc$adjp),day4_sig_features_nc$adjp,main="day4 strains",xlab="features",ylab="adj. pval")
points(1:length(day4_sig_features$adjp),day4_sig_features$adjp,col="red")
points(1:length(day4_sig_features_v1$adjp),day4_sig_features_v1$adjp,col="blue")

plot(1:length(day4_sig_features_nc_r$adjp),day4_sig_features_nc_r$adjp,main="day4 runday",xlab="features",ylab="adj. pval")
points(1:length(day4_sig_features_r$adjp),day4_sig_features_r$adjp,col="red")
points(1:length(day4_sig_features_r_v1$adjp),day4_sig_features_r_v1$adjp,col="blue")

plot(1:length(day12_sig_features_nc$adjp),day12_sig_features_nc$adjp,main="day12 strains",xlab="features",ylab="adj. pval")
points(1:length(day12_sig_features$adjp),day12_sig_features$adjp,col="red")
points(1:length(day12_sig_features_v1$adjp),day12_sig_features_v1$adjp,col="blue")

plot(1:length(day12_sig_features_nc_r$adjp),day12_sig_features_nc_r$adjp,main="day12 runday",xlab="features",ylab="adj. pval")
points(1:length(day12_sig_features_r$adjp),day12_sig_features_r$adjp,col="red")
points(1:length(day12_sig_features_r_v1$adjp),day12_sig_features_r_v1$adjp,col="blue")

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend(x = "top",inset = 0,
       legend = c("Non-corrected", "PC(R2<0.5)","Top 5 PC'S"),   text.col=c("black","red","blue"), horiz = TRUE)
dev.off()
