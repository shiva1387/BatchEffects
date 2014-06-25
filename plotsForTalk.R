##Linear model and Analysis of distance
##Author: Shiv
##Version :18-06-14

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


##### Linear models

#day4
#RawData
lm_day4_raw<-compute_linearModel_batchEffect(t(ScaleData(ms_data_day4_nonzero)),SampleGroup_day4,RunDay_day4)
lm_runday_d4_pval_raw<-compute.r2.pval(lm_day4_raw,"r2")#extracts the first item in the list
lm_strain_d4_pval_raw<-compute.r2.pval(lm_day4_raw,"pval")

#BatchCorrected
lm_day4_corrected<-compute_linearModel_batchEffect(t(batch_corrected_mat_d4),SampleGroup_day4,RunDay_day4)
lm_runday_d4_pval_corrected<-compute.r2.pval(lm_day4_corrected,"r2")#extracts the first item in the list
lm_strain_d4_pval_corrected<-compute.r2.pval(lm_day4_corrected,"pval")

lm_day4_pvalues<-cbind(lm_runday_d4_pval_raw,lm_strain_d4_pval_raw,lm_runday_d4_pval_corrected,lm_strain_d4_pval_corrected)
write.table(as.data.frame(lm_day4_pvalues),"day4_linearModel_pvalues.txt",quote=FALSE)

pdf("lm_day4.pdf",height=8,width=8)
par(mfrow=c(2,2))
hist(lm_strain_d4_pval_raw,main="Raw data Strain",ylim=c(0,14000))
hist(lm_runday_d4_pval_raw,main="Raw data Runday",ylim=c(0,14000))
hist(lm_strain_d4_pval_corrected,main="After correction Strain",ylim=c(0,14000))
hist(lm_runday_d4_pval_corrected,main="After correction Runday",ylim=c(0,14000))
#plot(lm_strain_d4_nonzero_loadings_pval_corrected,lm_runday_d4_nonzero_loadings_pval_corrected)
dev.off()


#day12
#RawData
lm_day12_raw<-compute_linearModel_batchEffect(t(ScaleData(ms_data_day12_nonzero)),SampleGroup_day12,RunDay_day12)
lm_runday_d12_pval_raw<-compute.r2.pval(lm_day12_raw,"r2")#extracts the first item in the list
lm_strain_d12_pval_raw<-compute.r2.pval(lm_day12_raw,"pval")

#BatchCorrected
lm_day12_corrected<-compute_linearModel_batchEffect(t(batch_corrected_mat_d12),SampleGroup_day12,RunDay_day12)
lm_runday_d12_pval_corrected<-compute.r2.pval(lm_day12_corrected,"r2")#extracts the first item in the list
lm_strain_d12_pval_corrected<-compute.r2.pval(lm_day12_corrected,"pval")
lm_day12_pvalues<-cbind(lm_runday_d12_pval_raw,lm_strain_d12_pval_raw,lm_runday_d12_pval_corrected,lm_strain_d12_pval_corrected)
write.table(as.data.frame(lm_day12_pvalues),"day12_linearModel_pvalues.txt",quote=FALSE)

pdf("lm_day12.pdf",height=8,width=8)
par(mfrow=c(2,2))
hist(lm_strain_d12_pval_raw,main="Raw data Strain",ylim=c(0,11000))
hist(lm_runday_d12_pval_raw,main="Raw data Runday",ylim=c(0,11000))
hist(lm_strain_d12_pval_corrected,main="After correction Strain",ylim=c(0,11000))
hist(lm_runday_d12_pval_corrected,main="After correction Runday",ylim=c(0,11000))
#plot(lm_strain_d12_nonzero_loadings_pval_corrected,lm_runday_d12_nonzero_loadings_pval_corrected)
dev.off()

### Analysis of distance

library(preprocessCore)
library(vegan)

data<-ms_data_day4_nonzero
processed_data<-normalize.quantiles(as.matrix(data),copy=TRUE)
processed_data<-scale(data,center=T,scale=T)
processed_data<-processed_data-min(processed_data)
###DAY12
a<-adonis(t(processed_data)~ RunDay_day12, method = "euclidean", perm=999)
# Call:
#   adonis(formula = t(processed_data) ~ RunDay_day12, permutations = 999,      method = "bray") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_day12   2    1.1987 0.59935  14.844 0.20519  0.001 ***
#   Residuals    115    4.6433 0.04038         0.79481           
# Total        117    5.8420                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Call:
#   adonis(formula = t(processed_data) ~ RunDay_day12, permutations = 999,      method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# RunDay_day12   2     48187   24093  10.051 0.1488  0.001 ***
#   Residuals    115    275657    2397         0.8512           
# Total        117    323844                 1.0000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###DAY4
# a<-adonis(t(processed_data)~ RunDay_day4, method = "bray", perm=999)
# Call:
#   adonis(formula = t(processed_data) ~ RunDay_day4, permutations = 999,      method = "bray") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_day4   3    1.0299 0.34329  9.4705 0.19016  0.001 ***
#   Residuals   121    4.3861 0.03625         0.80984           
# Total       124    5.4159                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Call:
#   adonis(formula = t(processed_data) ~ RunDay_day4, permutations = 999,      method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_day4   3     36079 12026.4  4.8135 0.10662  0.001 ***
#   Residuals   121    302318  2498.5         0.89338           
# Total       124    338397                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

### Residual variance
#Functions

calculatePC_ResidualVariance<-function(data_matrix) {
fit_mzbysam_princomp<-compute_pca(data_matrix,"scale")
residual_variance_mzbysam_princomp<-fit_mzbysam_princomp$sdev^2/sum(fit_mzbysam_princomp$sdev^2)
resid_var<-round(100-sum(residual_variance_mzbysam_princomp[1:i])*100,2)
residual_variance_mzbysam_princomp<-cbind(1:length(residual_variance_mzbysam_princomp),residual_variance_mzbysam_princomp)
return(residual_variance_mzbysam_princomp)
}

#Plot pc and residual variance

residual_variance_day4_mzbysam_princomp<-calculatePC_ResidualVariance(ms_data_day4_nonzero)
write.table(residual_variance_day4_mzbysam_princomp,"Residual_variance_day4.txt",quote=FALSE)
residual_variance_day12_mzbysam_princomp<-calculatePC_ResidualVariance(ms_data_day12_nonzero)
write.table(residual_variance_day12_mzbysam_princomp,"Residual_variance_day12.txt",quote=FALSE)

PC_Choice<-4 #number of pc's 7 for day4 and 4 for day12 removed
# for( i in 1:nrow(residual_variance_day12_mzbysam_princomp)) {
# ResidVar<-paste0("Residual Variance =", round(100-sum(residual_variance_day12_mzbysam_princomp[1:i,2])*100,2))
# residual_variance_day12_mzbysam_princomp[i,1]<-0
# png(paste0('day12/',i,".png"))
# plot(1:nrow(residual_variance_day12_mzbysam_princomp),residual_variance_day12_mzbysam_princomp[,2],pch=16,main=ResidVar,ylab="Percentage explained",xlab="PC's",type='n')
# points(1:nrow(residual_variance_day12_mzbysam_princomp),residual_variance_day12_mzbysam_princomp[,2],pch=16,col=ifelse(residual_variance_day12_mzbysam_princomp[,1]<1,"#00000033","black"))
# legend(x="top",inset=0,legend=c(paste0("PC's removed=", paste(seq(1:i),collapse=","))))
# dev.off()
# }

##To create a movie
#system("./c:/Program Files/ImageMagick-6.8.9-Q16/convert.exe -delay 0.5 *.png plot.mpg")


