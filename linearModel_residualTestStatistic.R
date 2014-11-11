## Analyzing linear model--Extracting and generating test statistics on the residuals from a model built on runday
## Shiv 10 Nov 2014

### Libraries

library(Hmisc)

### Functions

compute_linearModel_ext_resid<-function(data_matrix,StrainId,RunDayId) { #dependent.factor1 is Strain id(sample groups) and dependent.factor2 is RunDay 
  lm_pca_scores<-apply(data_matrix,2, function(x) {
    lm_val<-lm(x~ as.factor(RunDayId))
    #p.val_runday<-anova(lm_val)$'Pr(>F)'[1]
    lm_val_strain<-lm(lm_val$residuals~ as.factor(StrainId))
    #p.val_strain<-anova(lm_val_strain)$'Pr(>F)'[1]
    return(lm_val$residuals)
  })
}


#######################


lm_day4_raw<-compute_linearModel_ext_resid(t(ScaleData(ms_data_day4_nonzero)),SampleGroup_day4,RunDay_day4)
data_matrix_d4<-as.matrix(t(lm_day4_raw))
classlabel_factor_d4<-as.numeric(as.factor(SampleGroup_day4))-1
dataset_sig_features_day4<-mt.maxT(data_matrix_d4,classlabel_factor_d4,test="f",side="abs",fixed.seed.sampling="y",B=1000,nonpara="n")
# head(dataset_sig_features_day4)
# index  teststat  rawp  adjp
# 504.288785798976@458.052   8314 123.06919 0.001 0.001
# 325.235718390465@412.964   5012 115.36232 0.001 0.001
# 504.304757321015@448.017   8318 112.54592 0.001 0.001
# 395.242420519267@446.756   6412  97.02215 0.001 0.001
# 924.566737966965@697.4745 12985  88.85149 0.001 0.001
# 925.569988480991@697.451  13004  88.58353 0.001 0.001
metab_sig_day4<-cbind(data_matrix_d4[id.sig_dataset_day4,],round(dataset_sig_features_day4$adjp[dataset_sig_features_day4$index %in% id.sig_dataset_day4],5))
id.sig_dataset_day4 <- sort(dataset_sig_features_day4[dataset_sig_features_day4$adjp < 0.05,c(1)]) #getting the column which provides index of rows satisying the condition
length(id.sig_dataset_day4)

lm_day12_raw<-compute_linearModel_ext_resid(t(ScaleData(ms_data_day12_nonzero)),SampleGroup_day12,RunDay_day12)
data_matrix_d12<-as.matrix(t(lm_day12_raw))
classlabel_factor_d12<-as.numeric(as.factor(SampleGroup_day12))-1
dataset_sig_features_day12<-mt.maxT(data_matrix_d12,classlabel_factor_d12,test="f",side="abs",fixed.seed.sampling="y",B=1000,nonpara="n")
id.sig_dataset_day12 <- sort(dataset_sig_features_day12[dataset_sig_features_day12$adjp < 0.05,c(1)]) #getting the column which provides index of rows satisying the condition
metab_sig_day12<-cbind(data_matrix_d12[id.sig_dataset_day12,],round(dataset_sig_features_day12$adjp[dataset_sig_features_day12$index %in% id.sig_dataset_day12],5))
length(id.sig_dataset_day12)

# adjusted p values from matrix
#for day4 3208
#for day12 2555

############ Num and denum of test statistic

lm_day4_raw_numdenum<-mt.teststat.num.denum(as.matrix(t(lm_day4_raw)),classlabel_factor_d4,test="f",nonpara="n")
lm_day4_raw_num_values_s<-lm_day4_raw_numdenum[,1]
lm_day4_raw_denum_values_s<-lm_day4_raw_numdenum[,2]

lm_day12_raw_numdenum<-mt.teststat.num.denum(as.matrix(t(lm_day12_raw)),classlabel_factor_d12,test="f",nonpara="n")
lm_day12_raw_num_values_s<-lm_day12_raw_numdenum[,1]
lm_day12_raw_denum_values_s<-lm_day12_raw_numdenum[,2]

### Correlation as a function of test statistic (Rodgers et al, 1988) # Coded by Shiv March 7, 2014
#--function to calculate the r2 value associated with an F-statistic
#--F is an F-statistic, k is the number of groups, N is the total number of samples
r2.from.Fstat<-function(F,k,N){(F*(k-1))/(F*(k-1)+(N-k))}

# compute fvalue and r2
dataset_fvalue_d4<-mt.teststat(data_matrix_d4,classlabel_factor_d4,test="f",nonpara="n")
dataset_r2_d4<-r2.from.Fstat(dataset_fvalue_d4,length(unique(classlabel_factor_d4)),ncol(data_matrix_d4))

dataset_fvalue_d12<-mt.teststat(data_matrix_d12,classlabel_factor_d12,test="f",nonpara="n")
dataset_r2_d12<-r2.from.Fstat(dataset_fvalue_d12,length(unique(classlabel_factor_d12)),ncol(data_matrix_d12))

mz_rt_day4 <- strsplit(rownames(ms_data_day4_nonzero), "\\@")
mz_day4<-sapply(mz_rt_day4 , function (x) if(length(x) == 2) x[1] else as.character(NA))
rt_day4<-sapply(mz_rt_day4 , function (x) if(length(x) == 2) x[2] else as.character(NA))

mz_rt_day12 <- strsplit(rownames(ms_data_day12_nonzero), "\\@")
mz_day12<-sapply(mz_rt_day12 , function (x) if(length(x) == 2) x[1] else as.character(NA))
rt_day12<-sapply(mz_rt_day12 , function (x) if(length(x) == 2) x[2] else as.character(NA))

sigfeat_day4_bc_s<-strsplit(rownames(metab_sig_day4), "\\@")
sigfeat_day4_bc_mz_s<-sapply(sigfeat_day4_bc_s , function (x) if(length(x) == 2) x[1] else as.character(NA))
sigfeat_day4_bc_rt_s<-sapply(sigfeat_day4_bc_s , function (x) if(length(x) == 2) x[2] else as.character(NA))
tot_mz_s_day4<-length(sigfeat_day4_bc_mz_s)

sigfeat_day12_bc_s<-strsplit(rownames(metab_sig_day12), "\\@")
sigfeat_day12_bc_mz_s<-sapply(sigfeat_day12_bc_s , function (x) if(length(x) == 2) x[1] else as.character(NA))
sigfeat_day12_bc_rt_s<-sapply(sigfeat_day12_bc_s , function (x) if(length(x) == 2) x[2] else as.character(NA))
tot_mz_s_day12<-length(sigfeat_day12_bc_mz_s)

day4_numdenum_s<-lm_day4_raw_denum_values_s/lm_day4_raw_num_values_s
day4_numdenum_s_ind<-match(sigfeat_day4_bc_s,mz_rt_day4); day4_numdenum_s_ind<-day4_numdenum_s_ind[!is.na(day4_numdenum_s_ind)]

day12_numdenum_s<-lm_day12_raw_denum_values_s/lm_day12_raw_num_values_s
day12_numdenum_s_ind<-match(sigfeat_day12_bc_s,mz_rt_day12); day12_numdenum_s_ind<-day12_numdenum_s_ind[!is.na(day12_numdenum_s_ind)]

tmp1<-day4_numdenum_s[day4_numdenum_s_ind]
tmp1<-tmp1[tmp1<0.1]

tmp2<-day12_numdenum_s[day12_numdenum_s_ind]
tmp2<-tmp2[tmp2<0.1]

#plotting test statistic vs no-of-sig-features strain
sigfeat_day4_s<-rep('nosig',length(rt_day4))
sigfeat_day4_s[day4_numdenum_s_ind]<-'sig'

sigfeat_day12_s<-rep('nosig',length(rt_day12))
sigfeat_day12_s[day12_numdenum_s_ind]<-'sig'

df1_d4s<-data.frame(lm_day4_raw_denum_values_s,sigfeat_day4_s)
colnames(df1_d4s)<-c('denom','sigfeat')
df1_d4s$bins <- cut2(df1_d4s$denom, c(0.025,0.05,0.1,0.5,1,1.5,2,3))
df2_d4s<-as.data.frame.matrix(table(df1_d4s$bins,df1_d4s$sigfeat))
#colnames(df2_d4s)<-c('bin','nosig','sig')
df2_d4s$bin<-gsub('\\[|\\)|\\]', '', rownames(df2_d4s))
splits<-strsplit(df2_d4s$bin, ",")
df2_d4s$newBin<-as.numeric(sapply(splits , function (x) if(length(x) == 2) x[2] else as.character(NA)))

df1_d12s<-data.frame(lm_day12_raw_denum_values_s,sigfeat_day12_s)
colnames(df1_d12s)<-c('denom','sigfeat')
df1_d12s$bins <- cut2(df1_d12s$denom, c(0.025,0.05,0.1,0.5,1,1.5,2,3))
df2_d12s<-as.data.frame.matrix(table(df1_d12s$bins,df1_d12s$sigfeat))
#colnames(df2_d12s)<-c('bin','nosig','sig')
df2_d12s$bin<-gsub('\\[|\\)|\\]', '', rownames(df2_d12s))
splits<-strsplit(df2_d12s$bin, ",")
df2_d12s$newBin<-as.numeric(sapply(splits , function (x) if(length(x) == 2) x[2] else as.character(NA)))

#plot

pdf("linearModel_testStatistic.pdf",height=12,width=12)

par(mfrow=c(4,2))

plot(rt_day4,mz_day4,pch=1,main="day 4 Strain",ylab="m/z",xlab=paste0("rt(s)","\n","no of features=",tot_mz_s_day4),col="#00000033")
points(sigfeat_day4_bc_rt_s,sigfeat_day4_bc_mz_s,pch=1,col="red")

plot(rt_day12,mz_day12,pch=1,main="day 12 Strain",ylab="m/z",xlab=paste0("rt(s)","\n","no of features=",tot_mz_s_day12),col="#00000033")
points(sigfeat_day12_bc_rt_s,sigfeat_day12_bc_mz_s,pch=1,col="red")

plot(rt_day4,mz_day4,pch=1,main="day 4 Strain-Test statistic",ylab="m/z",xlab=paste0("rt(s)","\n","no of features=",length(tmp1)),col="#00000033")
points(sigfeat_day4_bc_rt_s,sigfeat_day4_bc_mz_s,pch=1,col=ifelse(day4_numdenum_s[day4_numdenum_s_ind]<0.1,"red","#00000033"))

plot(rt_day12,mz_day12,pch=1,main="day 12 Strain-Test statistic",ylab="m/z",xlab=paste0("rt(s)","\n","no of features=",length(tmp2)),col="#00000033")
points(sigfeat_day12_bc_rt_s,sigfeat_day12_bc_mz_s,pch=1,col=ifelse(day12_numdenum_s[day12_numdenum_s_ind]<0.1,"red","#00000033"))

plot(dataset_fvalue_d4,dataset_r2_d4,pch=1,main="day 4 Strain-F val VS r2 statistic",ylab="r2",xlab="fvalue",col=ifelse(sigfeat_day4_s=="sig","red","#00000033"))
legend(x="topright",inset=0, legend=c("sig","non-sig"),col=c("red","black"),pch=1)  

plot(dataset_fvalue_d12,dataset_r2_d12,pch=1,main="day 12 Strain-F val VS r2 statistic",ylab="r2",xlab="fvalue",col=ifelse(sigfeat_day12_s=="sig","red","#00000033"))
legend(x="topright",inset=0, legend=c("sig","non-sig"),col=c("red","black"),pch=1)  

plot(log(df2_d4s$newBin),df2_d4s$nosig,type='l', main="Strain-test statistic variance Vs Sig features", ylim=c(0,max(df2_d4s$nosig,df2_d4s$sig)))
lines(log(df2_d4s$newBin),df2_d4s$sig,col="red")
legend(x="topright",inset=0, c("sig","non-sig"),col=c("red","black"),lty=1)    

plot(log(df2_d12s$newBin),df2_d12s$nosig,type='l', main="Strain-test statistic variance Vs Sig features", ylim=c(0,max(df2_d12s$nosig,df2_d12s$sig)))
lines(log(df2_d12s$newBin),df2_d12s$sig,col="red")
legend(x="topright",inset=0, c("sig","non-sig"),col=c("red","black"),lty=1)    

dev.off()