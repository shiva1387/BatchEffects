## Analyzing linear model--Extracting and generating test statistics on the residuals from a model built on runday
## Shiv 10 Nov 2014
## Modified on 04-Dec-2014
### Libraries

library(Hmisc)

### Functions

compute_linearModel_ext_resid<-function(data_matrix,StrainId,RunDayId) { #dependent.factor1 is Strain id(sample groups) and dependent.factor2 is RunDay 
  lm_pca_scores<-apply(data_matrix,2, function(x) {
    lm_val<-lm(x~ as.factor(RunDayId))
    p.val_runday<-anova(lm_val)$'Pr(>F)'[1]
    lm_val_strain<-lm(lm_val$residuals~ as.factor(StrainId))
    #p.val_strain<-anova(lm_val_strain)$'Pr(>F)'[1]
    return(list(pvalue=p.val_runday,resids=lm_val$residuals))
    #return(list(p.val_runday,lm_val$residuals)) #
  })
} # Modified to return both p values(run day) and residuals

extract_pval<-function(pval_residuals) {
  runday_pval<-rep(999,length(pval_residuals))
  for(i in 1:length(pval_residuals))
  {
    runday_pval[i]<-pval_residuals[[i]][1]
  }
 return(runday_pval)
}

extract_residual<-function(pval_residuals) {
  strain_residual<-vector("list", length(pval_residuals)) #list()
  for(i in 1:length(pval_residuals))
  {
    strain_residual[i]<-pval_residuals[[i]][2]
  }
  strain_residual_df<-do.call(rbind.data.frame, strain_residual)
  return(strain_residual_df)
}


#######################


lm_day4_raw<-compute_linearModel_ext_resid(t(ScaleData(ms_data_day4_nonzero)),SampleGroup_day4,RunDay_day4)
p.val_runday_d4<-as.numeric(extract_pval(lm_day4_raw))
data_matrix_d4<-extract_residual(lm_day4_raw)
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
id.sig_dataset_day4 <- sort(dataset_sig_features_day4[dataset_sig_features_day4$adjp < 0.05,c(1)]) #getting the column which provides index of rows satisying the condition
metab_sig_day4<-cbind(data_matrix_d4[id.sig_dataset_day4,],round(dataset_sig_features_day4$adjp[dataset_sig_features_day4$index %in% id.sig_dataset_day4],5))
length(id.sig_dataset_day4)

lm_day12_raw<-compute_linearModel_ext_resid(t(ScaleData(ms_data_day12_nonzero)),SampleGroup_day12,RunDay_day12)
p.val_runday_d12<-as.numeric(extract_pval(lm_day12_raw))
data_matrix_d12<-extract_residual(lm_day12_raw)
classlabel_factor_d12<-as.numeric(as.factor(SampleGroup_day12))-1
dataset_sig_features_day12<-mt.maxT(data_matrix_d12,classlabel_factor_d12,test="f",side="abs",fixed.seed.sampling="y",B=1000,nonpara="n")
id.sig_dataset_day12 <- sort(dataset_sig_features_day12[dataset_sig_features_day12$adjp < 0.05,c(1)]) #getting the column which provides index of rows satisying the condition
metab_sig_day12<-cbind(data_matrix_d12[id.sig_dataset_day12,],round(dataset_sig_features_day12$adjp[dataset_sig_features_day12$index %in% id.sig_dataset_day12],5))
length(id.sig_dataset_day12)


#### functions
calculatePOverlaps<-function(runday,strain,threshold)
{
  common_05<-intersect(which(runday<threshold),which(strain<threshold))
  ronly_05<-setdiff(which(runday<threshold),common_05)
  sonly_05<-setdiff(which(strain<threshold),common_05)
  return(cbind(length(common_05),length(ronly_05),length(sonly_05),threshold))
}

binseq<-seq(0,1,0.05)
get.hist.data<-function(x,b){hist(x,b,plot=F)$counts}

########## day4
d4.pval.runday.hd<-get.hist.data(p.val_runday_d4,binseq)
d4.pval.strain.hd<-get.hist.data(dataset_sig_features_day4$adjp,binseq)
d12.pval.runday.hd<-get.hist.data(p.val_runday_d12,binseq)
d12.pval.strain.hd<-get.hist.data(dataset_sig_features_day12$adjp,binseq)
hd.xaxis<-seq(0,1,0.05)[-21]

range12<-max(c(d4.pval.runday.hd,d4.pval.strain.hd,d12.pval.runday.hd,d12.pval.strain.hd))

pdf("linearModel_permutation_all.pdf",width=8,height=8)
par(mfrow=c(2,2))
plot(hd.xaxis,d4.pval.strain.hd,type="s",ylim=c(100,range12),log="y",las=1,ylab="Number of mass features",xlab="P-value",main="Day4 Strain")
lines(hd.xaxis,rep((0.05*13443),20),lty=2)
plot(hd.xaxis,d4.pval.runday.hd,type="s",ylim=c(100,range12),log="y",las=1,ylab="Number of mass features",xlab="P-value",main="Day4 Runday")
lines(hd.xaxis,rep((0.05*13443),20),lty=2)
plot(hd.xaxis,d12.pval.strain.hd,type="s",ylim=c(100,range12),log="y",las=1,ylab="Number of mass features",xlab="P-value",main="Day12 Strain")
lines(hd.xaxis,rep((0.05*10687),20),lty=2)
plot(hd.xaxis,d12.pval.runday.hd,type="s",ylim=c(100,range12),log="y",las=1,ylab="Number of mass features",xlab="P-value",main="Day12 Runday")
lines(hd.xaxis,rep((0.05*10687),20),lty=2)
dev.off()

#pval threshold
#day4
features_d4_01<-calculatePOverlaps(p.val_runday_d4,dataset_sig_features_day4$adjp,0.01)
features_d4_05<-calculatePOverlaps(p.val_runday_d4,dataset_sig_features_day4$adjp,0.05)
features_d4_1<-calculatePOverlaps(p.val_runday_d4,dataset_sig_features_day4$adjp,0.1)

features_d4_pvalThreshold<-rbind(features_d4_01,features_d4_05,features_d4_1)
colnames(features_d4_pvalThreshold)<-c("Common","Runday only","Strain only","threshold")
rownames(features_d4_pvalThreshold)<-c("pval< .01","pval< .05","pval< .1")


#ronly.mz_rt <- strsplit(rownames(ms_data_total), "\\@")
# mz<-sapply(mz_rt , function (x) if(length(x) == 2) x[1] else as.character(NA))
# rt<-sapply(mz_rt , function (x) if(length(x) == 2) x[2] else as.character(NA))

metadata_melt_d4<-melt(features_d4_pvalThreshold,id="threshold")
metadata_melt_d4<-metadata_melt_d4[metadata_melt_d4[,2]!="threshold",]
features_d4_pvalThreshold_plot<-ggplot(metadata_melt_d4,aes(x=Var2,value,fill=Var1))+ geom_point (aes(color=Var1,shape = Var2),size=4) + coord_cartesian(ylim = c(0, 10000)) +
  theme_bw() + theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
                     strip.text.x = element_text(size=12, face="bold"),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black')) + xlab("Factor") + ylab("Number of mass features") + ggtitle("Day4- Linear model")

#day12
features_d12_01<-calculatePOverlaps(p.val_runday_d12,dataset_sig_features_day12$adjp,0.01)
features_d12_05<-calculatePOverlaps(p.val_runday_d12,dataset_sig_features_day12$adjp,0.05)
features_d12_1<-calculatePOverlaps(p.val_runday_d12,dataset_sig_features_day12$adjp,0.1)

features_d12_pvalThreshold<-rbind(features_d12_01,features_d12_05,features_d12_1)
colnames(features_d12_pvalThreshold)<-c("Common","Runday only","Strain only","threshold")
rownames(features_d12_pvalThreshold)<-c("pval< .01","pval< .05","pval< .1")


#ronly.mz_rt <- strsplit(rownames(ms_data_total), "\\@")
# mz<-sapply(mz_rt , function (x) if(length(x) == 2) x[1] else as.character(NA))
# rt<-sapply(mz_rt , function (x) if(length(x) == 2) x[2] else as.character(NA))

metadata_melt_d12<-melt(features_d12_pvalThreshold,id="threshold")
metadata_melt_d12<-metadata_melt_d12[metadata_melt_d12[,2]!="threshold",]
features_d12_pvalThreshold_plot<-ggplot(metadata_melt_d12,aes(x=Var2,value,fill=Var1))+ geom_point (aes(color=Var1,shape = Var2),size=4) + coord_cartesian(ylim = c(0, 10000)) +
  theme_bw() + theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
                     strip.text.x = element_text(size=12, face="bold"),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black')) + xlab("Factor") + ylab("Number of mass features") + ggtitle("Day12- Linear model")

pdf("pvalueDistributionThreshold_linearmodel.pdf",height=8,width=6)
multiplot(features_d4_pvalThreshold_plot, features_d12_pvalThreshold_plot) #)#,   pc1_corrected, pc2_corrected, pc3_corrected cols=2)
dev.off()

######## day12

pdf("linearModel_permutation.pdf",width=8,height=6)
par(mfrow=c(1,2))
hist(dataset_sig_features_day4$adjp,main=paste0("Exponential phase","\n","After correction Strain"),xlab="adjusted p values",ylim=c(0,8000))
hist(dataset_sig_features_day12$adjp,main=paste0("Stationary phase","\n","After correction Strain"),xlab="adjusted p values",ylim=c(0,8000))
#plot(lm_strain_d4_nonzero_loadings_pval_corrected,lm_runday_d4_nonzero_loadings_pval_corrected)
dev.off()

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

plot(dataset_fvalue_d12,dataset_r2_d12)

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

pdf("linearModel_testStatisticV1.pdf",height=12,width=12)

#par(mfrow=c(4,2))
par(mfrow=c(2,2))

# plot(rt_day4,mz_day4,pch=1,main="day 4 Strain",ylab="m/z",xlab=paste0("rt(s)","\n","no of features=",tot_mz_s_day4),col="#00000033")
# points(sigfeat_day4_bc_rt_s,sigfeat_day4_bc_mz_s,pch=1,col="red")
# 
# plot(rt_day12,mz_day12,pch=1,main="day 12 Strain",ylab="m/z",xlab=paste0("rt(s)","\n","no of features=",tot_mz_s_day12),col="#00000033")
# points(sigfeat_day12_bc_rt_s,sigfeat_day12_bc_mz_s,pch=1,col="red")
# 
# plot(rt_day4,mz_day4,pch=1,main="day 4 Strain-Test statistic",ylab="m/z",xlab=paste0("rt(s)","\n","no of features=",length(tmp1)),col="#00000033")
# points(sigfeat_day4_bc_rt_s,sigfeat_day4_bc_mz_s,pch=1,col=ifelse(day4_numdenum_s[day4_numdenum_s_ind]<0.1,"red","#00000033"))
# 
# plot(rt_day12,mz_day12,pch=1,main="day 12 Strain-Test statistic",ylab="m/z",xlab=paste0("rt(s)","\n","no of features=",length(tmp2)),col="#00000033")
# points(sigfeat_day12_bc_rt_s,sigfeat_day12_bc_mz_s,pch=1,col=ifelse(day12_numdenum_s[day12_numdenum_s_ind]<0.1,"red","#00000033"))

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