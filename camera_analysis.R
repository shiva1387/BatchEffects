isotopeData<-read.table('algae_4setblanks_021213_mzrt_camera.txt',header=TRUE,sep='\t')
#isotopeData$mzrt<-paste0(round(as.numeric(isotopeData$mz),5),'@',isotopeData$rt)
isotopeData$mz<-round(as.numeric(isotopeData$mz),5)
isotopeData$rt<-round(as.numeric(isotopeData$mz),0)
isotopeData$mzrt<-paste0(isotopeData$mz,'@',isotopeData$rt)
nonSingleton<-duplicated(isotopeData$pcgroup) | duplicated(isotopeData$pcgroup, fromLast = TRUE)
mzrt_nonSingleton<-isotopeData[which(nonSingleton==TRUE),]
#23775

################

sigfeat_day12_bc_s<-strsplit(rownames(day12_nonzero_sigfeat_s_matrix[[4]]), "\\@")
#9044
sigfeat_day12_bc_mz_s<-sapply(sigfeat_day12_bc_s , function (x) if(length(x) == 2) x[1] else as.character(NA))
sigfeat_day12_bc_rt_s<-sapply(sigfeat_day12_bc_s , function (x) if(length(x) == 2) x[2] else as.character(NA))
sigfeat_day12_bc_mz_s<-round(as.numeric(sigfeat_day12_bc_mz_s),5)
sigfeat_day12_bc_rt_s<-round(as.numeric(sigfeat_day12_bc_rt_s),0)
sigfeat_day12_bc_mzrt_s<-paste0(sigfeat_day12_bc_mz_s,'@',sigfeat_day12_bc_rt_s)
sigfeat_day12_bc<-data.frame(sigfeat_day12_bc_mzrt_s,sigfeat_day12_bc_mz_s,sigfeat_day12_bc_rt_s)
#mzratio<-round(as.numeric(sigfeat_day12_bc_mz_s),5)

significant_nonSingleton_day12<-match(sigfeat_day12_bc_mzrt_s,mzrt_nonSingleton$mzrt); 
significant_nonSingleton_day12<-significant_nonSingleton_day12[!is.na(significant_nonSingleton_day12)]
### 59

sigfeat_day4_bc_s<-strsplit(rownames(day4_nonzero_sigfeat_s_matrix[[7]]), "\\@")
#3979
sigfeat_day4_bc_mz_s<-sapply(sigfeat_day4_bc_s , function (x) if(length(x) == 2) x[1] else as.character(NA))
sigfeat_day4_bc_rt_s<-sapply(sigfeat_day4_bc_s , function (x) if(length(x) == 2) x[2] else as.character(NA))
sigfeat_day4_bc_mz_s<-round(as.numeric(sigfeat_day4_bc_mz_s),5)
sigfeat_day4_bc_rt_s<-round(as.numeric(sigfeat_day4_bc_rt_s),0)
sigfeat_day4_bc_mzrt_s<-paste0(sigfeat_day4_bc_mz_s,'@',sigfeat_day4_bc_rt_s)
sigfeat_day4_bc<-data.frame(sigfeat_day4_bc_mzrt_s,sigfeat_day4_bc_mz_s,sigfeat_day4_bc_rt_s)

# significant_nonSingleton_day4<-match(sigfeat_day4_bc_mzrt_s,mzrt_nonSingleton$mzrt); 
# significant_nonSingleton_day4<-significant_nonSingleton_day4[!is.na(significant_nonSingleton_day4)]
# #20

################ MZ bins

#functions to get matching mz and rt
cmpMZ <- function(mz, mzlist, mz_cutoff){abs(mz-mzlist) <= mz_cutoff}
#if there is more than one match, change function to pick the closest!
cmpRT <- function(ret, retlist,ret_cutoff){abs(ret-retlist) <= ret_cutoff}


#setting tolerance limits--To write as a function!

##day4
mz_tolerance=10 #ppm
rt_tolerance=5 #secs
for(i in 1:length(sigfeat_day4_bc_s))
{
  
  ppm_tolerance<-(trunc(sigfeat_day4_bc_mz_s[i])*mz_tolerance)/10^6 
  match_mz<-cmpMZ(sigfeat_day4_bc_mz_s[i],mzrt_nonSingleton$mz,ppm_tolerance)
  match_mz_ind<-which(match_mz) #obtaining indices only for the significant features
  match_ret<-cmpRT(sigfeat_day4_bc_rt_s[i],mzrt_nonSingleton$rt[match_mz_ind],rt_tolerance)
  mz_rt_match<-match_mz_ind[which(match_ret)]
  matched_value<-mzrt_nonSingleton[mz_rt_match,"mzrt"];
  matched_value=ifelse(exists("matched_value"),matched_value)
  #print(c(i,sigfeat_day4_bc_mz_s[i],sigfeat_day4_bc_rt_s[i],mz_rt_match,matched_value))
  sigfeat_day4_bc[i,4]<-matched_value
}
sigfeat_day4_bc_nonsingletons<-sigfeat_day4_bc[!is.na(sigfeat_day4_bc$V4),]
sigfeat_day4_nonsingletonList_mz<-strsplit(sigfeat_day4_bc_nonsingletons$V4, "\\@")
#writing mz rt to a file
sigfeat_day4_nonsingletonList_mzrt<-do.call(rbind.data.frame, sigfeat_day4_nonsingletonList_mz)
colnames(sigfeat_day4_nonsingletonList_mzrt)<-c("mz","rt(secs)")
sigfeat_day4_nonsingletonList_mzrt$mz<-round(as.numeric(as.character(sigfeat_day4_nonsingletonList_mzrt$mz)),5)
write.table(sigfeat_day4_nonsingletonList_mzrt,"sigfeat_day4_nonsingletonList_mzrt.txt",quote=FALSE,row.names=FALSE)
#writing mz alone to a file
sigfeat_day4_nonsingletonList_mz<-sapply( sigfeat_day4_nonsingletonList_mz, function (x) if(length(x) == 2) x[1] else as.character(NA))
sigfeat_day4_nonsingletonList_mz<-round(as.numeric(sigfeat_day4_nonsingletonList_mz),5)
write.table(sigfeat_day4_nonsingletonList_mz,"sigfeat_day4_nonsingletonList_mz.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
#2263 features found after removing singletons


###day12
mz_tolerance=10 #ppm
rt_tolerance=5 #secs
for(i in 1:length(sigfeat_day12_bc_s))
{
  
  ppm_tolerance<-(trunc(sigfeat_day12_bc_mz_s[i])*mz_tolerance)/10^6 
  match_mz<-cmpMZ(sigfeat_day12_bc_mz_s[i],mzrt_nonSingleton$mz,ppm_tolerance)
  match_mz_ind<-which(match_mz) #obtaining indices only for the significant features
  match_ret<-cmpRT(sigfeat_day12_bc_rt_s[i],mzrt_nonSingleton$rt[match_mz_ind],rt_tolerance)
  mz_rt_match<-match_mz_ind[which(match_ret)]
  matched_value<-mzrt_nonSingleton[mz_rt_match,"mzrt"];
  matched_value=ifelse(exists("matched_value"),matched_value)
  #print(c(i,sigfeat_day12_bc_mz_s[i],sigfeat_day12_bc_rt_s[i],mz_rt_match,matched_value))
  sigfeat_day12_bc[i,4]<-matched_value
}
sigfeat_day12_bc_nonsingletons<-sigfeat_day12_bc[!is.na(sigfeat_day12_bc$V4),]
sigfeat_day12_nonsingletonList_mz<-strsplit(sigfeat_day12_bc_nonsingletons$V4, "\\@")
#writing mz rt to a file
sigfeat_day12_nonsingletonList_mzrt<-do.call(rbind.data.frame, sigfeat_day12_nonsingletonList_mz)
colnames(sigfeat_day12_nonsingletonList_mzrt)<-c("mz","rt(secs)")
sigfeat_day12_nonsingletonList_mzrt$mz<-round(as.numeric(as.character(sigfeat_day12_nonsingletonList_mzrt$mz)),5)
write.table(sigfeat_day12_nonsingletonList_mzrt,"sigfeat_day12_nonsingletonList_mzrt.txt",quote=FALSE,row.names=FALSE)
#writing mz alone to a file
sigfeat_day12_nonsingletonList_mz<-sapply( sigfeat_day12_nonsingletonList_mz, function (x) if(length(x) == 2) x[1] else as.character(NA))
sigfeat_day12_nonsingletonList_mz<-round(as.numeric(sigfeat_day12_nonsingletonList_mz),5)
write.table(sigfeat_day12_nonsingletonList_mz,"sigfeat_day12_nonsingletonList_mz.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
#3396 features found after removing singletons

################ Inflation statistic singletons Vs non-singletons
# here i=4 is the best case for day12, aft batch effect removal
i=4
a<-as.numeric(as.factor(SampleGroup_day12))-1
#b<-mt.teststat.num.denum(as.matrix(svd_day12_nonzero[[4]]),a,test="f",nonpara="n")
sigfeat_day12_bc_s<-strsplit(rownames(day12_nonzero_sigfeat_s_matrix[[i]]), "\\@")

day12_numdenum_s<-day12_nonzero_denum_values_s[,i]

#plotting test statistic vs no-of-sig-features strain
sigfeat_day12_s<-rep('nosig',length(rt_day12))
sigfeat_day12_s[significant_nonSingleton_day12]<-'nonsingle'

df1_s<-data.frame(as.numeric(day12_nonzero_denum_values_s[,i]),sigfeat_day12_s)
colnames(df1_s)<-c('denom','sigfeat')
df1_s$bins <- cut2(df1_s$denom, c(0.025,0.05,0.1,0.5,1,1.5,2,3))
df2_s<-as.data.frame.matrix(table(df1_s$bins,df1_s$sigfeat))
#colnames(df2_s)<-c('bin','nosig','sig')
df2_s$bin<-gsub('\\[|\\)|\\]', '', rownames(df2_s))
splits<-strsplit(df2_s$bin, ",")
df2_s$newBin<-as.numeric(sapply(splits , function (x) if(length(x) == 2) x[2] else as.character(NA)))

plot(log(df2_s$newBin),df2_s$nosig,type='l', main="Strain-test statistic variance Vs Sig features", ylim=c(0,max(df2_s$nosig,df2_s$sig)))
lines(log(df2_s$newBin),df2_s$sig,col="red")
legend(x="topright",inset=0, c("sig","non-sig"),col=c("red","black"),lty=1)    
