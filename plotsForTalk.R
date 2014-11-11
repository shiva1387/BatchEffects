##Linear model and Analysis of distance, PC axis Versus Time Stamp
##Author: Shiv
##Version :04-08-14

### Library

library('ggplot2')
library('RColorBrewer')

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

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#### Reading data

RunDay_TimeStamp_d12<-read.table("RunDay_TimeStamp_day12.txt",header=TRUE,sep='\t')
RunDay_TimeStamp_d12<-as.data.frame(RunDay_TimeStamp_d12)
rownames(RunDay_TimeStamp_d12)<-RunDay_TimeStamp_d12[,]
RunDay_TimeStamp_d12<-RunDay_TimeStamp_d12[colnames(ms_data_day12_nonzero),]

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

### P adjusted 

lm_strain_d4_pval_raw_adj<-round(p.adjust(lm_strain_d4_pval_raw, "BH"),3)
lm_runday_d4_pval_raw_adj<-round(p.adjust(lm_runday_d4_pval_raw, "BH"),3)
lm_strain_d4_pval_corrected_adj<-round(p.adjust(lm_strain_d4_pval_corrected, "BH"),3)
lm_runday_d4_pval_corrected_adj<-round(p.adjust(lm_runday_d4_pval_corrected, "BH"),3)

pdf("lm_day4_padj.pdf",height=8,width=8)
par(mfrow=c(2,2))
hist(lm_strain_d4_pval_raw_adj,main="Raw data Strain",ylim=c(0,14000))
hist(lm_runday_d4_pval_raw_adj,main="Raw data Runday",ylim=c(0,14000))
hist(lm_strain_d4_pval_corrected_adj,main="After correction Strain",ylim=c(0,14000))
hist(lm_runday_d4_pval_corrected_adj,main="After correction Runday",ylim=c(0,14000))
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

data<-ms_data_day12_nonzero
processed_data<-normalize.quantiles(as.matrix(data),copy=TRUE)
processed_data<-scale(data,center=T,scale=T)
processed_data<-processed_data-min(processed_data)

####
### The following section presents analysis of distance for the full dataset(including zeroes)
### The dataset was used after scaling using the above method

ms_data_day12_s<-ScaleData(ms_data_day12)
ms_data_day4_s<-ScaleData(ms_data_day4)

###DAY12
# > a<-adonis(t(ms_data_day12_s)~ RunDay_day12, method = "bray", perm=999)
# > a
# 
# Call:
#   adonis(formula = t(ms_data_day12_s) ~ RunDay_day12, permutations = 999,      method = "bray") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_day12   2    1.6373 0.81867  14.859 0.20535  0.001 ***
#   Residuals    115    6.3362 0.05510         0.79465           
# Total        117    7.9735                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# > a<-adonis(t(ms_data_day12_s)~ RunDay_day12, method = "euclidean", perm=999)
# > a
# Call:
#   adonis(formula = t(ms_data_day12_s) ~ RunDay_day12, permutations = 999,      method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_day12   2    324496  162248    9.76 0.14511  0.001 ***
#   Residuals    115   1911727   16624         0.85489           
# Total        117   2236223                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###DAY4
# > a<-adonis(t(ms_data_day4_s)~ RunDay_day4, method = "bray", perm=999)
# > a
# 
# Call:
#   adonis(formula = t(ms_data_day4_s) ~ RunDay_day4, permutations = 999,      method = "bray") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_day4   3    1.3791 0.45971  9.8792 0.19675  0.001 ***
#   Residuals   121    5.6304 0.04653         0.80325           
# Total       124    7.0096                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# > a<-adonis(t(ms_data_day4_s)~ RunDay_day4, method = "euclidean", perm=999)
# > a
# Call:
#   adonis(formula = t(ms_data_day4_s) ~ RunDay_day4, permutations = 999,      method = "euclidean") 
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# RunDay_day4   3    203530   67843  5.3115 0.11637  0.001 ***
#   Residuals   121   1545519   12773         0.88363           
# Total       124   1749049                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

####
### The following section presents analysis of distance for the complete.dataset(after removing features which had zeroes)
###


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

############################################################################################################
## Estimate correlation with PCOA and within runday-time (Modified by Shiv on 29-07-14)

compute_pca<-function(dataset,preprocess_method) {
  #dataset<-as.matrix(dataset) #For batch corrected matrix, log transformation is not possible
  dataset<-as.matrix(log(dataset)) #log transform the data using natural log
  if(preprocess_method=="norm") {
    processed_data<-normalize.quantiles(as.matrix(dataset),copy=TRUE)
  }   else  {
    processed_data<-scale(dataset,center=T,scale=T)
    processed_data<-processed_data-min(processed_data)
  }
  pca_results <- princomp(processed_data,cor=F,scores=T) ### IMP: choose quantile normalized or scaled data
  return(pca_results)
}

compute_linearModel<-function(results.from.pca,dependent.factor) { #dependent.factor is either RunDay_day4 or 12 (or) SampleGroup_day4 or 12 
  lm_pca_scores<-apply(results.from.pca$loadings,2, function(x) {
    lm_val<-lm(x~ as.factor(dependent.factor)) 
    lm_cor<-summary(lm_val)
    p.val<-anova(lm_val)$'Pr(>F)'[1]
    return(list(lm_cor$r.squared,p.val))
  })
}

#function to extract r2 value from list containing r2 and p.val returned from linear model
compute.r2.pval<-function(linearmodel_list,r2.pval) {
  if(r2.pval=="r2") { #WARNING:code implictly assumes r2 is in the first column and p.val in the second
    return (sapply(linearmodel_list, function(x){as.numeric(x[1])}))
  } else{
    return (sapply(linearmodel_list, function(x){as.numeric(x[2])}))
  }
}

#For batch corrected matrix, log transformation is not possible
#### IMPORTANT --COMPUTE FUNTION AGAIN IF CHANGING BETWEEN BATCH CORRECTED MATRIX AND RAW DATA
fit_day4_RunDay_TimeStamp_mzbysam_princomp<-compute_pca(ms_data_day4_nonzero,"scale") #Check compute_pca function  to verify that log transformation is active before running

# batch_corrected_mat_d12
# batch_corrected_mat_d4
# ms_data_day12_nonzero
# ms_data_day4_nonzero

# lm_pca_scores_RunDay_TimeStamp_day12_nonzero_loadings<-compute_linearModel(fit_day12_RunDay_TimeStamp_mzbysam_princomp,as.factor(RunDay_TimeStamp_d12$Time))
# lm_pca_scores_RunDay_TimeStamp_nonzero_loadings_r.sq<-compute.r2.pval(lm_pca_scores_RunDay_TimeStamp_day12_nonzero_loadings,"r2")
# lm_pca_scores_RunDay_TimeStamp_nonzero_loadings_pval<-compute.r2.pval(lm_pca_scores_RunDay_TimeStamp_day12_nonzero_loadings,"pval")

### plotting PC1 Vs Time Stamp

##reading in date time information

day4_dateTime<-read.table("day4_dateTime.txt",row.names=1,header=TRUE)
day4_dateTime<-day4_dateTime[colnames(ms_data_day4_nonzero),]
day4_dateTime$dateTime <- strptime(paste(day4_dateTime$date, day4_dateTime$time), "%m/%d/%Y %H:%M:%S")
day4_dateTime$date <- as.Date(as.character(day4_dateTime$date),format="%m/%d/%Y")
day4_dateTime$time <- as.POSIXct(as.character(day4_dateTime$time),format="%H:%M:%S")

fit_day4_RunDay_TimeStamp_mzbysam_princomp_loadings<-fit_day4_RunDay_TimeStamp_mzbysam_princomp$loadings

plot(day4_dateTime$dateTime,fit_day4_RunDay_TimeStamp_mzbysam_princomp_loadings[,1],xaxt="n")
#axis.POSIXct(1,at=day4_dateTime$dateTime,labels=format(day4_dateTime$dateTime,"%b %d"),las=2)
axis.POSIXct(1, at=day4_dateTime$dateTime, labels=format(day4_dateTime$dateTime, "%H:%M %b %d"))

# ggplot 
test<-data.frame(Date = day4_dateTime$dateTime,value=fit_day4_RunDay_TimeStamp_mzbysam_princomp_loadings[,3], SampleGroup=SampleGroup_day4)
#Choosing color
colourCount = length(unique(test$SampleGroup))
getPalette = colorRampPalette(brewer.pal(10,"Paired"))
set.seed(1) #important to set seed so that we obtain the same shapes for strains all the time
pch_types<-c(15, 16, 17, 18, 25, 8)
pch_values<-sample(pch_types, 22, replace = TRUE)

##plotting
plot1<- ggplot(test, aes(Date, value,colour= factor(SampleGroup), shape = factor(SampleGroup))) + geom_point(size=2) 
plot2 <- plot1 + scale_colour_manual('Strains', values=getPalette(colourCount)) + scale_shape_manual('Strains',values=pch_values)
plot3<- plot2+ theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1,size=12),axis.text.y=element_text(size=12),axis.title.y=element_blank(),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black'))+ xlab("Run TimeStamp") + ylab("PC 3 loadings") +  ggtitle("Raw PC3 Vs TimeStamp- day4") ####RAW  Batch effect removed


pdf("ms_data_d12_raw_PC3VsTimeStamp.pdf",height=8,width=16, paper="a4")
plot3
dev.off()

# pc1_raw<-plot3
# pc2_raw<-plot3
# pc3_raw<-plot3

# pc1_corrected<-plot3
# pc2_corrected<-plot3
# pc3_corrected<-plot3

pdf("ms_data_d4_multiplotVsTimeStamp_raw.pdf",height=8,width=12)
multiplot(pc1_raw, pc2_raw, pc3_raw) #)#,   pc1_corrected, pc2_corrected, pc3_corrected cols=2)
dev.off()


##############
## To extract the hours and minutes and catergorize them into 4 batches
timeBins<-read.table('timeBins.txt',header=FALSE)
timeBins<-as.data.frame(timeBins)
timeBins_split <- strsplit(as.character(timeBins$V1), ":")
timeBins$time<-sapply(timeBins_split , function (x) if(length(x) == 3) paste0(x[1],'.',x[2]) else as.character(NA)) #Extracting only the HRS and MINS
timeBins$time<-as.numeric(timeBins$time)
timeBins$bin<-ifelse(timeBins$time >= 0 & timeBins$time <= 6,'A', #first 6 hours
                     ifelse(timeBins$time > 6 & timeBins$time <= 12,'B', #second 6 hours
                            ifelse(timeBins$time > 12 & timeBins$time <= 18,'C','D'))) #third and fourth 6 hours
write.table(timeBins,"timeBins_mod.txt",quote=FALSE,sep="\t")


