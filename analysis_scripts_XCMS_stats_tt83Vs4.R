#--analysis log for preprocessing of Shiv TT8 experiment--
#--started on 16 May 2013--
#-- Modified by Shiv for TT8
#--Author: Rohan
#--Comparing tt8 4 bioreps with 3 bioreps
#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 

library(ggplot2)
library(reshape)
library(xcms)
library(stringr)
library(vegan)
library(ape)
library(gplots)
library(AStream)
library(marray)
library(qvalue)


#############
# User      #
# Specific  #
# Variables #
#############

directory<-"D:/Shiv_data/Shiv_NUS/PhD_Project/Integrating_Omics/Analysis_R/MetabolomicsSoftwareComparison/Data/TT8_RawData_Metabolomics/Agilent/MZXML_PEAK/"; 
#Contains path for .tsv file

#directory<-"D:/Shiv_data/Shiv_NUS/PhD_Project/Integrating_Omics/Analysis_R/MetabolomicsSoftwareComparison/Data/TT8_RawData_Metabolomics/Agilent/MZXML_MS1/RemovedNoisyBiorep/WS/"
#Contains path for raw mzxml files

# The "PATH", Remember to use "/" and not "/" in the path
num_cores<-8; # The number of CPUs available, be aware of memory issue

setwd(directory)

####### reading in the mass spec data 

mzfilename<-"XCMS_TT8_21052013_PEAK.tsv"
filename<-gsub(".tsv","",mzfilename)

ms_data_total<-read.table(mzfilename,sep="\t",header=T,check.names=FALSE)
xsannotated<-ms_data_total
str(ms_data_total)

xsannotated<-ms_data_total

no_of_info_cols<-9
data_start_index<-no_of_info_cols+1;
data_end_index<-dim(ms_data_total)[2];

ms_data<-ms_data_total[,data_start_index:data_end_index]
names(ms_data)
dim(ms_data)
xsannotated1<-as.matrix(ms_data)


ms_data_avg_biorep<-byapply(ms_data,3,rowMeans)
colnames(ms_data_avg_biorep)<-c("TT8-B1","TT8-B2","TT8-B3","WS-B1","WS-B2","WS-B3")
ms_data_avg_biorep<-as.data.frame(ms_data_avg_biorep)

ms_data_tt8_ws<-byapply(ms_data_avg_biorep,3,rowMeans)
colnames(ms_data_tt8_ws)<-c("TT8","WS")
ms_data_tt8_ws<-as.data.frame(ms_data_tt8_ws)

#--load in Rasmus' annotated version of the data...
#xsannotated<-read.table(file="xsannotated.csv",as.is=T,header=T,blank.lines.skip=T,sep=",")
#> dim(xsannotated)
#[1] 2260   72
#--check against file..there are 2261 including the header. Fine.

#--these are the samples, they are grouped in 3's by technical replicates...

#> colnames(xsannotated1)

#--most important thing to plot up is the m/z-by-r/t map, not coding size of peak in this case...
#plot(xsannotated[,"mz"],xsannotated[,"rt"],xlab="m/z",ylab="rt",las=1,main=filename)
#dev.print(pdf,file="mass-features-plot.pdf")

#--and plot each "marginal"...
#par(mfrow=c(2,2))
#plot(xsannotated[,"rt"],apply(xsannotated1,1,sum),xlab="rt",ylab="Summed counts",las=1,type="p",pch=3)
#plot(xsannotated[,"mz"],apply(xsannotated1,1,sum),xlab="m/z",ylab="Summed counts",las=1,type="p",pch=3)
#plot(xsannotated[,"rt"],apply(xsannotated1,1,sum),xlab="rt",ylab="Summed counts",las=1,type="p",pch=3,log="y")
#plot(xsannotated[,"mz"],apply(xsannotated1,1,sum),xlab="m/z",ylab="Summed counts",las=1,type="p",pch=3,log="y")
#par(mfrow=c(1,1))

#--make a log10 version...
#xsannotated1.log10<-log10(xsannotated1+1)

#--boxplots...
# boxplot(xsannotated1.log10,xlab="Timepoints",ylab="Log2 feature counts",las=1,xaxt="n")
# abline(v=seq(1,18,3)-0.5)
# axis(side=1,at=seq(2,18,3),as.character(1:6))
# dev.print(pdf,file="boxplots.pdf")

#--measures of variability...
#--overall standard deviation across time for each metabolite...
tt8_4br<-log10(ms_data+1) # XCMS_TT8_29042013.tsv
tt8_3br<-log10(ms_data1+1) # XCMS_TT8_14052013.tsv
par(mfrow=c(2,1))
tt8_12<-log10(tt8_4br[,1:12]+1)
tt8_9<-log10(tt8_3br[,1:9]+1)
plot(apply(tt8_12,1,mean),apply(tt8_12,1,sd),pch=1,col=2,xlab="mean within bioreps",ylab="sd within bioreps")
points(apply(tt8_9,1,mean),apply(tt8_9,1,sd),pch=1,col=3)
legend("topleft", c("tt8_4brs","tt8_3brs"), cex=0.8, pch=1 ,col=c(2,3))

ws_12<-log10(tt8_4br[,13:24]+1)
ws_9<-log10(tt8_3br[,10:18]+1)
plot(apply(ws_12,1,mean),apply(ws_12,1,sd),pch=1,col=2,xlab="mean within bioreps",ylab="sd within bioreps")
points(apply(ws_9,1,mean),apply(ws_9,1,sd),pch=1,col=3)
legend("topleft", c("ws_4brs","ws_3brs"), cex=0.8, pch=1 ,col=c(2,3))


#[1] 6660
#3659 (<0.005)

##################

#xsannotated1.log10.row.sd<-apply(xsannotated1.log10,1,sd)

#--between time-point variability...
#--need a vector of memberships of timepoints...
#labels.for.tr.4br<-sort(rep(1:2,12))
#labels.for.tr.3br<-sort(rep(1:2,9))
#here 1 is tt8 and 2 is ws

#--write a function to compute these...
#mean.by.tr<-function(x,cl){tapply(x,INDEX=cl,FUN=mean)}
#sd.by.tr<-function(x,cl){tapply(x,INDEX=cl,FUN=sd)}
#xsannotated1.log10.row.sd.of.tp.means<-apply(t(apply(xsannotated1.log10,1,mean.by.tr,cl=labels.for.tr)),1,sd)
#xsannotated1.log10.row.mean.of.tp.sd<-apply(t(apply(xsannotated1.log10,1,sd.by.tr,cl=labels.for.tr)),1,mean)


tt8_4br.log10.row.sd.of.tp.means<-apply(t(apply(tt8_4br,1,mean.by.tr,cl=labels.for.tr.4br)),1,sd)
tt8_4br.log10.row.mean.of.tp.sd<-apply(t(apply(tt8_4br,1,sd.by.tr,cl=labels.for.tr.4br)),1,mean)


tt8_3br.log10.row.sd.of.tp.means<-apply(t(apply(tt8_3br,1,mean.by.tr,cl=labels.for.tr.3br)),1,sd)
tt8_3br.log10.row.mean.of.tp.sd<-apply(t(apply(tt8_3br,1,sd.by.tr,cl=labels.for.tr.3br)),1,mean)

plot(apply(tt8_4br,1,mean),tt8_4br.log10.row.sd.of.tp.means,pch=1,col=2,,xlab="Mean read count (log10)",ylab="Between treatment variability")
points(apply(tt8_3br,1,mean),tt8_3br.log10.row.sd.of.tp.means,pch=1,col=3)
legend("topleft", c("tt8_4br","tt8_3br"), cex=0.8, pch=1 ,col=c(2,3))

t1<-tt8_3br[1,]
s1<-sum(t1[1,1:9])/9
s2<-sum(t1[1,10:18])/9
s3<-c(s1,s2)
s4<-sd(s3)
apply(t(apply(t1,1,mean.by.tr,cl=labels.for.tr.3br)),1,sd)

### mean vairance relationships

tt8_4br.log10.row.inter.to.intra.tp.var<-tt8_4br.log10.row.sd.of.tp.means/tt8_4br.log10.row.mean.of.tp.sd
tt8_3br.log10.row.inter.to.intra.tp.var<-tt8_3br.log10.row.sd.of.tp.means/tt8_3br.log10.row.mean.of.tp.sd

par(mfrow=c(2,2))
plot(apply(tt8_4br,1,mean),apply(tt8_4br,1,sd),pch=1,col=2,ylim=c(0,2.6),xlab="Mean read count (log10)",ylab="SD",las=1)
points(apply(tt8_3br,1,mean),apply(tt8_3br,1,sd),pch=1,col=3)
legend("topright", c("tt8_4br","tt8_3br"), cex=0.8, pch=1 ,col=c(2,3))

plot(apply(tt8_4br,1,mean),tt8_4br.log10.row.sd.of.tp.means,ylim=c(0,2.5),pch=1,col=2,xlab="Mean read count (log10)",ylab="Between treatment variability",las=1)
points(apply(tt8_3br,1,mean),tt8_3br.log10.row.sd.of.tp.means,pch=1,col=3)

plot(apply(tt8_4br,1,mean),tt8_4br.log10.row.mean.of.tp.sd,ylim=c(0,2.5),pch=1,col=2,xlab="Mean read count (log10)",ylab="Within treatment variability",las=1)
points(apply(tt8_3br,1,mean),tt8_3br.log10.row.mean.of.tp.sd,pch=1,col=3)

plot(apply(tt8_4br,1,mean),tt8_4br.log10.row.inter.to.intra.tp.var,log="y",las=1,pch=1,col=2,xlab="Mean read count (log10)",ylab="Between:within treatment variability",las=1)
points(apply(tt8_3br,1,mean),tt8_3br.log10.row.inter.to.intra.tp.var,pch=1,col=3)
#abline(h=1,lty=2)
par(mfrow=c(1,1))

########### ttest and fdr

tt8_4br.p.value = apply(tt8_4br, 1, function(x) { t.test(x[1:12], x[13:24]) $p.value } )
tt8_3br.p.value = apply(tt8_3br, 1, function(x) { t.test(x[1:9], x[10:18]) $p.value } )
# Check the first few brain ones to make sure the first one agrees with our single-gene command

tt8_4br.p.value[1:5]
tt8_4br.fdr.pvals = p.adjust(tt8_4br.p.value, method="fdr")
tt8_3br.fdr.pvals = p.adjust(tt8_3br.p.value, method="fdr")

# features pval less than 0.05

length(tt8_3br.fdr.pvals[tt8_3br.fdr.pvals<0.05])
#[1] 5784
#3433(<0.005)

# str(intersect(which(tt8_4br.fdr.pvals<0.05),which(tt8_4br.log10.row.inter.to.intra.tp.var>5)))
# int [1:87] 162 173 182 250 655 689 715 889 1015 1082 ...
# 
# head(intersect(which(tt8_3br.fdr.pvals<0.05),which(tt8_3br.log10.row.inter.to.intra.tp.var>5)))
# int [1:112] 69 171 181 190 261 391 432 707 741 768 ...

tt8_diff_ions<-ms_data_total[intersect(which(tt8_3br.fdr.pvals<0.05),which(tt8_3br.log10.row.inter.to.intra.tp.var>5)),]

heatmap.2(as.matrix(tt8_3br[intersect(which(tt8_3br.fdr.pvals<0.05),which(tt8_3br.log10.row.inter.to.intra.tp.var>5)),]),dendrogram="column",trace="none",mar=c(10,10),col=greenred, xlab = NULL)

#--compute ratio of these...
#xsannotated1.log10.row.inter.to.intra.tp.var<-xsannotated1.log10.row.sd.of.tp.means/xsannotated1.log10.row.mean.of.tp.sd

#--plot this against our other data...
# #> par()$mar
# [1] 5.1 4.1 4.1 2.1
# par(mfrow=c(3,2))
# par(mar=c(5.1,6.1,4.1,2.1))
# plot(xsannotated[,"rt"],apply(xsannotated1,1,sum),xlab="rt",ylab="",las=1,type="p",pch=3,cex.axis=1)
# mtext("Counts",side=2,at=1.25e+09,line=4,cex=0.80)
# plot(xsannotated[,"mz"],apply(xsannotated1,1,sum),xlab="m/z",ylab="",las=1,type="p",pch=3)
# mtext("Counts",side=2,at=1.25e+09,line=4,cex=0.80)
# plot(xsannotated[,"rt"],apply(xsannotated1,1,sum),xlab="rt",ylab="",las=1,type="p",pch=3,log="y")
# mtext("Counts (log10)",side=2,at=1.25e+07,line=4,cex=0.80)
# plot(xsannotated[,"mz"],apply(xsannotated1,1,sum),xlab="m/z",ylab="",las=1,type="p",pch=3,log="y")
# mtext("Counts (log10)",side=2,at=1.25e+07,line=4,cex=0.80)
# plot(xsannotated[,"rt"],xsannotated1.log10.row.inter.to.intra.tp.var,xlab="rt",ylab="",las=1,type="p",pch=3,log="y")
# abline(h=1,lty=2)
# mtext("Between:within treatment variability",side=2,at=5,line=4,cex=0.80)
# plot(xsannotated[,"mz"],xsannotated1.log10.row.inter.to.intra.tp.var,xlab="m/z",ylab="",las=1,type="p",pch=3,log="y")
# abline(h=1,lty=2)
# mtext("Between:within treatment variability",side=2,at=5,line=4,cex=0.80)
# par(mfrow=c(1,1))
# dev.print(pdf,file="marginal-counts-plots.pdf")

#--what about mean-variance relationships?
#par(mfrow=c(2,2))
#plot(apply(xsannotated1.log10,1,mean),apply(xsannotated1.log10,1,sd),ylim=c(0,2.6),xlab="Mean read count (log10)",ylab="SD",las=1)
#plot(apply(xsannotated1.log10,1,mean),xsannotated1.log10.row.sd.of.tp.means,ylim=c(0,2.5),xlab="Mean read count (log10)",ylab="Between treatment variability",las=1)
#plot(apply(xsannotated1.log10,1,mean),xsannotated1.log10.row.mean.of.tp.sd,ylim=c(0,2.5),xlab="Mean read count (log10)",ylab="Within treatment variability",las=1)
#plot(apply(xsannotated1.log10,1,mean),xsannotated1.log10.row.inter.to.intra.tp.var,log="y",las=1,xlab="Mean read count (log10)",ylab="Between:within treatment variability",las=1)
#abline(h=1,lty=2)
#par(mfrow=c(1,1))
#dev.print(pdf,file="mean-variance-relationships.pdf")
 
#--try clustering samples...
#plot(hclust(dist(t(xsannotated1.log10)),method="average"),labels=labels.for.tr,las=1,ann=F)
#dev.print(pdf,file="naive-hclust.pdf")

#--this is hard to actually to assess how inter-related the samples are...try using a cophenetic appraoch instead?
#--make a cophenetic matrix from the class label vectors...

#tp.membership.mat<-matrix(NA,nrow=length(labels.for.tr),ncol=length(labels.for.tr))
#for(i in (1:length(labels.for.tr)))
#{
# for(j in (1:length(labels.for.tr)))
# {
#  tp.membership.mat[i,j]<-(labels.for.tr[i]-labels.for.tr[j])
# }
#}

#--write a function to process this entire concept...
compute.distance.analysis<-function(datamat,tpmmat,method="euclidean",rows,...)
{
 d<-datamat[rows,]
 dmat<-as.matrix(dist(t(d),method=method))
 res<-cor(as.vector(dmat),as.vector(abs(tpmmat)),method="spearman")
 return(res)
}
 
#--compute across the range of between:within time point variation...
compute.cor.series<-function(datamat,tpmmat,method="euclidean",testvec,sseq)
{
 #--compute correlation between 2 distance matrices across a range of thresholds
 res1=NULL
 res2=NULL
 
 for(k in (1:length(sseq)))
 {
  tmp1<-compute.distance.analysis(datamat=datamat,tpmmat=tpmmat,method=method,rows=which(testvec>sseq[k]))
  tmp2<-length(which(testvec>sseq[k]))
  res1<-c(res1,tmp1)
  res2<-c(res2,tmp2)
 }
 res<-cbind(threshold=sseq,cor=res1,no.features=res2)
 return(res)
}
 
#matCorExample<-compute.cor.series(xsannotated1.log10,tp.membership.mat,"euclidean",log10(xsannotated1.log10.row.inter.to.intra.tp.var),seq(-0.43,1.57,0.05)) 
#par(mfrow=c(2,1))
#plot(10^matCorExample[,1],matCorExample[,2],type="b",pch=3,log="x",las=1)
#abline(v=1,lty=2)
#plot(10^matCorExample[,1],matCorExample[,3],type="b",pch=3,log="xy",las=1)
#par(mfrow=c(1,1))

#--plot this up as a 2-by-2 panel plot...
# m<-t(matrix(c(1:4),2,2))
# image(as.matrix(dist(t(xsannotated1.log10),method="euclidean")),las=1,xaxt="n",yaxt="n",xlab="Samples (time points)",ylab="Samples (time points)")
# axis(side=1,at=(seq(1,18,3)+0.5)/18,labels=as.character(1:6),cex.axis=0.50)
# axis(side=2,at=(seq(1,18,3)+0.5)/18,labels=as.character(1:6),cex.axis=0.50,las=1)
# image(abs(tp.membership.mat),las=1,xaxt="n",yaxt="n",xlab="Samples (time points)",ylab="Samples (time points)")
# axis(side=1,at=(seq(1,18,3)+0.5)/18,labels=as.character(1:6),cex.axis=0.50)
# axis(side=2,at=(seq(1,18,3)+0.5)/18,labels=as.character(1:6),cex.axis=0.50,las=1)
# plot(abs(tp.membership.mat),as.matrix(dist(t(xsannotated1.log10[which(log10(xsannotated1.log10.row.inter.to.intra.tp.var)>log10(5)),]),method="euclidean")),pch=3,las=1,xlab="Difference in time",ylab="Difference in metabolite profile")
# plot(10^matCorExample[,1],matCorExample[,2],type="b",pch=3,log="x",las=1,xlab="Between:within time point variability",ylab="Mantel correlation statistic")
# dev.print(pdf,file="distance-matrix-correlation-analysis.pdf")

#--this looks much better... 
#plot(hclust(dist(t(xsannotated1.log10[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),])),method="average"),labels=labels.for.tr,las=1,ann=F) 
#dev.print(pdf,file="tp-var-subsetted-hclust.pdf")

#--which features are these?
#plot(xsannotated[,"mz"],xsannotated[,"rt"],xlab="m/z",ylab="rt",las=1,pch=16,col="lightgrey")
#points(xsannotated[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),"mz"],xsannotated[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),"rt"],pch=16,col=1)
#dev.print(pdf,file="mass-features-plot-reproducible.pdf")

#--but are these really biological informative? this suggests they are rather boring?
#matplot(scale(t(xsannotated1.log10[which(xsannotated1.log10.row.inter.to.intra.tp.var>5)[1:10],]),T,F),col=1,type="l",lty=1)
#--plot some examples...
#plot(labels.for.tr,xsannotated1.log10[which(xsannotated1.log10.row.inter.to.intra.tp.var>5)[300],],pch=16,ylim=c(0,8))
#lines(tapply(xsannotated1.log10[which(xsannotated1.log10.row.inter.to.intra.tp.var>5)[300],],INDEX=labels.for.tr,FUN=mean))

#--have a look at the whole data set?
#--average within data points...
xsannotated1.ave.within.tp<-t(apply(xsannotated1,1,mean.by.tr,cl=labels.for.tr))
xsannotated1.log10.ave.within.tp<-t(apply(xsannotated1.log10,1,mean.by.tr,cl=labels.for.tr))
matplot(scale(t(xsannotated1.log10.ave.within.tp),T,F),lty=1,col=1,type="l",xaxt="n",las=1,ylab="Double centered mass-feature log10 counts")
axis(side=1,at=1:2,as.character(1:2))
abline(v=13,lty=2)
dev.print(pdf,file="matplot-dc-profiles-across-treatment.pdf")

#--time for some PCA methinks...
#--to emphasise changes between mass features over time, we should probably use double-centered data.
#--simply visualise changes using icon plots...
double.center.data.matrix<-function(m)
{
 rM<-rowMeans(m)
 cM<-colMeans(m)
 mM<-mean(as.vector(m))

 res<-matrix(NA,nrow(m),ncol(m))
 for(i in (1:nrow(m)))
 {
  for(j in (1:ncol(m)))
  {
   res[i,j]<-m[i,j]-rM[i]-cM[j]+mM
  }
 }
 
 rownames(res)<-rownames(m)
 colnames(res)<-colnames(m)
 
 return(res)
 
}

#--compute double-centered matrix...
#xsannotated1.log10.dc<-double.center.data.matrix(xsannotated1.log10)
#xsannotated1.log10.ave.within.tp.dc<-double.center.data.matrix(xsannotated1.log10.ave.within.tp)

#> range(colMeans(xsannotated1.log10.dc))
#[1] -7.763701e-16  5.639540e-17
#> range(rowMeans(xsannotated1.log10.dc))
#[1] -8.289665e-16  5.921189e-17

#--compute SVD...
#--for 6 sample case...
#xsannotated1.log10.dc.svd<-svd(xsannotated1.log10.dc)
#--get an approximation using pc1 and pc2...
#xsannotated1.log10.dc.svd.pc12<-cbind(xsannotated1.log10.dc.svd$u[,1]*xsannotated1.log10.dc.svd$d[1],xsannotated1.log10.dc.svd$u[,2]*xsannotated1.log10.dc.svd$d[2])
#xsannotated1.log10.dc.svd.pc23<-cbind(xsannotated1.log10.dc.svd$u[,2]*xsannotated1.log10.dc.svd$d[2],xsannotated1.log10.dc.svd$u[,3]*xsannotated1.log10.dc.svd$d[3])
#--and for 6 sample case (averaged within timepoints)
#xsannotated1.log10.ave.within.tp.dc.svd<-svd(xsannotated1.log10.ave.within.tp.dc)
#--get an approximation using pc1 and pc2...
#xsannotated1.log10.ave.within.tp.dc.svd.pc12<-cbind(xsannotated1.log10.ave.within.tp.dc.svd$u[,1]*xsannotated1.log10.ave.within.tp.dc.svd$d[1],xsannotated1.log10.ave.within.tp.dc.svd$u[,2]*xsannotated1.log10.ave.within.tp.dc.svd$d[2])
#xsannotated1.log10.ave.within.tp.dc.svd.pc23<-cbind(xsannotated1.log10.ave.within.tp.dc.svd$u[,2]*xsannotated1.log10.ave.within.tp.dc.svd$d[2],xsannotated1.log10.ave.within.tp.dc.svd$u[,3]*xsannotated1.log10.ave.within.tp.dc.svd$d[3])

#--try this but for rows that meet our variability criteria...
#xsannotated1.log10.ave.within.tp.dc.with.var.criteria.svd<-svd(xsannotated1.log10.ave.within.tp.dc[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),])
#xsannotated1.log10.ave.within.tp.dc.with.var.criteria.svd.pc12<-cbind(xsannotated1.log10.ave.within.tp.dc.with.var.criteria.svd$u[,1]*xsannotated1.log10.ave.within.tp.dc.with.var.criteria.svd$d[1],xsannotated1.log10.ave.within.tp.dc.with.var.criteria.svd$u[,2]*xsannotated1.log10.ave.within.tp.dc.with.var.criteria.svd$d[2])

#--are these consistent with one another?
#--how to test this?!?!?! common principal components?

#--calculate variance and total percent variance per component...
compute.variance.per.pc<-function(svdObj){(svdObj$d^2)/nrow(svdObj$u)}
compute.cum.per.total.var<-function(svdObj){cumsum(svdObj$d^2)/sum(svdObj$d^2)*100}

#xsannotated1.log10.dc.svd.cumpertotalvar<-compute.cum.per.total.var(xsannotated1.log10.dc.svd)
#> xsannotated1.log10.dc.svd.cumpertotalvar[1:2]
#[1] 18.03405 32.16553

#xsannotated1.log10.ave.within.tp.dc.svd.cumpertotalvar<-compute.cum.per.total.var(xsannotated1.log10.ave.within.tp.dc.svd)
#> xsannotated1.log10.ave.within.tp.dc.svd.cumpertotalvar[1:2]
#[1] 29.09748 44.56719

#xsannotated1.log10.ave.within.tp.dc.with.var.criteria.svd.cumpertotalvar<-compute.cum.per.total.var(xsannotated1.log10.ave.within.tp.dc.with.var.criteria.svd)
#> xsannotated1.log10.ave.within.tp.dc.with.var.criteria.svd.cumpertotalvar[1:2]
#[1] 40.65605 61.56712

#--plot up an example PCA...
#plot(xsannotated1.log10.ave.within.tp.dc.svd.pc12,pch=16,main="Double centered, all mass features",xlab="PC1",ylab="PC2",las=1)
#abline(h=0,lty=2)
#abline(v=0,lty=2)
#dev.print(pdf,file="pca1.pdf")

#--do some icon plots..
#code in ge-biplot.R file
#GE.Icon(xsannotated1.log10.ave.within.tp.dc.svd.pc12,"Double centered, all mass features",xsannotated1.log10.ave.within.tp.dc,xlab="PC1",ylab="PC2",las=1)
#abline(h=0,lty=2)
#abline(v=0,lty=2)
#dev.print(pdf,file="icon1.pdf")

#GE.Icon(xsannotated1.log10.ave.within.tp.dc.svd.pc12[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),],"Double centered, variance ratio filtered",xsannotated1.log10.ave.within.tp.dc[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),],xlab="PC1",ylab="PC2",las=1)
#abline(h=0,lty=2)
#abline(v=0,lty=2)
#dev.print(pdf,file="icon2.pdf")

#GE.Icon(xsannotated1.log10.ave.within.tp.dc.with.var.criteria.svd.pc12,"Double centered, variance ratio filtered (only)",xsannotated1.log10.ave.within.tp.dc[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),],xlab="PC1",ylab="PC2",las=1)
#abline(h=0,lty=2)
#abline(v=0,lty=2)
#dev.print(pdf,file="icon3.pdf")

#--seems to be informative, but we still have the problem of how to select samples...
#--try using the polar coordinate concept you developed in late 2011--

pca.convert.to.polar<-function(m)
{
 #--takes a 2-column matrix with projections of data onto 1st and 2pc and computes polar coordinates ('m')
 #--might be useful for displaying results of dense icon plots?
 #--'res' is an 2-column matrix with phi and dist as columns
 res=NULL 
 
 #--first compute Euclidean distance from origin...
 edToOrigin<-sqrt(apply(apply(m,2,FUN=function(x){x^2}),1,sum))
 res<-cbind(phi=(atan2(m[,2],m[,1])*180/pi)+180,dist=edToOrigin)  
 return(res)
}


#--convert above objects to polar...
#xsannotated1.log10.dc.svd.pc12.polar<-pca.convert.to.polar(xsannotated1.log10.dc.svd.pc12)
#xsannotated1.log10.dc.svd.pc23.polar<-pca.convert.to.polar(xsannotated1.log10.dc.svd.pc23)
#xsannotated1.log10.ave.within.tp.dc.svd.pc12.polar<-pca.convert.to.polar(xsannotated1.log10.ave.within.tp.dc.svd.pc12)
#xsannotated1.log10.ave.within.tp.dc.svd.pc23.polar<-pca.convert.to.polar(xsannotated1.log10.ave.within.tp.dc.svd.pc23)

#--do some plots...
#plot(xsannotated1.log10.dc.svd.pc12.polar,pch=16,col=1,las=1,main="Double centered, all mass features",xlab="Angle to +ve x-axis (phi)",ylab="Distance to origin")
#--and highlight selected mass features...
#plot(xsannotated1.log10.dc.svd.pc12.polar,pch=16,col="lightgrey",las=1,main="Double centered, all mass features",xlab="Angle to +ve x-axis (phi)",ylab="Distance to origin")
#points(xsannotated1.log10.dc.svd.pc12.polar[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),],pch=16,col=1)

#--try some experiments with cross-correlation distance...only use 20 sample (within-time-point averaged) case for this analysis.

#--construct a vector to use as a "bait"
#--in this case, I'll use a step response to model DO across time (it turns "on" at time point 13/20)
#sim.oxygen.profile<-rep(0,6)
#sim.oxygen.profile[4:6]<-1

#--for each row of 'xsannotated1.log10' and/or 'xsannotated1.log10.dc', compute cross-correlation with this profile, and extract max(abs(ccf()$acf)) and ccf()$lag[which.max(abs(ccf()$acf))]
compute.ccf.against.bait.signal<-function(bait,mat)
{
 #--'bait' is a numeric vector of length M
 #--'mat' is a data matrix of dimension N-by-M
 #--compute ccf between 'bait' and each row of 'mat'
 #--return the following vectors as list elements: 'ccf.max'=max(abs(ccf()$acf)) and 'ccf.maxlag'=ccf()$lag[which.max(abs(ccf()$acf))]
 res<-list(ccf.max=NULL,ccf.maxlag=NULL)
 
 for(i in (1:nrow(mat)))
 {
  curccf<-ccf(bait,mat[i,],plot=F)
  mind<-which.max(abs(curccf$acf))
  res$ccf.max<-c(res$ccf.max,abs(curccf$acf)[mind])
  res$ccf.maxlag<-c(res$ccf.maxlag,curccf$lag[mind])
 } 

 return(res)
 
}

#--try some plots...
#plot((xsannotated1.log10.ave.within.tp.dc[14,]-min(xsannotated1.log10.ave.within.tp.dc[14,]))/(max(xsannotated1.log10.ave.within.tp.dc[14,])-min(xsannotated1.log10.ave.within.tp.dc[14,])),type="l",las=1,xlab="Time",ylab="Normalised signal")
#legend(1,1,c("Mass feature 14","Ideal DO"),ncol=1,lty=c(1,2),bty="n")
#points(sim.oxygen.profile,type="s",lty=2)
#dev.print(pdf,file="cross-cor-example.pdf")

#do.ccf.xsannotated1.log10.ave.within.tp<-compute.ccf.against.bait.signal(sim.oxygen.profile,xsannotated1.log10.ave.within.tp)
#do.ccf.xsannotated1.log10.ave.within.tp.dc<-compute.ccf.against.bait.signal(sim.oxygen.profile,xsannotated1.log10.ave.within.tp.dc)

#--how much of this is non-random?
#--permute columns of 'xsannotated1.log10' etc and repeat a few times...
permute.columns<-function(x){x[,sample(1:ncol(x),ncol(x),F)]}

perform.randon.ccf.analysis<-function(bait,m,B)
{
 res<-list(r.ccf.max=NULL,r.ccf.maxlag=NULL)
 for(i in (1:B))
 {
  print(i)
  mr<-permute.columns(m)
  curccfVec<-compute.ccf.against.bait.signal(bait=bait,mat=mr)
  res$r.ccf.max<-cbind(res$r.ccf.max,curccfVec$ccf.max)
  res$r.ccf.maxlag<-cbind(res$r.ccf.maxlag,curccfVec$ccf.maxlag)
}

return(res)

}

#--compute randomised versions...
#random100.do.ccf.xsannotated1.log10.ave.within.tp<-perform.randon.ccf.analysis(sim.oxygen.profile,xsannotated1.log10.ave.within.tp,100)
# shiv changes from do.ccf.xsannotated1.log10.ave.within.tp to xsannotated1.log10.ave.within.tp
#--and ranks of observed vs. random cases...
find.rank.of.obs.vs.null<-function(obs,nullMat)
{
 res=NULL
 for(k in (1:length(obs)))
 {
  tmp<-length(which(nullMat[k,]>=obs[k]))
  res<-c(res,(tmp+1)/(ncol(nullMat)+1))
 }
 return(res)
}

#--this are probably P-values under reasonable definition?
#--so use R/qvalue to compute q-values for each mass feature
#pval.do.ccf.xsannotated1.log10.ave.within.tp<-find.rank.of.obs.vs.null(do.ccf.xsannotated1.log10.ave.within.tp$ccf.max,random100.do.ccf.xsannotated1.log10.ave.within.tp$r.ccf.max)
#qobj.do.ccf.xsannotated1.log10.ave.within.tp<-qvalue(pval.do.ccf.xsannotated1.log10.ave.within.tp)

#> qobj.do.ccf.xsannotated1.log10.ave.within.tp$pi0
#[1] 0.5313482

#--plot this up as follows...
#matplot(t(xsannotated1.log10.ave.within.tp.dc[which(qobj.do.ccf.xsannotated1.log10.ave.within.tp$qvalue<0.05),]),type="l",lty=1,col=1)

#--and do a PCA on these...('ccfsig')
#xsannotated1.log10.ave.within.tp.dc.ccfsig.svd<-svd(xsannotated1.log10.ave.within.tp.dc[which(qobj.do.ccf.xsannotated1.log10.ave.within.tp$qvalue<0.05),])
#> str(xsannotated1.log10.ave.within.tp.dc.ccfsig.svd)
#List of 3
# $ d: num [1:20] 16.9 9.89 4.11 3.04 2.92 ...
# $ u: num [1:506, 1:20] -0.01067 -0.02972 0.00235 -0.03269 -0.01902 ...
# $ v: num [1:20, 1:20] 0.02992 -0.00254 -0.01006 -0.04877 -0.06364 ...

#xsannotated1.log10.ave.within.tp.dc.ccfsig.svd.pc12<-cbind(xsannotated1.log10.ave.within.tp.dc.ccfsig.svd$u[,1]*xsannotated1.log10.ave.within.tp.dc.ccfsig.svd$d[1],xsannotated1.log10.ave.within.tp.dc.ccfsig.svd$u[,2]*xsannotated1.log10.ave.within.tp.dc.ccfsig.svd$d[2])

#--and an icon plot...
#GE.Icon(xsannotated1.log10.ave.within.tp.dc.ccfsig.svd.pc12,"Double centered, all mass features",xsannotated1.log10.ave.within.tp.dc[which(qobj.do.ccf.xsannotated1.log10.ave.within.tp$qvalue<0.05),],xlab="PC1",ylab="PC2",las=1)

#--how many of these are in our >5 set?
#intersect(which(qobj.do.ccf.xsannotated1.log10.ave.within.tp$qvalue<0.05),which(xsannotated1.log10.row.inter.to.intra.tp.var>5))
#> str(intersect(which(qobj.do.ccf.xsannotated1.log10.ave.within.tp$qvalue<0.05),which(xsannotated1.log10.row.inter.to.intra.tp.var>5)))
# int [1:149] 19 28 82 128 150 157 173 177 206 223 ...

#--plots for abstract...
#par(mfrow=c(2,2))
#--plot 1...
#plot(xsannotated[,"mz"],xsannotated[,"rt"],xlab="m/z",ylab="rt",las=1,pch=16,col="lightgrey")
#mtext("A",side=3,at=0,line=1,cex=1.2)
#points(xsannotated[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),"mz"],xsannotated[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),"rt"],pch=16,col=1)
#--plot 2...
#GE.Icon(xsannotated1.log10.ave.within.tp.dc.with.var.criteria.svd.pc12,"",xsannotated1.log10.ave.within.tp.dc[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),],xlab="PC1",ylab="PC2",las=1)
#abline(h=0,lty=2)
#abline(v=0,lty=2)
#mtext("B",side=3,at=-9.5,line=1,cex=1.2)
#--plot 3...
#plot(10^matCorExample[,1],matCorExample[,2],type="b",pch=3,log="x",las=1,xlab="Between:within time point variability",ylab="Mantel correlation statistic")
#mtext("C",side=3,at=0.2,line=1,cex=1.2)
#--plot 4...
#plot(do.ccf.xsannotated1.log10.ave.within.tp$ccf.max,qobj.do.ccf.xsannotated1.log10.ave.within.tp$qval,pch=3,log="y",xlim=c(0,1),las=1,xlab="Cross-correlation metric",ylab="False discovery rate (q-value)")
#abline(h=0.05,lty=2)
#mtext("D",side=3,at=-0.15,line=1,cex=1.2)
#par(mfrow=c(1,1))
#dev.print(pdf,file="kirkegaard-abstract.pdf")

#--write a function to use polar coordinate PCA plots to better advantage...
#--plot 3 panels, the first is the augmented phi-vs-r plot, showing the region selected.
#--the second is the r/t-m/z plane, with the relevant points highlighted
#--the third is a matplot of the metabolite profiles

augmented.pca.analysis<-function(massFeatureDataMat,massFeatureCountsMat,polarPCA,massFeatureSetOfInterest,phi1,phi2,rmin=0)
{
 #--see notes above
 #--'massFeatureDataMat' is a 2-column matrix containing r/t and m/z values for all included features
 #--'massFeatureCountsMat' is a data matrix of mass features by samples
 #--'polarPCA' is a 2 column matrix containing projections of the data onto principal component axes, expressed in polar coordinates (phi,r)
 #--'massFeatureSetOfInterest' is an index vector of mass features you're particularly interested in having a look at
 #--'phi1' and 'phi2' are angle that fix the interval to display in panels 2 and 3

 res=NULL
 
 #--first check that there are some features of interest in the interval phi1-phi2
 indToUse1<-massFeatureSetOfInterest[intersect(which(polarPCA[massFeatureSetOfInterest,1]>phi1),which(polarPCA[massFeatureSetOfInterest,2]>rmin))]
 indToUse2<-massFeatureSetOfInterest[intersect(which(polarPCA[massFeatureSetOfInterest,1]<phi2),which(polarPCA[massFeatureSetOfInterest,2]>rmin))]
 indToUse<-intersect(indToUse1,indToUse2)
 if(length(indToUse)==0)
 {
  stop("No features of interest in this phi interval. Please review\n")
 }
 
 par(mfrow=c(3,1))
 #--plot the first panel...
 plot(polarPCA,pch=16,col="lightgrey",xlab="phi",ylab="r",las=1,main=paste("phi=",phi1,"-",phi2,"| r >",rmin))
 points(polarPCA[massFeatureSetOfInterest,],pch=16,col=1)
 points(polarPCA[indToUse,],pch=16,col=2)
 abline(h=rmin,lty=2)
 #--and the second...
 plot(massFeatureDataMat,pch=16,col="lightgrey",xlab="m/z",ylab="rt",las=1)
 points(massFeatureDataMat[massFeatureSetOfInterest,],pch=16,col=1)
 points(massFeatureDataMat[indToUse,],pch=16,col=2)
 #--and the third...
 matplot(t(massFeatureCountsMat[indToUse,]),col=1,type="l",lty=1,las=1)
 abline(v=13,lty=2)
 
 res<-indToUse
 return(res)
} 
   
#--we'll need the following object...
#xsannotated.rt.and.mz.mat<-xsannotated[,c("mz","rt")]
#--example...returns indices for further analysis
#> augmented.pca.analysis(xsannotated.rt.and.mz.mat,xsannotated1.log10.ave.within.tp.dc,xsannotated1.log10.ave.within.tp.dc.svd.pc12.polar,which(xsannotated1.log10.row.inter.to.intra.tp.var>5),35,45,rmin=0)
#[1]  418  848  874  912  925  941 1822 2153

#--added on 9 April 2013--
#--output table for AStream analysis...make 2 versions, one with technical replicates averages, one with them in-place.
#--see http://www.urr.cat/AStream/AStream.html for notes:
#astreamData1<-data.frame(xsannotated[,"mz"],xsannotated[,"rt"],xsannotated1)
#colnames(astreamData1)<-c("mz","rt",paste("t",paste(labels.for.tr,rep(1:3,6),sep="."),sep=""))
#astreamSamples1<-paste("t",paste(labels.for.tr),sep="")
#astreamList1<-list(data=astreamData1,class=astreamSamples1)

#--these are the within time point averaged data...
#astreamData2<-data.frame(xsannotated[,"mz"],xsannotated[,"rt"],xsannotated1.ave.within.tp)
#colnames(astreamData2)<-c("mz","rt",paste("t",1:6,sep=""))
#astreamSamples2<-rep("t",6)
#astreamList2<-list(data=astreamData2,class=data.frame(class=astreamSamples2))

#--run datanorm with default parameters...
#astream2.datanorm<-data.norm(astreamList2)

#> str(astream2.datanorm)
#List of 7
# $ data       : num [1:2260, 1:20] 48234 233123 420153 93445 518451 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : NULL
#  .. ..$ : chr [1:20] "t1" "t2" "t3" "t4" ...
# $ class      : chr [1:20] "t1" "t2" "t3" "t4" ...
# $ mz         : num [1:2260] 100 101 101 101 101 ...
# $ rt         : num [1:2260] 25.7 938 1005.9 1006.8 46.5 ...
# $ feature.set: num [1:35771, 1:6] 2 2 2 2 2 2 2 2 2 2 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : NULL
#  .. ..$ : chr [1:6] "Index feature 1" "Index feature 2" "Correlation" "%samples" ...
# $ outliers   : int(0) 
# $ scores     :List of 2
#  ..$ int_0: num [1:20] -0.536 -0.893 -0.179 -0.893 0.536 ...
#  ..$ int_1: num [1:20] 0.1341 0.0103 0.0407 0.2241 0.3098 ...

#--and 'head' the feature.set element...
#> head(astream2.datanorm$feature.set)
#     Index feature 1 Index feature 2 Correlation %samples RTdiff  M/Zdiff
#[1,]               2             638   0.8188348      100 0.9430 122.5536
#[2,]               2             649   0.8526627      100 0.8070 123.5532
#[3,]               2             659   0.7682488      100 0.8040 124.5328
#[4,]               2             661   0.8143427      100 0.9880 124.5504
#[5,]               2             794   0.8215341      100 0.8010 147.5488
#[6,]               2             977   0.8213097      100 0.8015 180.5409

#--then run isotope search...
#astream2.datanorm.isotopesearch<-isotope.search(astream2.datanorm,mz.tol=3e-3)
#ISOTOPE ANALYSIS:
#	- M/Z tolerance:	0.003u
#	- C12-isotope m/z:	12u
#	- C13#1-isotope m/z:	13.00335u
#	- C13#2-isotope m/z:	14.0067u
#	- 13C#1/12C-isotope:	143 feature pairs
#	- 13C#2/12C-isotope:	3 feature pairs
#	- Total number:	146 isotope patterns

#> dim(astream2.datanorm.isotopesearch$isotopes)
#[1] 146   6

#--default adduct search...
#astream2.datanorm.isotopesearch.adductsearch<-adduct.search(astream2.datanorm.isotopesearch)

#ADDUCT ANALYSIS:
#	- M/Z tolerance:	+/-0.003u
#	- 0 [M-CO]adducts.
#	- 0 [M-H2O]adducts.
#	- 0 [M-NH3]adducts.
#	- 11 [M+NH4]adducts.
#	- 15 [M+Na]adducts.
#	- 0 [M+H+Na]adducts.
#	- 12 [M+K]adducts.
#	- 0 [M+H+K]adducts.
#	- 1 [M+NaCOOH]adducts.
#	- 12 [M+NaHCOOH]adducts.

#> printresults("test.txt",astream2.datanorm.isotopesearch,astream2.datanorm.isotopesearch.adductsearch,mzsort=T)

#--print these out at two different ranges...
#printresults("test.txt",astream2.datanorm.isotopesearch,astream2.datanorm.isotopesearch.adductsearch,mzsort=T)
#printresults("test1.txt",astream2.datanorm.isotopesearch,astream2.datanorm.isotopesearch.adductsearch,mzsort=T,range=0.001)
#--now see manual analysis of test.txt -> test.xlsx and test1.txt -> test1.xlsx
#--read back in annotated version of test1.xlsx
#test1.annotated<-read.table(file="test1-annotated.txt",header=T,sep="\t")
#> dim(test1.annotated)
#[1] 33 11

#--make a per-compund version of these data for further annotation in Metlin...
#--for each detected compound, list in a row of a table...
#test1.annotatedTable=NULL
#for(k in (1:nrow(test1.annotated)))
#{
# cat(k)
# curmf<-test1.annotated[k,"IndexFeature"]
# curmass<-as.character(test1.annotated[k,"MassSubmitted"])
# cururl<-test1.annotated[k,"Metlin"]
# curMetlinMolId<-as.character(test1.annotated[k,"Metlin.1"])
# tmpchar<-strsplit(curMetlinMolId,"|",fixed=T)[[1]]
# tmp<-data.frame(mass.feature=rep(curmf,length(tmpchar)),mass.submitted=rep(curmass,length(tmpchar)),url=rep(cururl,length(tmpchar)),metlin=tmpchar)
# test1.annotatedTable<-rbind(test1.annotatedTable,tmp)
#} 
 
#--write this file...
#write.table(test1.annotatedTable,file="test1.annotatedTable.txt",col.names=T,row.names=F,sep="\t")
 
#--have annotated this with data from individual metabolite pages....
#--read in the mass features that have KEGG CPD ids (we'll need to annotate adducts into this...)
#test1.annotatedTable.kegg.cpd<-read.table(file="test1.annotatedTable-kegg-cpd.txt",header=T,sep="\t")
#save(test1.annotatedTable.kegg.cpd,file="test1.annotatedTable.kegg.cpd.rda")
#--use this in ~/Documents/projects/scelse/ulu.pandan.0313/connections.analysis.pilot to analyse occurence in batch1 samples...

#--how many of these features are in our selected sets?
#test1.annotatedTable.unique.mf<-unique(test1.annotatedTable[,1])

#--check these against both sets of features...
#> intersect(test1.annotatedTable.unique.mf,which(qobj.do.ccf.xsannotated1.log10.ave.within.tp$qvalue<0.05))
#[1]   19 1644
#> intersect(test1.annotatedTable.unique.mf,which(xsannotated1.log10.row.inter.to.intra.tp.var>5))
# [1]   19   66  119  293  541  542  581  667  741  925  966 1147 1546 1973

#--Rohan started here again on May 3 2013--
#--do a mass-by-rt plot...to illustrate the location of these...
#plot(xsannotated[,"mz"],xsannotated[,"rt"],xlab="m/z",ylab="rt",las=1,pch=16,col="lightgrey")
#text(xsannotated[test1.annotatedTable.unique.mf,"mz"],xsannotated[test1.annotatedTable.unique.mf,"rt"],as.character(test1.annotatedTable.unique.mf),cex=0.5)
#dev.print(pdf,file="annotated-compunds-mass-feature-plot.pdf")

#plot(xsannotated[,"mz"],xsannotated[,"rt"],xlab="m/z",ylab="rt",las=1,pch=16,col="lightgrey")
#points(xsannotated[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),"mz"],xsannotated[which(xsannotated1.log10.row.inter.to.intra.tp.var>5),"rt"],col=1,pch=16)

#plot(xsannotated[,"mz"],xsannotated[,"rt"],xlab="m/z",ylab="rt",las=1,pch=16,col="lightgrey")
#points(xsannotated[which(qobj.do.ccf.xsannotated1.log10.ave.within.tp$qvalue<0.05),"mz"],xsannotated[which(qobj.do.ccf.xsannotated1.log10.ave.within.tp$qvalue<0.05),"rt"],col=1,pch=16)

#--do some plots...
matplot(t(xsannotated1.log10.ave.within.tp.dc[intersect(test1.annotatedTable.unique.mf,which(qobj.do.ccf.xsannotated1.log10.ave.within.tp$qvalue<0.05)),]),type="l",lty=c(1,2),las=1,col=1,
            ylab="Normalised mass feature intensity",xlab="Time point",xlim=c(0,21))
abline(v=13,lty=4,col="lightgrey")
text(c(20,20),t(xsannotated1.log10.ave.within.tp.dc[intersect(test1.annotatedTable.unique.mf,which(qobj.do.ccf.xsannotated1.log10.ave.within.tp$qvalue<0.05)),])[20,],c(19,1644),pos=4,cex=0.75)

#--what do the "good variance" features look like?
#matplot(t(xsannotated1.log10.ave.within.tp.dc[intersect(test1.annotatedTable.unique.mf,which(xsannotated1.log10.row.inter.to.intra.tp.var>5)),]),type="l",lty=1,las=1,col=1)
#abline(v=13,lty=2)

#--try on a PC plot...
#plot(xsannotated1.log10.ave.within.tp.dc.svd.pc12,main="Double centered, all mass features",xlab="PC1",ylab="PC2",las=1,col="lightgrey",pch=16)
#points(xsannotated1.log10.ave.within.tp.dc.svd.pc12[which(qobj.do.ccf.xsannotated1.log10.ave.within.tp$qvalue<0.03),],pch=16,col=1)
#abline(h=0,lty=2)
#abline(v=0,lty=2)

#--Rohan finished editing here on May 3 2013--

