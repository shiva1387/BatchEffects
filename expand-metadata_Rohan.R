#--generate a binary matrix which maps how the levels of these factors map onto samples....as a prelude to statistical analysis
expand.metadata<-function(metadataMat)
{
 #--'metdataMat' is a sample-by-variable matrix of metadata: generally categorical variables
 #--expand the levels within each variable, and construct a submatrix which shows distribution of level across samples.
 #--repeat for all variable, building up a (levels in variable)-by-sample matrix
 
 res=NULL
 
 for(variable in (1:ncol(metadataMat)))
 {
  curvariable<-metadataMat[,variable]
  curlevels<-unique(as.character(curvariable))
  curvariablename<-colnames(metadataMat)[variable]
  print(curlevels)
  curmat=NULL
  for(level in (1:length(curlevels)))
  {
   curcol<-rep(0,nrow(metadataMat))
   ind<-which(curvariable==curlevels[level])
   curcol[ind]<-1
   curmat<-rbind(curmat,curcol)
  }
  colnames(curmat)<-rownames(metadataMat)
  rownames(curmat)<-paste(curvariablename,curlevels,sep=":")
  res<-rbind(res,curmat)
 }
 return(res)
}

#--in this example, 'metadataMat1' is e.g.
#> head(metadata1[,c(2,6:9)])
#              Day Day.of.extraction Volume.reduction Day.of.volume.reduction
#II_T2S_I1_A_1   1           5/10/13          0.1 ml                  5/10/13
#II_T2S_I1_A_2   1           5/10/13          0.1 ml                  5/10/13
#II_T2S_I1_B_1   1           5/10/13          0.1 ml                  5/10/13
#II_T2S_I1_B_2   1           5/10/13          0.1 ml                  5/10/13
#II_T2S_I2_A_1   1           5/10/13          0.1 ml                  5/10/13
#II_T2S_I2_A_2   1           5/10/13          0.1 ml                  5/10/13
#              Data.Acq..Date.Time
#II_T2S_I1_A_1             5/14/13
#II_T2S_I1_A_2             5/14/13
#II_T2S_I1_B_1             5/14/13
#II_T2S_I1_B_2             5/14/13
#II_T2S_I2_A_1             5/14/13
#II_T2S_I2_A_2             5/14/13

a<-expand.metadata(metadata1[,c(2,6:9)])
#--offset rows of "a" with integers show we can plot easily...
construct.layer.rug<-function(bm)
{
 #--'bm' is a binary matrix you want to plot in a "layered rug" format (e.g. to see group memberships etc)
 res=NULL
 offsetMat<-matrix(rep(1:nrow(bm),ncol(bm)),nrow(bm),ncol(bm))
 res<-bm*offsetMat
 res[res<1e-15]<-NA
 rownames(res)<-rownames(bm)
 colnames(res)<-colnames(bm)
 return(res)
}

aclr<-construct.layer.rug(a)

#--combine with boxplot etc
#> par()$mar
#[1] 5.1 4.1 4.1 2.1

par(mfrow=c(2,1))
par(mar=c(1.1,9.1,4.1,2.1))
boxplot(log10(mfdata1+1),xaxt="n",las=1,ylab="Mass feature intensity (log10+1)")
par(mar=c(5.1,9.1,0.1,2.1))
matplot(t(aclr),type="p",pch=15,col=1,yaxt="n",ylab="",xlab="Samples")
axis(side=2,at=1:nrow(aclr),labels=rownames(aclr),las=1,cex.axis=0.65)
par(mfrow=c(1,1)) 

