###########################################
# Data Analysis in R-Malaysian algae data #
###########################################
# Author(s): Shiv
# Version: 28092013 
# Input: ".tsv" file from XCMS 
# Software: XCMS
# File Location :
# Output: 
# Modified By :Shivshankar Umashankar 

#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 

library(ggplot2)
library(data.table)
library(xcms)
library(stringr)
library(vegan)
library(ape)
library(gplots)


#############
# User      #
# Specific  #
# Variables #
#############

directory<- "F:/Vinay's Algae Data/Aug2013/RawData/MZXML/15th MZxml"
#Contains path for .tsv file
# The "PATH", Remember to use "/" and not "/" in the path

num_cores<-8; # The number of CPUs available, be aware of memory issue

setwd(directory)
getwd()

#################
# Get filenames #
#################


#mzxmlfiles<-list.files(directory,pattern=".mzXML",full.names=TRUE) # Get filenames with a given pattern and include subdirectories 

####### reading in the mass spec data 

mzfilename<-"Algae_22strains_indreps_blanks_matrix.tsv"
ms_data_total<-read.table(mzfilename,sep="\t",header=T,check.names=FALSE,row.names=1)
str(ms_data_total)

################Formatting column names

names(ms_data_total)<-gsub('\\[', '', names(ms_data_total))
names(ms_data_total)<-gsub('\\]', '', names(ms_data_total))
a<-names(ms_data_total)
b<-gsub('\\(raw)', '', a)
c<-gsub('\\, ','\\-', b)
names(ms_data_total)<-c
rm(a,b,c)

#Only for full data -22 strains
#ms_data_total_strains<-ms_data_total[, -grep('ACN', names(ms_data_total))] #removing blanks
#ms_data_total_strains<-ms_data_total_strains[, -grep('matri_', names(ms_data_total_strains))] #removing matrix

metadata<-read.table("Algae_22strains_indreps_blanks_matrix_metadata.txt",check.names=FALSE,sep="\t")

#ms_data<-ms_data_total_strains[,1:262]
#a<-names(ms_data_total_strains)

# Assigning sample groups

a<-names(ms_data_total)
#strsplit(a[3],"_")
#gsub("^([^_]*_[^_]*)_.*$", "\\1", a[3])
#b<-paste(strsplit(a,"_")[[1]][1:2],collapse = "_")

SampleGroups<-sapply(a, function(x) paste(strsplit(x,"_")[[1]][1:2],collapse = "_"))
SampleGroups<-as.vector(SampleGroups)

SampleGroups<-gsub('b1','',SampleGroups)
SampleGroups<-gsub('b2','',SampleGroups)
SampleGroups<-gsub('b3','',SampleGroups)

#write.table(SampleGroups0,"SampleGroups.txt",row.name=FALSE,quote=F)
#SampleGroups1<-read.table("SampleGroups.txt",header=F)

#no_of_info_cols<-1
#data_start_index<-no_of_info_cols+1;
#data_end_index<-dim(ms_data_total)[2];

ms_data<-ms_data_total[,1:311]
names(ms_data)
dim(ms_data)
#ms_data<-log(ms_data)

#ms_data_day4<-ms_data[, grep('D4', names(ms_data))] 
#ms_data_day12<-ms_data[, grep('D12', names(ms_data))] 

SampleGroup<-SampleGroups[1:311]

#name_list <- strsplit(SampleGroup, "_")
#GrowthStage<-sapply(name_list , function (x) if(length(x) == 2) x[1] else as.character(NA))
#StrainId<-sapply(name_list , function (x) if(length(x) == 2) x[2] else as.character(NA))
#StrainName<-as.factor(names(ms_data))
#StrainName0<-names(ms_data)

#######################################
###### Exploratory data analysis ######
#######################################

##########################
# Histogram of mz and rt #
##########################

mz_rt <- strsplit(rownames(ms_data_total), "\\@")
mz<-sapply(mz_rt , function (x) if(length(x) == 2) x[1] else as.character(NA))
rt<-sapply(mz_rt , function (x) if(length(x) == 2) x[2] else as.character(NA))

png("algae_417reps_retcor_MPP_sorted.png")
plot(rt,mz,pch=1,main="algae_417reps_retcor_MPP_sorted",ylab="m/z",xlab="rt")
dev.off()

### Measure of variation ########

png('boxplot.png',width=8000)
boxplot(log(ms_data[227:232]),range=0,xlab="Samples",ylab="Log2 feature counts",las=1,cex.axis=0.7)
dev.off()

## Generating boxplot and cv for each sample group

groups<-unique(SampleGroup)

for(i in 1:length(groups))
{
png(paste(groups[i],".png"),width=500,height=500)
ms_data_grp1<-ms_data[, grep(groups[i], names(ms_data))]
par(mar=c(10,2,2,2))
boxplot(log(ms_data_grp1),range=0,ylab="Log2 feature counts",las=2,cex.axis=1)
if(ncol(ms_data_grp1)==6)
{
text(3.5, 10,font=3, paste("cv:", round(cv(ms_data_grp1[,1]), 2), "   cv:", round(cv(ms_data_grp1[,2]), 2),"   cv:", round(cv(ms_data_grp1[,3]), 2),
                  "   cv:", round(cv(ms_data_grp1[,4]), 2),"   cv:", round(cv(ms_data_grp1[,5]), 2),"   cv:", round(cv(ms_data_grp1[,6]), 2)
                  ),cex=1)
}
else if(ncol(ms_data_grp1)==4)
{
  text(3.5, 10,font=3, paste("cv:", round(cv(ms_data_grp1[,1]), 2), "   cv:", round(cv(ms_data_grp1[,2]), 2),"   cv:", round(cv(ms_data_grp1[,3]), 2),
                             "   cv:", round(cv(ms_data_grp1[,4]), 2)
  ),cex=1)
}
else if(ncol(ms_data_grp1)==7)
{
  text(3.5, 10,font=3, paste("cv:", round(cv(ms_data_grp1[,1]), 2), "   cv:", round(cv(ms_data_grp1[,2]), 2),"   cv:", round(cv(ms_data_grp1[,3]), 2),
                             "   cv:", round(cv(ms_data_grp1[,4]), 2),"   cv:", round(cv(ms_data_grp1[,5]), 2),"   cv:", round(cv(ms_data_grp1[,6]), 2),
                             "   cv:", round(cv(ms_data_grp1[,7]), 2)
  ),cex=1)
}
dev.off()
graphics.off()
font=1
}

#plot(t(log1p(ms_data))~SampleGroup,range=0,xlab="Samples",ylab="Log2 feature counts",las=1,col=SampleGroup)

# plotting coefficient of variations among samples(cols)

#Coefficient of variation
png("cv_samples.png",width=4000)
plot(1:ncol(ms_data),apply(ms_data,2,cv)) # For each column 
dev.off()

# for each feature in each sample group

ms_data_tst<-data.table(t(log(ms_data)))
ms_data_tst1<-cbind(ms_data_tst,as.factor(SampleGroup))
lapply(ms_data_tst1,class) #checking col classes
cv1 <- function(x) (sd(x)/mean(x)) * 100

cv_mz_feature2<-t(ms_data_tst1[,lapply(.SD,cv1),by=V2])
#lapply(cv_mz_feature,class)
cv_mz_feature<-as.data.frame(cv_mz_feature2[2:nrow(cv_mz_feature2),])
#cv_mz_feature<-as.data.frame(lapply(cv_mz_feature,as.numeric))

colnames(cv_mz_feature)<-unique(SampleGroup)

write.table(cv_mz_feature,"mz_features_cv.txt",quote=FALSE,sep="\t")
mz_grp_cv<-read.table("mz_features_cv.txt",sep='\t',header=TRUE,row.name=1)
mz_grp_cv[is.na(mz_grp_cv)]<-0


####

internalStandard<-mz_grp_cv[4473,]

#mz_grp_cv<-as.data.frame(mz_grp_cv)
#colnames(mz_grp_cv)<-d

mz_grp_cv<-data.frame(cbind(mz,rt,mz_grp_cv))
for (i in seq(3,ncol(mz_grp_cv),4))
{
  mz<-as.vector(mz_grp_cv$mz)
  rt<-as.vector(mz_grp_cv$rt)
  rt<-as.numeric(rt)
  print (i)
  y<-mz_grp_cv[,i]
  rt_cv<-cbind(as.character(rt),as.numeric(y))
  rt_cv<-data.frame(rt,mz_grp_cv[,i],mz_grp_cv[,i+1],mz_grp_cv[,i+2],mz_grp_cv[,i+3])
  rt_cv<-rt_cv[order(rt_cv[,1]),]
  
  ## now plot scatter plot
  png(paste("hist-coefofvar/", i, "_mz.png", sep = ""),height=800,width=1200)
  par(mfrow=c(2,2))
  plot(mz,mz_grp_cv[,i],main=names(mz_grp_cv)[i],xlab="mz",ylab="cv of log(counts)")
  plot(mz,mz_grp_cv[,i+1],main=names(mz_grp_cv)[i+1],xlab="mz",ylab="cv of log(counts)")
  plot(mz,mz_grp_cv[,i+2],main=names(mz_grp_cv)[i+2],xlab="mz",ylab="cv of log(counts)")
  plot(mz,mz_grp_cv[,i+3],main=names(mz_grp_cv)[i+3],xlab="mz",ylab="cv of log(counts)")
  
  ## now plot histogram
  # hist(mz_grp_cv[,i],main=names(mz_grp_cv)[i],ylab="number of features",xlab="cv of log(counts)")
  # hist(mz_grp_cv[,i+1],main=names(mz_grp_cv)[i+1],ylab="number of features",xlab="cv of log(counts)")
  # hist(mz_grp_cv[,i+2],main=names(mz_grp_cv)[i+2],ylab="number of features",xlab="cv of log(counts)")
  # hist(mz_grp_cv[,i+3],main=names(mz_grp_cv)[i+3],ylab="number of features",xlab="cv of log(counts)")
  
  dev.off()
  
  png(paste("CoefficientofVariation/", i, "_rt.png", sep = ""),height=800,width=1200)
  par(mfrow=c(2,2))
  plot(rt_cv[,1],rt_cv[,2],main=names(mz_grp_cv)[i],xlab="rt",ylab="cv of log(counts)")
  plot(rt_cv[,1],rt_cv[,3],main=names(mz_grp_cv)[i+1],xlab="rt",ylab="cv of log(counts)")
  plot(rt_cv[,1],rt_cv[,4],main=names(mz_grp_cv)[i+2],xlab="rt",ylab="cv of log(counts)")
  plot(rt_cv[,1],rt_cv[,5],main=names(mz_grp_cv)[i+3],xlab="rt",ylab="cv of log(counts)")
  dev.off()
}
graphics.off()
# rm(rt_cv)

#plot(mz,as.numeric(mz_grp_cv[,4]),main=names(mz_grp_cv)[4],xlab="mz",ylab="cv of log trans counts")
#plot(mz,as.numeric(mz_grp_cv[,5]),main=names(mz_grp_cv)[5],col="green")
#plot(mz,as.numeric(mz_grp_cv[,6]),main=names(mz_grp_cv)[6],col="red")
#plot(mz,as.numeric(mz_grp_cv[,7]),main=names(mz_grp_cv)[7],col="blue")

#rownames(cv_mz_feature)<-rownames(ms_data_total)

#abline(v=seq(1,24,3)-0.5)
#axis(side=1,at=seq(2,24,3),as.character(1:8))
dev.off()


### Counting number of zeroes for each m/z feature in each group

zero_count<-function(x) sum(x == 0) 
zero_mz_feature2<-t(ms_data_tst1[,lapply(.SD,zero_count),by=V2])
zero_mz_feature<-as.data.frame(zero_mz_feature2[2:nrow(zero_mz_feature2),])
colnames(zero_mz_feature)<-unique(SampleGroup)

write.table(zero_mz_feature,"mz_features_zeroes.txt",quote=FALSE,sep="\t")
mz_grp_zero<-read.table("mz_features_zeroes.txt",sep='\t',header=TRUE,row.name=1)

mz_grp_zero<-data.frame(cbind(mz,rt,mz_grp_zero))

for (i in seq(3,ncol(mz_grp_cv),2))
{
  mz<-as.vector(mz_grp_cv$mz)
  rt<-as.vector(mz_grp_cv$rt)
 ## now plot scatter plot
  png(paste("noofZeros/", i, "_mz.png", sep = ""),height=800,width=1200)
  par(mfrow=c(2,2))
  plot(mz,mz_grp_zero[,i],main=names(mz_grp_zero)[i],xlab="mz",ylab="No of zeros")
  hist(mz_grp_zero[,i],main=names(mz_grp_zero)[i],xlab="no of zeroes",ylab="Features", breaks=6)
  plot(mz,mz_grp_zero[,i+1],main=names(mz_grp_zero)[i+1],xlab="mz",ylab="No of zeros")
  hist(mz_grp_zero[,i+1],main=names(mz_grp_zero)[i+1],xlab="no of zeroes",ylab="Features", breaks=6)
  dev.off()
}
graphics.off()

#CoefofVar-Vs-Zeros

for (i in seq(3,ncol(mz_grp_cv),2))
{
  mz<-as.vector(mz_grp_cv$mz)
  rt<-as.vector(mz_grp_cv$rt)
  ## now plot scatter plot
  png(paste("CoefofVar-Vs-Zeros/", i, ".png", sep = ""),height=800,width=1200)
  par(mfrow=c(2,2))
  
  cor1<-cor(mz_grp_zero[,i],mz_grp_cv[,i], use="all.obs", method="pearson")
  plot(mz_grp_zero[,i],mz_grp_cv[,i],main=names(mz_grp_zero)[i],xlab="No of zeroes",ylab="Coefficient of variation")
  text(1, 180, paste("R2 =", round(cor1, 3)),cex=1)
  hist(mz_grp_cv[,i],main=names(mz_grp_zero)[i],xlab="Coefficient of variation",ylab="No of Features")
  
  cor2<-cor(mz_grp_cv[,i+1],mz_grp_zero[,i+1], use="all.obs", method="pearson")  
  plot(mz_grp_zero[,i+1],mz_grp_cv[,i+1],main=names(mz_grp_zero)[i+1],xlab="No of zeroes",ylab="Coefficient of variation")
  text(1, 180, paste("R2 =", round(cor2, 3)),cex=1)
  hist(mz_grp_cv[,i+1],main=names(mz_grp_zero)[i+1],xlab="Coefficient of variation",ylab="No of Features")
  
  dev.off()
}
graphics.off()

######### Converting values based on 50% zero count per sample group

new_ms_data<-ms_data_tst1[, lapply(.SD, function(v) { 
  len <- length(v)
  if((sum(v==0)/len)>0.5) rep(0,len) else v
}), by=V2]

rownames(new_ms_data)<-colnames(ms_data)
new_ms_data<-as.data.frame(t(new_ms_data))

###### removing replicates(samples with large variance) from the dataset

rem <- c('D4_094_b1_r002','D4_104_b3_r002','D4_245b1_r001','D4_254b1_r001','D4_325_b1_r002')

new_ms_data1<-new_ms_data[, !names(new_ms_data) %in% rem]

names(new_ms_data1)



#df[, lapply(.SD, function(v) { 
#  len <- length(v)
#  if((sum(v==0)/len)>0.5) rep(0L,len) else v
#}), by="Group", .SDcols=c("r1","r2","r3")]


########################################
###### Multivariate data analysis ######
########################################

#heatmap.2(as.matrix(mz_grp_cv),dendrogram="column",trace="none",mar=c(10,10),Rowv=FALSE,labRow= NULL, ylab=NULL)

########ttest

## ttest between columns 

combos <- combn(ncol(ms_data),2)

adply(combos, 2, function(x) {
  test <- t.test(ms_data_test[, x[1]], ms_data_test[, x[2]],paired=F,var.equal=F,p.adjust="BH")

  out <- data.frame("var1" = colnames(ms_data)[x[1]]
                    , "var2" = colnames(ms_data[x[2]])
                    , "t.value" = sprintf("%.3f", test$statistic)
                    ,  "df"= test$parameter
                    ,  "p.value" = sprintf("%.3f", test$p.value)
                    )
  return(out)

})


#t.test(log(ms_data[,1]),log(ms_data[,5]),paired=F,var.equal=F,p.adjust="BH")
#pairwise_ttest_strains<-pairwise.t.test(t(ms_data_test),StrainName_test,p.adj="BH", paired=F, var.equal=T)
#pairwise.wilcox.test(t(ms_data_test),p.adj="BH", paired=F, var.equal=T)

## ttest on rows
t.list <- apply(ms_data,1,function(x){t.test(x[1:6],x[61:66],paired=TRUE,var.equal=F,p.adjust=BH)$p.value}) 
a.list <- apply(ms_data,1,function(x){wilcox.test(x[1:42],x[43:84],paired=TRUE,var.equal=F,p.adj=BH,exact=F)$p.value})
 
id.sig <- which(a.list < 0.01 );
metab.sig<-names(id.sig)
ttest.table <- cbind(ms_data[metab.sig, ])

#fold change
fc_val<-apply(ttest.table,1,function(x){foldchange(mean(x[1:42]),mean(x[43:84]))})
fc_sig<- which(fc_val>10);

fc.metab.sig<-names(fc_sig)
significant.metab <- cbind(ms_data[fc.metab.sig,],fc_val[fc_sig])

write.table(significant.metab,"comparison_highVslow_DAY12.txt",sep='\t',col.names=NA,quote=FALSE)

# Writing to a file

write.table(pairwise_ttest_strains$p.value,"paiwise_day4.txt",sep='\t',col.names=NA,quote=FALSE)
day4_strains<-read.table("paiwise_day4.txt",sep='\t',header=TRUE,row.names=1)

#apply(d, 2, function(d) {t.test(x = d[,1], y = d[,2])$p.value})

########## Clustering ######

d <- dist(t(new_ms_data1), method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
#pdf('hclust.pdf',width=800)
png('hclust_reps_zeroCor.png',width=8000)
plot (fit)
dev.off()

clusDendro<-as.dendrogram(fit,hang=2)
labelColors<-c("black","red","blue")
clusMember<-rep(1,length(names(new_ms_data1)))
clusMember[grep("D12",names(new_ms_data1))]<-2
clusMember[grep("D4",names(new_ms_data1))]<-3
names(clusMember)<-names(new_ms_data1)

colLab <- function(n)
{
  if(is.leaf(n)) {
    a <- attributes(n)
    # clusMember - a vector designating leaf grouping
    # labelColors - a vector of colors for the above grouping
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}

clusDendro<-dendrapply(clusDendro, colLab)

png('hclust_reps_blanks_matrix_rem5rep.png',width=8000,height=1200)
plot(clusDendro)
dev.off()
  
#### PCOA

bray_grps.D <- vegdist(t(new_ms_data1[,25:281]), "bray")#calculating distance matrix using bray curtis
res_grps <- pcoa(bray_grps.D)
res_grps$values
pdf(file="PCOA_samples.pdf")
biplot(res_grps)
dev.off()

pdf(file="PCOA_ordplot_zeroCor_rem5rep.pdf",width=12, height=12, paper="a4r")
pc<-cmdscale(bray_grps.D, k=10, eig=TRUE, add=TRUE, x.ret =TRUE) 
#PCoA.res<-capscale(bray_grps.D~1,distance="bray") 
#Create ordination plot     
fig<-ordiplot(scores(pc)[,c(1,2)], type="t", main="PCoA Samples",cex=0.5)
dev.off()

x11()
ordiplot (scores(pc)[,c(1,2)], display = 'sp', type = 'n',main="PCoA Samples", cex=0.5)
points(scores(pc)[,c(1,2)], col = clusMember , pch = clusMember )
legend("bottomleft", legend = c("Day4","Day12"), pch = 1:2,col = c("blue","red"))


#### PCA

fit <- princomp(new_ms_data, scale=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
screeplot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

fit <- prcomp(new_ms_data,scale=TRUE)
mycolors <- c("red", "green", "blue", "magenta", "black")
plot(fit$x, type="n"); text(fit$x, rownames(fit$x), cex=0.5)

library(psych)
fit <- principal(log1p(ms_data), nfactors=5, rotate="varimax")
fit


####nn-graphs

ms_graph_data_strains<-t(log1p(ms_data))

rownames(ms_graph_data_strains)<-rownames(ms_graph_data_strains)
ms.d_str<-as.matrix(dist(ms_graph_data_strains,"euclidean"))
ms.1nn_str<-define.1nn.graph(ms.d_str)

pdf(file="NN_graph_strains.pdf",paper="a4r")
plot(define.1nn.graph(ms.d_str),layout=layout.fruchterman.reingold,vertex.label=V(ms.1nn_str)$name,vertex.label.cex=0.8,main="Strains")
dev.off()


################PERMANOVA

#ms_data_total$DayStrain<-paste0(ms_data_total$Day,ms_data_total$Strain)
a<-adonis(t(log(ms_data+1)) ~ GrowthStage*StrainId+RunDay, strata=StrainId  , data=ms_data,method = "bray", perm=999)
b1<-metaMDS(log(ms_data+1),k=2)

plot(b)
ordiplot(b,type="n")
orditorp(b,display="species",col="red",air=0.01)
orditorp(b,display="sites",cex=0.01,air=0.01)

distance<-vegdist(log(t(ms_data+1)), method="bray")
model1<-betadisper(distance, SampleGroup)
permute_data<-permutest(model1, pairwise = TRUE,permutations=999)
posthoc_test<-TukeyHSD(model1)

write.table(permute_data$pairwise,"permute_data_day12",quote=F)
write.table(posthoc_test$group,"posthoc_test_day12",quote=F)

###############Plotting TIC's

#combined

getTICs(files=mzxmlfiles,pdfname=paste("xcms_agilent_TIC_algae_15th.pdf", sep = ""),rt="raw")

#ws

mzxmlfiles1<-mzxmlfiles[1:6] 
getTICs(files=mzxmlfiles1,pdfname=paste("15th_ACN_Blank.pdf", sep = ""),rt="raw")

#tt8
tt8_mzxmlfiles<-mzxmlfiles[-grep('WS',mzxmlfiles)] 
getTICs(files=tt8_mzxmlfiles,pdfname=paste("xcms_agilent_TIC_TT8_140513.pdf", sep = ""),rt="raw")

######### ggplot ########

ggplot(ms_data_nd_new, aes_string(x=rt,y=n)) + geom_line(aes(colour = series)) + facet_grid(series ~ .)+theme(axis.ticks = element_blank(), axis.text.y = element_blank())
lapply(list("00111","12922"), function(i) ggplot(i,aes(x=rt,value))+geom_point())
log_blanks<-log1p(blanks)
log_blanks$rt<-ms_data_total[,5]
log_blanks<- melt(log_blanks,  id = 'rt', variable_name = 'series')
ggplot(log_blanks, aes(x=rt,value)) + geom_line(aes(colour = series)) +geom_rug(sides="b")

######### TIC overlays #############################################

getTIC <- function(file,rtcor=NULL) {
     object <- xcmsRaw(file)
     cbind(if (is.null(rtcor)) object@scantime else rtcor, rawEIC(object,mzrange=range(object@env$mz))$intensity) 
}

##
##  overlay TIC from all files in current folder or from xcmsSet, create pdf
##
getTICs <- function(xcmsSet=NULL,files=NULL, pdfname="TICs.pdf",rt=c("raw","corrected")) {
  if (is.null(xcmsSet)) {
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                      "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
    if (is.null(files))
        files <- getwd()
    info <- file.info(files)
    listed <- list.files(files[info$isdir], pattern = filepattern,
                          recursive = TRUE, full.names = TRUE)
    files <- c(files[!info$isdir], listed)
  } else {
    files <- filepaths(xcmsSet)
  }

  N <- length(files)
  TIC <- vector("list",N)

  for (i in 1:N) {
      cat(files[i],"\n")
      if (!is.null(xcmsSet) && rt == "corrected")
        rtcor <- xcmsSet@rt$corrected[[i]] else 
          rtcor <- NULL
      TIC[[i]] <- getTIC(files[i],rtcor=rtcor)
  }

  pdf(pdfname,w=16,h=10)
      cols <- rainbow(N)
      lty = 1:N
      pch = 1:N
      xlim = range(sapply(TIC, function(x) range(x[,1])))
      ylim = range(sapply(TIC, function(x) range(x[,2])))
      plot(0, 0, type="n", xlim = xlim, ylim = ylim, main = "Total Ion Chromatograms", xlab = "Retention Time", ylab = "TIC")
      for (i in 1:N) {
      tic <- TIC[[i]]
      points(tic[,1], tic[,2], col = cols[i], pch = pch[i], type="l")
      }
      legend("topright",paste(basename(files)), col = cols, lty = lty, pch = pch)
  dev.off()
  
  invisible(TIC)
}

######### Averaging samples # ref http://stackoverflow.com/questions/10704344/averaging-every-16-columns-in-r

# Average function 
byapply <- function(x, by, fun, ...)
{
    # Create index list
    if (length(by) == 1)
    {
        nc <- ncol(x)
        split.index <- rep(1:ceiling(nc / by), each = by, length.out = nc)
    } else # 'by' is a vector of groups
    {
        nc <- length(by)
        split.index <- by
    }
    index.list <- split(seq(from = 1, to = nc), split.index)

    # Pass index list to fun using sapply() and return object
    sapply(index.list, function(i)
            {
                do.call(fun, list(x[, i], ...))
            })
}


################################## Nearest neighbours algorithm #########################
# Author(s): Shiv
# Version :19092013

#--code for constructing and analysing nearest neighbour graphs
#--started by RW on 9 Decvember 2012--
#--in this instance, we're restricting attention to 1-nn graphs

#--take a (named) distance matrix, for each element, find it's nearest neighbour, construct a graph based on non-redundant edges
define.1nn.graph<-function(Dmat)
{
 #--take a (named) distance matrix, Dmat, for each element, find it's nearest neighbour, construct a graph based on non-redundant edges
 #--Dmat should be class 'matrix' not 'dist'
 #--returns a graph
 
 res=NULL

 #--define nodes and targets... 
 nodes<-colnames(Dmat)
 nnodes<-rep("",nrow(Dmat))
 val<-rep("",nrow(Dmat))

 #--remove diagonal...
 ndDmat<-Dmat
 diag(ndDmat)<-NA

 #--go through set of nodes and find nearest neighbour...
 for(i in (1:length(nodes)))
 {
  curcol<-ndDmat[,nodes[i]]
  nnodes[i]<-nodes[which.min(curcol)]
  val[i]<-sort(curcol, FALSE)[1]
 }

 #--make edgelist...and define a directed graph..
 el<-cbind(nodes,nnodes)
 el1<-cbind(nodes,nnodes,val)
 res<-graph.edgelist(el)

 return(res)
 
}

#--example
#a<-rnorm(100,0,1)
#b<-matrix(a,20,5)
#rownames(b)<-letters[1:nrow(b)]
#--make a distance matrix...
#bd<-as.matrix(dist(b))
#--and compute 1-nn graph
#b1.1nn<-define.1nn.graph(bd)
#--and plot...
#plot(define.1nn.graph(bd),layout=layout.fruchterman.reingold)

#--Rohan finished editing here on 9 December 2012--



######### Metabolite Identification ################################

aracyc<-read.table("D:/Shiv_data/Shiv_NUS/PhD_Project/Integrating_Omics/Analysis_R/MetabolomicsSoftwareComparison/Aracyc/aracyc_20120827_uniqueMass.txt",sep="\t",header=FALSE)
str(aracyc)

# mzlist<-read.table("D:/Shiv_data/Shiv_NUS/PhD_Project/Integrating_Omics/Analysis_R/MetabolomicsSoftwareComparison/Data/TT8_RawData_Metabolomics/Agilent/MZXML_MS1/RemovedNoisyBiorep/mzlist_fc2p05.txt",header=FALSE)
# 
# mzlist<-mzlist[order(mzlist$V1),]
# ### Finding the closest value
# 
# #Modified by Shiv --from The R Book Second edition pg 47
# 
# aracyc_cpds<-aracyc[,1]
# mz_match <- list()
# 
# ## Function 
# 
# closest<-function(xv,sv)
# {
#   #xv[which(abs(xv-sv)==min(abs(xv-sv)))]
#   ifelse(min(abs(xv-sv))<0.5,xv[which(abs(xv-sv)==min(abs(xv-sv)))],0)
# }
# 
# for(i in 1:length(mzlist))
# {
#   cpd<-closest(aracyc_cpds,mzlist[i])
#   #print(mzlist[i,1])
#   print(paste(mzlist[i],cpd,sep=" "))
#   if(cpd >0)
#   {
#     mz_match[[i]]<-cbind(aracyc[which(aracyc$V1==cpd),],mzlist[i])#aracyc[cpd,]
#   }
# }
# 
# mz_id <- do.call("rbind",mz_match) 
# colnames(mz_id)<-c('AracycMZ_Mass','Metabolite_name','MZ')
# mz_id<-mz_id[,c('MZ','AracycMZ_Mass','Metabolite_name')]
# 
# write.table(mz_id,"mz_id.txt",quote=FALSE,sep="\t",row.name=FALSE)




