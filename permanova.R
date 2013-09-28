########################
# load mzXML data in R-Arabidopsis data#
########################
# Author(s): Shiv
# Version: 13052013 
# Input: ".tsv" file from TT8 orbitrap data 
# Software: XCMS
# File Location :D:\Shiv_data\Shiv_NUS\PhD_Project\Integrating_Omics\Analysis_R\MetabolomicsSoftwareComparison\Data\TT8_RawData_Metabolomics\Orbitrap\
# Output: 
# Modified By :Shivshankar Umashankar 

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


#############
# User      #
# Specific  #
# Variables #
#############

directory<-"F:/Vinay's Algae Data/Aug2013/Algae/test"; 
#Contains path for .tsv file

#directory<-"D:/Shiv_data/Shiv_NUS/PhD_Project/Integrating_Omics/Analysis_R/MetabolomicsSoftwareComparison/Data/TT8_RawData_Metabolomics/Agilent/MZXML_MS1/All/";
#Contains path for raw mzxml files

# The "PATH", Remember to use "/" and not "/" in the path
num_cores<-8; # The number of CPUs available, be aware of memory issue

setwd(directory)

#################
# Get filenames #
#################


mzxmlfiles<-list.files(directory,pattern=".mzXML",full.names=TRUE) # Get filenames with a given pattern and include subdirectories 

####### reading in the mass spec data 

mzfilename<-"Algae_day4_NONAVG.tsv"
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

# Assigning sample groups

a<-names(ms_data_total)
#strsplit(a[3],"_")
#gsub("^([^_]*_[^_]*)_.*$", "\\1", a[3])
#b<-paste(strsplit(a,"_")[[1]][1:2],collapse = "_")

SampleGroups<-sapply(a, function(x) paste(strsplit(x,"_")[[1]][1:2],collapse = "_"))
SampleGroups<-as.vector(SampleGroups)

SampleGroups<-gsub('b1','',SampleGroups)
SampleGroups<-gsub('b2','',SampleGroups)
SampleGroups<-gsub('b2','',SampleGroups)

#write.table(SampleGroups0,"SampleGroups.txt",row.name=FALSE,quote=F)
#SampleGroups1<-read.table("SampleGroups.txt",header=F)

#no_of_info_cols<-1
#data_start_index<-no_of_info_cols+1;
#data_end_index<-dim(ms_data_total)[2];

ms_data<-ms_data_total[,2:206]
names(ms_data)
dim(ms_data)
ms_data<-log(ms_data)

SampleGroup<-SampleGroups[2:206]

#name_list <- strsplit(names(ms_data), "-")
#GrowthStage<-sapply(name_list , function (x) if(length(x) == 2) x[1] else as.character(NA))
#StrainId<-sapply(name_list , function (x) if(length(x) == 2) x[2] else as.character(NA))
#StrainName<-as.factor(names(ms_data))
#StrainName0<-names(ms_data)

################PERMANOVA

#ms_data_total$DayStrain<-paste0(ms_data_total$Day,ms_data_total$Strain)
a<-adonis(t(log(ms_data+1)) ~ StrainName, data=ms_data_total,method = "bray", perm=999)
b1<-metaMDS(log(ms_data+1),k=2)

plot(b)
ordiplot(b,type="n")
orditorp(b,display="species",col="red",air=0.01)
orditorp(b,display="sites",cex=0.01,air=0.01)

distance<-vegdist(log(t(ms_data)), method="bray")
model1<-betadisper(distance, StrainName)
a<-permutest(model1, pairwise = TRUE,permutations=999)
b<-TukeyHSD(model1)

#ttest

pairwise_ttest_strains<-pairwise.t.test(t(ms_data_test),StrainName_test,p.adj="BH", paired=F, var.equal=T)

write.table(pairwise_ttest_strains$p.value,"paiwise_day4.txt",sep='\t',col.names=NA,quote=FALSE)
day4_strains<-read.table("paiwise_day4.txt",sep='\t',header=TRUE,row.names=1)

melt(pairwise_wilcoxon_strains[[3]])
d<-log(ms_data_test+1)
apply(d, 2, function(d) {t.test(x = d[,1], y = d[,2])$p.value})

##########################
# Histogram of mz and rt #
##########################

mz_list<-ms_data_total$mz
rt_list<-ms_data_total$rt

png("algae_417reps_retcor_MPP_sorted.png")
plot(ms_data_total$"Retention Time",ms_data_total$Mass,pch=1,main="algae_417reps_retcor_MPP_sorted",ylab="m/z",xlab="rt")
dev.off()
############## various plotting functions

##Plotting TIC's

#combined

getTICs(files=mzxmlfiles,pdfname=paste("xcms_agilent_TIC_algae_15th.pdf", sep = ""),rt="raw")

#ws

ws_mzxmlfiles<-mzxmlfiles[grep('WS',mzxmlfiles)] 
getTICs(files=ws_mzxmlfiles,pdfname=paste("xcms_agilent_TIC_WS_140513.pdf", sep = ""),rt="raw")

#tt8
tt8_mzxmlfiles<-mzxmlfiles[-grep('WS',mzxmlfiles)] 
getTICs(files=tt8_mzxmlfiles,pdfname=paste("xcms_agilent_TIC_TT8_140513.pdf", sep = ""),rt="raw")

### Measure of variation

png('boxplot.png',width=8000)

#boxplot(as.matrix(log(ms_data[,32:33])),range=0,xlab="Samples",ylab="Log2 feature counts",las=1,col=rainbow(9))
boxplot(t(log1p(ms_data))~SampleGroup,range=0,xlab="Samples",ylab="Log2 feature counts",las=1)
dev.off()

# plotting coefficient of variations among samples(cols)

plot(1:ncol(ms_data),apply(ms_data,2,cv))

#abline(v=seq(1,24,3)-0.5)
#axis(side=1,at=seq(2,24,3),as.character(1:8))
dev.off()

heatmap.2(as.matrix(log(day4_strains+1)),dendrogram="column",trace="none",mar=c(10,10),Rowv=FALSE,Colv=FALSE,xlab = NULL)

######### ggplot ########

ggplot(ms_data_nd_new, aes_string(x=rt,y=n)) + geom_line(aes(colour = series)) + facet_grid(series ~ .)+theme(axis.ticks = element_blank(), axis.text.y = element_blank())
lapply(list("00111","12922"), function(i) ggplot(i,aes(x=rt,value))+geom_point())
log_blanks<-log1p(blanks)
log_blanks$rt<-ms_data_total[,5]
log_blanks<- melt(log_blanks,  id = 'rt', variable_name = 'series')
ggplot(log_blanks, aes(x=rt,value)) + geom_line(aes(colour = series)) +geom_rug(sides="b")

########## Clustering ######

d <- dist(t(log(ms_data+1)), method = "euclidean") # distance matrix
fit <- hclust(d, method="average") 
#pdf('hclust.pdf',width=800)
plot (fit)
#dev.off()

clusDendro<-as.dendrogram(fit,hang=2)
labelColors<-c("red","blue")
clusMember<-rep(1,length(names(ms_data)))
clusMember[grep("D12",names(ms_data))]<-2
names(clusMember)<-names(ms_data)

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
plot(clusDendro,hang=0)

#### PCOA

bray_grps.D <- vegdist(t(log1p(ms_data)), "bray")#calculating distance matrix using bray curtis
res_grps <- pcoa(bray_grps.D)
res_grps$values
pdf(file="PCOA_samples.pdf")
biplot(res_grps)
dev.off()

pdf(file="PCOA_ordplot.pdf",width=12, heigh=12, paper="a4r")
pc<-cmdscale(bray_grps.D, k=10, eig=TRUE, add=TRUE, x.ret =TRUE)   
#Create ordination plot     
fig<-ordiplot(scores(pc)[,c(1,2)], type="t", main="PCoA Samples",cex=0.5)
dev.off()

x11()
ordiplot (scores(pc)[,c(1,2)], display = 'sp', type = 'n',main="PCoA Samples", cex=0.5)
points(scores(pc)[,c(1,2)], col = clusMember , pch = clusMember )
legend("bottomleft", legend = c("Day4","Day12"), pch = 1:2,col = c("blue","red"))


#### PCA

fit <- princomp(log1p(ms_data), cor=TRUE)
summary(fit) # print variance accounted for 
loadings(fit) # pc loadings 
plot(fit,type="lines") # scree plot 
fit$scores # the principal components
biplot(fit)

fit <- prcomp(t(log1p(ms_data)))
mycolors <- c("red", "green", "blue", "magenta", "black")
plot(fit$x, type="n"); text(fit$x, rownames(fit$x), cex=0.5)

library(psych)
fit <- principal(log1p(samples_avg), nfactors=5, rotate="varimax")
fit


####nn-graphs

ms_graph_data_strains<-t(log1p(ms_data))

rownames(ms_graph_data_strains)<-rownames(ms_graph_data_strains)
ms.d_str<-as.matrix(dist(ms_graph_data_strains,"euclidean"))
ms.1nn_str<-define.1nn.graph(Dmat)

pdf(file="NN_graph_strains.pdf",paper="a4r")
plot(define.1nn.graph(ms.d_str),layout=layout.fruchterman.reingold,vertex.label=V(ms.1nn_str)$name,vertex.label.cex=0.8,main="Strains")
dev.off()

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




