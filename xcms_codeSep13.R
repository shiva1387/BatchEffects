#To run xcms first export the library path as
#export LD_LIBRARY_PATH=/data/metagenome_data/netcdf/lib/
#shiv 20-Jan-2015
#edited date on 20-Jan-2015
#Then start R and load xcms library as usual

library(xcms)
library(CAMERA)

## Code runon 200115 with the following package versions

##R version 3.0.2 (2013-09-25)
##Platform: x86_64-unknown-linux-gnu (64-bit)

##locale:
## [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
## [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
## [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
## [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
## [9] LC_ADDRESS=C               LC_TELEPHONE=C            
##[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

##attached base packages:
##[1] parallel  stats     graphics  grDevices utils     datasets  methods  
##[8] base     

##other attached packages:
##[1] CAMERA_1.18.0      igraph_0.7.1       xcms_1.38.0        Biobase_2.22.0    
##[5] BiocGenerics_0.8.0 mzR_1.8.1          Rcpp_0.11.3       

##loaded via a namespace (and not attached):
## [1] acepack_1.3-3.3     cluster_1.15.3      codetools_0.2-9    
## [4] foreign_0.8-61      Formula_1.1-2       graph_1.40.1       
## [7] grid_3.0.2          Hmisc_3.14-5        lattice_0.20-29    
## [10] latticeExtra_0.6-26 nnet_7.3-8          RBGL_1.38.0        
## [13] RColorBrewer_1.0-5  rpart_4.1-8         splines_3.0.2      
## [16] stats4_3.0.2        survival_2.37-7     tools_3.0.2    

setwd('/data/metagenome_data/Desktop/Shiv_algae/Sep13_8strains/data/MZXML')

set1<-xcmsSet(nSlaves=44,method='centWave',ppm=30,peakwidth=c(5,60), prefilter=c(0,0),snthresh=6)
save(set1,file="set1_algae_8strainsSep13_200115_xcms1_38.rda")

set2 <- group(set1,bw=5,mzwid=0.015,minsamp=1,minfrac=0) 
set3 <- retcor(set2,method="obiwarp",plottype="none")
set4 <- group(set3,bw=5,mzwid=0.015,minsamp=1,minfrac=0)
save(set4,file="set4_algae_8strainsSep13_200115_xcms1_38.rda")

set5 <- fillPeaks(set4) 
save(set5,file="set5_algae_8strainsSep13_200115_xcms1_38.rda")

peaklist<-peakTable(set5,filebase="algae_8strainsSep13_200115_xcms1_38")
   

xsa<-xsAnnotate(set5)
xsaF<-groupFWHM(xsa, perfwhm=0.6)
xsaC<-groupCorr(xsaF)
xsaFI<-findIsotopes(xsaC)
xsaFA<-findAdducts(xsaFI, polarity='positive')
peaklist<-getPeaklist(xsaFA)
write.csv(peaklist, file='algae_8strainsSep13_200115_xcms1_38.csv')
