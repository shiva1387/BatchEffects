num_cores<-8; # The number of CPUs available, be aware of memory issue

library(xcms)

set1<-xcmsSet(nSlaves=44,method='centWave',ppm=30,peakwidth=c(5,60), prefilter=c(0,0),snthresh=6) 
set2 <- group(set1,bw=5,mzwid=0.015,minsamp=1,minfrac=0) 
set3 <- retcor(set2,method="obiwarp",plottype="none")
set4 <- group(set3,bw=5,mzwid=0.015,minsamp=1,minfrac=0)
set5 <- fillPeaks(set4) 
peaklist<-peakTable(set5,filebase="algae_blanks")
   
