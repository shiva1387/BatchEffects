
######

library(xcms)
library(gtools)
library(stringr)

### Groupling samples accoring to species
directory="F:/Vinay's Algae Data/Aug2013/Algae/MZXML";

mzxmlfiles<-list.files(directory,pattern=".mzXML",full.names=TRUE)

mzxmlfile_names<-mzxmlfiles

#mzxmlfile_names<-sub("MassProfileStage1/","",mzxmlfiles)

#mzxmlfile_names<-sub("DataSET","",mzxmlfile_names)

Blanks<-mzxmlfile_names[grep('blank',mzxmlfile_names)]

#write.table(Blanks,file="blank_file_names.txt",quote=F,row.names=F,col.names=F)

non_blanks<-mzxmlfile_names[-grep('blank', mzxmlfile_names)] 

#duplicates<-c('MassProfileStage1/DataSET1-vin001tp1tr1.mzXML','MassProfileStage1/DataSET2-vin002tp1tr1.mzXML')

#non_blanks<-non_blanks[!non_blanks %in% duplicates]

non_blanks<-mixedsort(mzxmlfile_names)

head(non_blanks)

getTICs(files=Blanks,pdfname="RAW_TICS/Blanks.pdf",rt="raw") 

############################################ plotting TICS ##############
i<-1

nm0<- paste(non_blanks[175:176],non_blanks[179:180],non_blanks[183:184])

nm<-unlist(strsplit(nm0, " "))

cat(nm,"\n")

getTICs(files=nm,pdfname=paste("../RAW_TICS/", "D12_001-15th.pdf", sep = ""),rt="raw")

#paste("MassProfileStage1/",Sample177,sep="")




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