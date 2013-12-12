#--generate a fake mass-feature by sample matrix with 2000 mass features and 10 columns...very fake example
#m<-matrix(runif(20000),nrow=2000,ncol=10)

#--Shiv: the columns would correspond to replicates in the analsyis
#--to compute the principal components, do the following...

#--centre 'm' so we can intepret PCA results as related to 'variance'
#m.pca.1<-princomp(m,cor=F,scores=T)

#> str(m.pca.1)
#List of 7
# $ sdev    : Named num [1:10] 0.301 0.299 0.298 0.292 0.291 ...
#  ..- attr(*, "names")= chr [1:10] "Comp.1" "Comp.2" "Comp.3" "Comp.4" ...
# $ loadings: loadings [1:10, 1:10] -0.0328 0.2486 -0.3168 0.2463 0.1944 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : NULL
#  .. ..$ : chr [1:10] "Comp.1" "Comp.2" "Comp.3" "Comp.4" ...
# $ center  : num [1:10] 0.5 0.507 0.492 0.506 0.499 ...
# $ scale   : num [1:10] 1 1 1 1 1 1 1 1 1 1
# $ n.obs   : int 2000
# $ scores  : num [1:2000, 1:10] -0.103 0.175 -0.154 0.343 0.427 ...
#  ..- attr(*, "dimnames")=List of 2
#  .. ..$ : NULL
#  .. ..$ : chr [1:10] "Comp.1" "Comp.2" "Comp.3" "Comp.4" ...
# $ call    : language princomp(x = m, cor = F, scores = T)
# - attr(*, "class")= chr "princomp"

#--do use the SVD to "do" PCA, we'll need the following:
#--first, center 'm' (but don't scale)
#ms<-scale(m,T,F)

#--compute the SVD of 'ms'
#ms.svd<-svd(ms)

#> str(ms.svd)
#List of 3
# $ d: num [1:10] 13.5 13.4 13.3 13.1 13 ...
# $ u: num [1:2000, 1:10] 0.00762 -0.01303 0.01143 -0.02547 -0.03171 ...
# $ v: num [1:10, 1:10] 0.0328 -0.2486 0.3168 -0.2463 -0.1944 ...

#--note that apart from an arbitary sign change, m.pca.1$loading is identical to ms.svd$v
#--the columns of these matrices are the "eigenvectors" of the covariance matrix of 'm'.
#plot(ms.svd$v,m.pca.1$loadings)
#abline(0,1,lty=2)

#--from ms.svd$d, we can calculate the variance captured by each component (eigenvectors) using this functions...
compute.variance.per.pc<-function(svdObj){(svdObj$d^2)/nrow(svdObj$u)}
#> compute.variance.per.pc(ms.svd)
# [1] 0.09056647 0.08947101 0.08865023 0.08546428 0.08464766 0.08380670
# [7] 0.08257358 0.08010863 0.07902905 0.07384159
 
#--and the cumulative variance across compoments using...
compute.cum.per.total.var<-function(svdObj){cumsum(svdObj$d^2)/sum(svdObj$d^2)*100}
#> compute.cum.per.total.var(ms.svd)
# [1]  10.80540  21.48010  32.05688  42.25355  52.35278  62.35168  72.20346
# [8]  81.76115  91.19003 100.00000

#--note that princomp outputs standard deviation but the we calculate variance using the svd version...so if take this into acount, the two are identical...
#plot(compute.variance.per.pc(ms.svd),m.pca.1$sdev^2)
#abline(0,1,lty=2)

#--and we can compute the "scores" from the svd using the following
ms.svd.scores<-ms.svd$u%*%diag(ms.svd$d)

#--and compare against 'm.pca.1', which are equal (aside from the arbitrary sign assignments on a given PC) 
#plot(ms.svd.scores,m.pca.1$scores)

