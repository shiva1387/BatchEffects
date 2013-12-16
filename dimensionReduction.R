#http://www.r-bloggers.com/introduction-to-feature-selection-for-bioinformaticians-using-r-correlation-matrix-filters-pca-backward-selection/
library(corrplot)
#corrplot: the library to compute correlation matrix.

datMy <- read.table("data.csv", header = TRUE)
#read the tab file using the read table function.


datMy.scale<- scale(datMy[2:ncol(datMy)],center=TRUE,scale=TRUE)

d<-t(ms_data_wo_outliers)
datMy.scale<- scale(d[2:ncol(d)],center=TRUE,scale=TRUE);
#scale all the features (from feature 2 bacause feature 1 is the predictor output)

corMatMy <- cor(datMy.scale)
#compute the correlation matrix

corrplot(corMatMy, order = "hclust")
#visualize the matrix, clustering features by correlation index