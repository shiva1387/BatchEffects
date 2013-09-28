library(plyr)

mzfilename<-"Comparison.txt"
ms_data_total<-read.table(mzfilename,sep="\t",header=T,check.names=FALSE,row.names=1)
str(ms_data_total)

ms_data_no_log<-ms_data_test

ms_data_test<-log(ms_data_test+1)

combos <- combn(ncol(ms_data_test),2)

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

t.test(log(ms_data[,1]),log(ms_data[,5]),paired=F,var.equal=F,p.adjust="BH")

t.test(ms_data[,1],ms_data[,2],paired=F,var.equal=F,p.adjust="BH")

boxplot(ms_data[,1],ms_data[,2])
x11()
boxplot(log10(ms_data[,1]),log(10ms_data[,2]))


plot(1:ncol(ms_data),apply(ms_data,2,cv))

wilcox.test(x, y,