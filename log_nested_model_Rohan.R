
#--a nested example...
scores<-c(25,29,14,11,11,6,22,18,17,20,5,2)
school<-factor(c("A","A","A","A","B","B","B","B","C","C","C","C"))
teacher<-factor(c(1,1,2,2,3,3,4,4,5,5,6,6))
teacher2<-factor(c(1,1,2,2,1,1,2,2,1,1,2,2))  # This is the way the data is coded for problems in the book

res1<-lm(scores~school+school/teacher)

#--take the anova(res1)
#res1.anova<-anova(res1)

#> res1.anova
#Analysis of Variance Table
#Response: scores
#               Df Sum Sq Mean Sq F value   Pr(>F)    
#school          2  156.5   78.25  11.179 0.009473 ** 
#school:teacher  3  567.5  189.17  27.024 0.000697 ***
#Residuals       6   42.0    7.00                     
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-- so if you take
#> res1.anova$F[1:2]
#[1] 11.17857 27.02381
#> res1.anova$P[1:2]
#[1] 0.0094725376 0.0006970135

#--you get the separate F--statistics and P-values, for each of the factors
