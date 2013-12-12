#http://talesofr.wordpress.com/2013/05/05/ridiculously-photogenic-factors-heatmap-with-p-values/


corfac <- data.frame(gender,var2,var3,var4,var5,var6)
summary(corfac)
class(env_parameters)
combos <- combn(ncol(env_parameters),2) # combinations without repetitions

combos <- expand.grid(rep(list(1:ncol(env_parameters)), 2 )) # combinations with repetitions
combos <- as.matrix(combos)
combos <- t(combos) # transpose matrix

mat1 <- adply(combos, 2, function(x) {
  #test <- chisq.test(env_parameters[, x[1]], env_parameters[, x[2]])
  test<-cor.test(env_parameters[, x[1]], env_parameters[, x[2]], use="all.obs", method="kendall", exact=FALSE) 
  
  out <- data.frame("Row" = colnames(env_parameters)[x[1]]
                    , "Column" = colnames(env_parameters[x[2]])
                    , "p.value" = round(test$p.value, 3)
                    ,  "statistic"= test$statistic
                    ,  "tau.value" = round(test$estimate,3)
  )
  return(out)
})

mat1.pval <- adply(combos, 2, function(x) {
  #test <- chisq.test(env_parameters[, x[1]], env_parameters[, x[2]])
  test<-cor.test(env_parameters[, x[1]], env_parameters[, x[2]], use="all.obs", method="kendall", exact=FALSE) 
  
  out <- data.frame("Row" = colnames(env_parameters)[x[1]]
                    , "Column" = colnames(env_parameters[x[2]])
                    , "tau.value" = round(test$estimate,3)
                    ,  "statistic"= test$statistic
                    ,  "p.value" = round(test$p.value, 3)
  )
  return(out)
})

mat2df <- cast(mat1, Row~Column) # Eureka!
mat2 <- as.matrix(mat2df)
head(mat2)

mat2df.pval <- cast(mat1.pval, Row~Column) # Eureka!
mat2.pval <- as.matrix(mat2df.pval)
head(mat2.pval)

mat1$stars <- cut(mat1$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
mat1$tauval.scal<-cut(mat1$tau.value,breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),include.lowest=TRUE, label=c("(-0.75,-1)","(-0.5,-0.75)","(-0.25,-0.5)","(0,-0.25)","(0,0.25)","(0.25,0.5)","(0.5,0.75)","(0.75,1)")) # the label for the legend

pdf('test.pdf')
ggplot(mat1, aes(Row, Column)) +
   geom_tile(aes(fill=mat1$tauval.scal),colour="white") +
  scale_fill_brewer(palette = "RdYlGn",name="Correlation") + geom_text(aes(label=mat1$stars), color="black", size=5.5) + 
  labs(y=NULL, x=NULL, fill="p.value") + 
  theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0))
dev.off()

p <- ggplot(aes(x=country, y=variable, fill=value), data=plot.data)
fig2 <- p + geom_tile() + scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") + 
  #   geom_text(aes(label=stars, color=value), size=8) + scale_colour_gradient(low="grey30", high="white", guide="none") +
  geom_text(aes(label=stars), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="z-value") + geom_vline(xintercept=2.5, size=1.5, color="grey50") + 
  theme_bw() + theme(axis.text.x=element_text(angle = -45, hjust = 0))
