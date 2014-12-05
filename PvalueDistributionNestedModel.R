d4.nested.pval.runday<-as.numeric(lm_day4_nested_results_pval$RunDay)
d4.nested.pval.strain<-as.numeric(lm_day4_nested_results_pval$Strain)

d12.nested.pval.runday<-as.numeric(lm_day12_nested_results_pval$RunDay)
d12.nested.pval.strain<-as.numeric(lm_day12_nested_results_pval$Strain)

#### functions
calculatePOverlaps<-function(dataset,threshold)
{
  common_05<-intersect(rownames(dataset)[which(dataset$RunDay<threshold)],
                       rownames(dataset)[which(dataset$Strain<threshold)])
  ronly_05<-setdiff(rownames(dataset)[which(dataset$RunDay<threshold)],common_05)
  sonly_05<-setdiff(rownames(dataset)[which(dataset$Strain<threshold)],common_05)
  return(cbind(length(common_05),length(ronly_05),length(sonly_05),threshold))
}

binseq<-seq(0,1,0.05)
get.hist.data<-function(x,b){hist(x,b,plot=F)$counts}

########## day4
d4.nested.pval.runday.hd<-get.hist.data(d4.nested.pval.runday,binseq)
d4.nested.pval.strain<-get.hist.data(d4.nested.pval.strain,binseq)
hd.xaxis<-(seq(0,1,0.05)+0.025)[-21]

range12<-max(c(d4.nested.pval.runday.hd,d4.nested.pval.strain))

par(mfrow=c(1,2))
plot(hd.xaxis,d4.nested.pval.strain,type="s",ylim=c(100,range12),log="y",las=1,ylab="Number of mass features",xlab="P-value",main="Strain effect")
lines(hd.xaxis,rep((0.05*13443),20),lty=2)
plot(hd.xaxis,d4.nested.pval.runday.hd,type="s",ylim=c(100,range12),log="y",las=1,ylab="Number of mass features",xlab="P-value",main="Runday effect")
lines(hd.xaxis,rep((0.05*13443),20),lty=2)

#pval threshold
features_d4_01<-calculatePOverlaps(lm_day4_nested_results_pval, 0.01)
features_d4_05<-calculatePOverlaps(lm_day4_nested_results_pval, 0.05)
features_d4_1<-calculatePOverlaps(lm_day4_nested_results_pval, 0.1)

features_d4_pvalThreshold<-rbind(features_d4_01,features_d4_05,features_d4_1)
colnames(features_d4_pvalThreshold)<-c("Common","Runday only","Strain only","threshold")
rownames(features_d4_pvalThreshold)<-c("pval< .01","pval< .05","pval< .1")


#ronly.mz_rt <- strsplit(rownames(ms_data_total), "\\@")
# mz<-sapply(mz_rt , function (x) if(length(x) == 2) x[1] else as.character(NA))
# rt<-sapply(mz_rt , function (x) if(length(x) == 2) x[2] else as.character(NA))

metadata_melt_d4<-melt(features_d4_pvalThreshold,id="threshold")
metadata_melt_d4<-metadata_melt_d4[metadata_melt_d4[,2]!="threshold",]
features_d4_pvalThreshold_plot<-ggplot(metadata_melt_d4,aes(x=Var2,value,fill=Var1))+ geom_point (aes(color=Var1,shape = Var2),size=4) + coord_cartesian(ylim = c(0, 12000)) +
  theme_bw() + theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
                                          strip.text.x = element_text(size=12, face="bold"),
                                          panel.grid.major.x = element_blank(), # to x remove gridlines
                                          panel.grid.major.y = element_blank(), # to y remove gridlines
                                          panel.border = element_blank(),  # remove top and right border
                                          panel.background = element_blank(),
                                          axis.line = element_line(color = 'black')) + xlab("Factor") + ylab("Number of mass features") + ggtitle("Day4- Nested linear model")
######## day12

d12.nested.pval.runday.hd<-get.hist.data(d12.nested.pval.runday,binseq)
d12.nested.pval.strain<-get.hist.data(d12.nested.pval.strain,binseq)
hd.xaxis<-(seq(0,1,0.05)+0.025)[-21]

range12<-max(c(d12.nested.pval.runday.hd,d12.nested.pval.strain))

par(mfrow=c(1,2))
plot(hd.xaxis,d12.nested.pval.strain,type="s",ylim=c(100,range12),log="y",las=1,ylab="Number of mass features",xlab="P-value",main="Strain effect")
lines(hd.xaxis,rep((0.05*10687),20),lty=2)
plot(hd.xaxis,d12.nested.pval.runday.hd,type="s",ylim=c(100,range12),log="y",las=1,ylab="Number of mass features",xlab="P-value",main="Runday effect")
lines(hd.xaxis,rep((0.05*10687),20),lty=2)

#pval threshold
features_d12_01<-calculatePOverlaps(lm_day12_nested_results_pval, 0.01)
features_d12_05<-calculatePOverlaps(lm_day12_nested_results_pval, 0.05)
features_d12_1<-calculatePOverlaps(lm_day12_nested_results_pval, 0.1)

features_d12_pvalThreshold<-rbind(features_d12_01,features_d12_05,features_d12_1)
colnames(features_d12_pvalThreshold)<-c("Common","Runday only","Strain only","threshold")
rownames(features_d12_pvalThreshold)<-c("pval< .01","pval< .05","pval< .1")


#ronly.mz_rt <- strsplit(rownames(ms_data_total), "\\@")
# mz<-sapply(mz_rt , function (x) if(length(x) == 2) x[1] else as.character(NA))
# rt<-sapply(mz_rt , function (x) if(length(x) == 2) x[2] else as.character(NA))

metadata_melt_d12<-melt(features_d12_pvalThreshold,id="threshold")
metadata_melt_d12<-metadata_melt_d12[metadata_melt_d12[,2]!="threshold",]
features_d12_pvalThreshold_plot<-ggplot(metadata_melt_d12,aes(x=Var2,value,fill=Var1))+ geom_point (aes(color=Var1,shape = Var2),size=4) + coord_cartesian(ylim = c(0, 12000)) +
  theme_bw() + theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),
                     strip.text.x = element_text(size=12, face="bold"),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black')) + xlab("Factor") + ylab("Number of mass features") + ggtitle("Day12- Nested linear model")


#
pdf("pvalueDistributionThreshold_nestedmodel.pdf",height=8,width=6)
multiplot(features_d4_pvalThreshold_plot, features_d12_pvalThreshold_plot) #)#,   pc1_corrected, pc2_corrected, pc3_corrected cols=2)
dev.off()

pdf("pvalueDistribution_nestedmodel.pdf",height=8,width=8)
par(mfrow=c(2,2))
plot(hd.xaxis,d4.nested.pval.strain,type="s",ylim=c(100,range12),log="y",las=1,ylab="Number of mass features",xlab="P-value",main="Day4 Strain effect")
lines(hd.xaxis,rep((0.05*13443),20),lty=2)
plot(hd.xaxis,d12.nested.pval.strain,type="s",ylim=c(100,range12),log="y",las=1,ylab="Number of mass features",xlab="P-value",main="Day12 Strain effect")
lines(hd.xaxis,rep((0.05*10687),20),lty=2)
plot(hd.xaxis,d4.nested.pval.runday.hd,type="s",ylim=c(100,range12),log="y",las=1,ylab="Number of mass features",xlab="P-value",main="Day4 Runday effect")
lines(hd.xaxis,rep((0.05*13443),20),lty=2)
plot(hd.xaxis,d12.nested.pval.runday.hd,type="s",ylim=c(100,range12),log="y",las=1,ylab="Number of mass features",xlab="P-value",main="Day12 Runday effect")
lines(hd.xaxis,rep((0.05*10687),20),lty=2)
dev.off()

