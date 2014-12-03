#########################################################
# Metabolite pathway analysis in R-Malaysian algae data #
#########################################################
# Author(s): Shiv
# Version: 19052014
# Input: ".txt" file 
# Modified By :Shivshankar Umashankar 
# Functions written here are used for analyzing pathway data by plaotting bar plots

#############
# Clear all #
#############
rm(list = ls()) # Clears the workspace
graphics.off() # Close all windows

##### Loading required packages 

library(ggplot2)
library(reshape2)
library(CCA)

#####
metabolite_pathways<-read.table("metabolitePathways.txt",header=TRUE,sep="\t")

metabolite_pathways_d12<-metabolite_pathways[,c(1,3,4,5)]

metabolite_pathways_d4<-metabolite_pathways[,c(1,3,6,7)]
metabolite_pathways_d4$Day4.FullData<-round((metabolite_pathways_d4$Day4.FullData/metabolite_pathways_d4$TotalMetabolites)*100,2)
metabolite_pathways_d4$Day4.Differential<-round((metabolite_pathways_d4$Day4.Differential/metabolite_pathways_d4$TotalMetabolites)*100,2)
metabolite_pathways_d4$TotalMetabolites<-round(100-(metabolite_pathways_d4$Day4.FullData+metabolite_pathways_d4$Day4.Differential))

test1<-melt(metabolite_pathways_d4)

p1<-ggplot(test1, aes(x=KEGG.ID, y=value, fill=variable)) + geom_bar(stat="identity") + scale_fill_manual(values=c("grey","#999999", "#E69F00"))+
  theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1,size=10),axis.text.y=element_text(size=18),axis.title.y=element_blank(),
                     panel.grid.major.x = element_blank(), # to x remove gridlines
                     panel.grid.major.y = element_blank(), # to y remove gridlines
                     panel.border = element_blank(),  # remove top and right border
                     panel.background = element_blank(),
                     axis.line = element_line(color = 'black'))
ggsave("day4_pathway.pdf", p1)

############## Analyzing enrichments

#Calculating hypergeometric score for each cazy class
#ref http://mengnote.blogspot.sg/2012/12/calculate-correct-hypergeometric-p.html
#ref http://stats.stackexchange.com/questions/72553/which-statistical-test-should-be-used-to-test-for-enrichment-of-gene-lists

## The enrichment analysis is limited to the total number of metabolites that were detected in these pathways

## For day4
head(metabolite_pathways_d12)
# KEGG.ID TotalMetabolites Day4.FullData Day4.Differential
# 1 cvr00906              115            24                13
# 2 cvr02010              122            15                 7
# 3 cvr00130               80            14                 5
# 4 cvr00260               51            13                 5
# 5 cvr00592               40            11                 6
# 6 cvr00960               68            12                 5


#total_metab<-sum(metabolite_pathways_d12[,2])
total_detected_metab<-sum(metabolite_pathways_d12[,3])
total_diff_metab<-sum(metabolite_pathways_d12[,4])

#enrichment_value<-rep(99,times=nrow(metabolite_pathways_d12))
enrichment_value_fisher<-rep(99,times=nrow(metabolite_pathways_d12))

for(i in 1:nrow(metabolite_pathways_d12))
{
  hitinSample<-metabolite_pathways_d12[i,4]
  hitinPop<-metabolite_pathways_d12[i,3]
  failinPop<-total_detected_metab-hitinPop
  sampSize<-total_diff_metab
  
  a<-matrix(c(hitinSample,hitinPop-hitinSample,
              sampSize-hitinSample,failinPop-sampSize+hitinSample),nrow=2,ncol=2)
  #enrichment_value[i]<-phyper(hitinSample-1,hitinPop,failinPop,sampSize,lower.tail=FALSE)
  
  fisher<-fisher.test(a,alternative="greater")
  enrichment_value_fisher[i]<-fisher$p.value
  
}

hitinSample<-metabolite_pathways_d12[1,4]
hitinPop<-metabolite_pathways_d12[1,3]
failinPop<-total_detected_metab-hitinPop
sampSize<-total_diff_metab

a<-matrix(c(hitinSample,hitinPop-hitinSample,
            sampSize-hitinSample,failinPop-sampSize+hitinSample),nrow=2,ncol=2)
rownames(a)<-c("sig","nonSig")
colnames(a)<-c("carotenoid","otherPathways")
#8,5,105,142
#13,11,100,136
#16,3,163,66
# fisher.test(matrix(c(3,3,176,66),nrow=2,ncol=2),alternative="greater")
# fisher.test(matrix(c(12,1,167,143),nrow=2,ncol=2),alternative="greater")
# fisher.test(a,alternative="greater")

enrichment_value_fisher_d12<-enrichment_value_fisher
enrichment_value_fisher_d4<-enrichment_value_fisher

metabolite_pathways_enrich<-cbind(metabolite_pathways,round(enrichment_value_fisher_d12,2),round(enrichment_value_fisher_d4,2))

#metabolite_pathways_d12_05<-metabolite_pathways_d12_enrich[metabolite_pathways_d12_enrich[,5]<0.05,]

write.table(metabolite_pathways_enrich,"metabolite_pathways_enrichFisher.txt",quote=FALSE,sep='\t',col.names=NA)



#### Analyzing correlations

mat_cor0<-matcor(as.matrix(mappedMetabolites_mean_d12),as.matrix(biochemData_d12))
mat_cor1<-as.matrix(mat_cor0$XYcor)
mat_cor_pos<-which( mat_cor1 > 0.85,arr.ind=TRUE) #getting the positive indices
mat_cor_pos_names<-cbind(rownames(mat_cor1)[mat_cor_pos[,1]],colnames(mat_cor1)[mat_cor_pos[,2]]) #getting the names for the indices
mat_cor_pos_names<-mat_cor_pos_names[!mat_cor_pos_names[,1]==mat_cor_pos_names[,2],] #removing circular references
# unique number of eentities involved in postie correlation freater than 0.85
#length(unique(c(mat_cor_pos_names[,1],mat_cor_pos_names[,2])))
#>DAY4--520
#>DAY12-497
write.table(mat_cor_pos_names,"matcor_pos_day12.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

mat_cor_neg<-which(mat_cor1 < -0.85,arr.ind=TRUE)
mat_cor_neg_names<-cbind(rownames(mat_cor1)[mat_cor_neg[,1]],colnames(mat_cor1)[mat_cor_neg[,2]])
mat_cor_neg_names<-mat_cor_neg_names[!mat_cor_neg_names[,1]==mat_cor_neg_names[,2],] #removing circular references
# unique number of eentities involved in postie correlation freater than 0.85
#length(unique(c(mat_cor_neg_names[,1],mat_cor_neg_names[,2])))
#>DAY4--31
#>DAY12--43
write.table(mat_cor_neg_names,"matcor_neg_day12.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

matcor_pos_neg<-rbind(mat_cor_pos_names,mat_cor_neg_names)
matcor_pos_neg<-matcor_pos_neg[complete.cases(matcor_pos_neg),]
write.table(matcor_pos_neg,"matcor_pos_neg_day12.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

#length(unique(c(matcor_pos_neg[,1],matcor_pos_neg[,2])))
#>DAY4-524
#>DAY12-515

Q<-qgraph(f)


# Correlations:

data(big5)
data(big5groups)

# Correlations:
Q <- qgraph(cor(big5), minimum = 0.25, cut = 0.4, vsize = 1.5, groups = big5groups, 
            legend = TRUE, borders = FALSE)
title("Big 5 correlations", line = 2.5)

# Venn diagram -Total metabolites
library(vennDiagram)
venn_d<-read.table("TotalMetabolites.txt",header=TRUE,sep="\t",row.name=1)
venn_d1<-as.data.frame(venn_d[,c(1,2,5,6)])
a<-vennCounts(venn_d1, include="both")
png("day4_12.png")
vennDiagram(a, include="both", mar=rep(1,4), cex=1.5,counts.col=heat.colors(9),circle.col="green")
dev.off()

#pheatmap
library('pheatmap')
pdf("totalMetabolites.pdf")
pheatmap(as.matrix(venn_d1),scale="none",col=(c("white","grey")),show_rownames=F,,fontsize_row=1)
dev.off()

## Getting strains which show the maximum deviation from exponential to stationary growth phase
a<-biochemData_d12-biochemData_d4[,1:6] #as biochemData_d4 has an extra column biochemData_d4[,7] for growthRate
a1<-abs(a)#Getting the absolute modulus of change
apply(a1,2,max)
# biomass biomassProduc   lipidProduc totalLipidCon  totalProtein  totalCarbCon 
# 511.9585714     0.2641983   120.1044073    34.6839925    34.4285945    28.8695324
b<-apply(a1,2,which.max)# getting row indices which hax the max value for each column
# biomass biomassProduc   lipidProduc totalLipidCon  totalProtein  totalCarbCon 
# 11             4            16            16             1            21 
# D12_187 D12_051 D12_254 D12_254 D12_001 D12_322
rownames(a1)[b] #getting the rownames for the rowindices







