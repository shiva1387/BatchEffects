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

#####
metabolite_pathways<-read.table("metabolitePathways.txt",header=TRUE,sep="\t")

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

