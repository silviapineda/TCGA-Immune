rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: GTEX blood analysis
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: 
###         
###
### Author: Silvia Pineda
### Date: May, 2019
############################################################################################
library(ggplot2)
library(Rtsne)
library(glmnet)
library(RColorBrewer)
library(pheatmap)
library(factoextra)
##Colors: brewer.pal(4, "Accent")
##"#BEAED4" (IGH) "#7FC97F" (IGK) "#FDC086" (IGL)
setwd("~/TCGA-Immune/")
load("Data/GTEX/MIXCR/GTEX_FullData.Rdata") 

############
## 1. Build the matrix with the clones by samples
data_IGK<-data_merge[which(data_merge$Gene=="IGK"),]
clone_type_IGK<-t(as.data.frame(unclass(table(data_IGK$V_J_lenghCDR3_CloneId,factor(data_IGK$Sample))))) 
data_IGL<-data_merge[which(data_merge$chainType=="IGL"),]
clone_type_IGL<-t(as.data.frame(unclass(table(data_IGL$V_J_lenghCDR3_CloneId,factor(data_IGL$Sample))))) 
data_IGH<-data_merge[which(data_merge$chainType=="IGH"),]
clone_type_IGH<-t(as.data.frame(unclass(table(data_IGH$V_J_lenghCDR3_CloneId,factor(data_IGH$Sample))))) 

##Build present vs no present
###Filter by clones that at least are share in 2 samples
clone_type_presence<-apply(clone_type_IGK,1,function(x) ifelse(x==0,0,1))
clone_type_filter<-clone_type_presence[which(rowSums(clone_type_presence)>1),] #
clone_type_filter_IGK<-clone_type_filter[,which(colSums(clone_type_filter)>0)] #2697  453

clone_type_presence<-apply(clone_type_IGL,1,function(x) ifelse(x==0,0,1))
clone_type_filter<-clone_type_presence[which(rowSums(clone_type_presence)>1),] #
clone_type_filter_IGL<-clone_type_filter[,which(colSums(clone_type_filter)>0)] #2497  453

clone_type_presence<-apply(clone_type_IGH,1,function(x) ifelse(x==0,0,1))
clone_type_filter<-clone_type_presence[which(rowSums(clone_type_presence)>1),] #
clone_type_filter_IGH<-clone_type_filter[,which(colSums(clone_type_filter)>0)] #3378  423


##Heatmap
tiff(paste0("Results/heatmap_common_clones_IGH_GTEx.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGH),border_color=F,show_rownames = F, show_colnames = F,
         color = c("white","#BEAED4"),breaks = c(0,0.9,1))
dev.off()
tiff(paste0("Results/heatmap_common_clones_IGK_GTEx.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGK),show_rownames = F, show_colnames = F,border_color=F,
         color = c("white","#7FC97F"),breaks = c(0,0.9,1))
dev.off()
tiff(paste0("Results/heatmap_common_clones_IGL_GTEx.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGL),border_color=F,show_rownames = F, show_colnames = F,
         color = c("white","#FDC086"),breaks = c(0,0.9,1))
dev.off()
