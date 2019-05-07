rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. Common clones analysis 
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Preparing Immune and Clinical for PAAD 
###         
###
### Author: Silvia Pineda
### Date: April, 2019
############################################################################################
library(ggplot2)
library(Rtsne)
library(glmnet)
library(RColorBrewer)
library(pheatmap)

setwd("~/TCGA-Immune/")

load("Data/PAAD/PAAD_FullData.Rdata")

###################################################
## Common clones across samples
#################################################
chain="TRGV"
############
## 1. Build the matrix with the clones by samples
data_qc_chain<-data_merge[which(data_merge$chainType==chain),]
clone_type<-t(as.data.frame(unclass(table(data_qc_chain$V_J_lenghCDR3_CloneId,factor(data_qc_chain$sample))))) 
##Build present vs no present
clone_type_presence<-apply(clone_type,1,function(x) ifelse(x==0,0,1))
###Filter by clones that at least are share in 2 samples
clone_type_filter<-clone_type_presence[which(rowSums(clone_type_presence)>1),]#32
id<-match(rownames(clone_type_filter),colnames(clone_type))
write.csv(clone_type[,id], file = paste0("Results/common_clones_",chain,".csv"))

pheatmap(t(clone_type_filter))

