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
###Prepare data to riun the python script using tha aa
data_merge$V_J_lenghCDR3aa = paste(data_merge$bestVGene, data_merge$bestJGene, nchar(as.character(data_merge$aaSeqCDR3)),sep="_")
data_clonesInference_aa<-data_merge[,c("seqID","sample","aaSeqCDR3","bestVGene","bestJGene","V_J_lenghCDR3aa")]
write.table(data_clonesInference_aa,file="Data/PAAD/data_clonesInference_aa.txt",row.names = F,sep="\t")

data_merge_aa<-read.csv("Data/PAAD/ClonesInfered_PAAD_aa.csv")
data_merge_aa$V_J_lenghCDR3aa_CloneId<-paste(data_merge_aa$V_J_lenghCDR3aa,data_merge_aa$CloneId,sep="_")
data_merge<-merge(data_merge,data_merge_aa,by="seqID")

#chain=c("IGHV","IGKV","IGLV")
############
## 1. Build the matrix with the clones by samples
data_qc_chain<-data_merge[which(data_merge$chainType=="TRAV" | 
                                  data_merge$chainType=="TRBV" |
                                data_merge$chainType=="TRDV" | 
                                  data_merge$chainType=="TRGV"),]
clone_type<-t(as.data.frame(unclass(table(data_qc_chain$V_J_lenghCDR3_CloneId,factor(data_qc_chain$sample))))) 
##Build present vs no present
clone_type_presence<-apply(clone_type,1,function(x) ifelse(x==0,0,1))
###Filter by clones that at least are share in 2 samples
clone_type_filter<-clone_type_presence[which(rowSums(clone_type_presence)>1),] #
clone_type_filter2<-clone_type_filter[,which(colSums(clone_type_filter)>0)] #2898 (IG) 180 samples # 555(TCR) 164 samples

id<-match(rownames(clone_type_filter),colnames(clone_type))
write.csv(clone_type[,id], file = "Results/common_clones_TCR.csv")

tiff(paste0("Results/heatmap_common_clones_TCR.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter2),border_color=F,color = colorRampPalette(c("white", "red"))(50))
dev.off()
