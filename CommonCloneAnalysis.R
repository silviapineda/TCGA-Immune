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
library(factoextra)
##Colors: brewer.pal(4, "Accent")
##"#BEAED4" (IGH) "#7FC97F" (IGK) "#FDC086" (IGL)
setwd("~/TCGA-Immune/")

load("Data/PAAD/PAAD_FullData.Rdata")

###################################################
## Common clones across samples
#################################################
###Prepare data to riun the python script using tha aa
#data_merge$V_J_lenghCDR3aa = paste(data_merge$bestVGene, data_merge$bestJGene, nchar(as.character(data_merge$aaSeqCDR3)),sep="_")
#data_clonesInference_aa<-data_merge[,c("seqID","sample","aaSeqCDR3","bestVGene","bestJGene","V_J_lenghCDR3aa")]
#write.table(data_clonesInference_aa,file="Data/PAAD/data_clonesInference_aa.txt",row.names = F,sep="\t")

#data_merge_aa<-read.csv("Data/PAAD/ClonesInfered_PAAD_aa.csv")
#data_merge_aa$V_J_lenghCDR3aa_CloneId<-paste(data_merge_aa$V_J_lenghCDR3aa,data_merge_aa$CloneId,sep="_")
#data_merge<-merge(data_merge,data_merge_aa,by="seqID")

#chain=c("IGHV","IGKV","IGLV")
############
## 1. Build the matrix with the clones by samples
#data_qc_chain<-data_merge[which(data_merge$chainType=="IGHV" | 
#                                  data_merge$chainType=="IGKV" |
#                                data_merge$chainType=="IGLV"),]
data_qc_chain_IGKV<-data_merge[which(data_merge$chainType=="IGKV"),]
clone_type_IGKV<-t(as.data.frame(unclass(table(data_qc_chain$V_J_lenghCDR3_CloneId,factor(data_qc_chain$sample))))) 
data_qc_chain_IGLV<-data_merge[which(data_merge$chainType=="IGLV"),]
clone_type_IGLV<-t(as.data.frame(unclass(table(data_qc_chain$V_J_lenghCDR3_CloneId,factor(data_qc_chain$sample))))) 
data_qc_chain_IGHV<-data_merge[which(data_merge$chainType=="IGHV"),]
clone_type_IGHV<-t(as.data.frame(unclass(table(data_qc_chain$V_J_lenghCDR3_CloneId,factor(data_qc_chain$sample))))) 

##Build present vs no present
clone_type_presence<-apply(clone_type_IGKV,1,function(x) ifelse(x==0,0,1))
clone_type_presence<-apply(clone_type_IGLV,1,function(x) ifelse(x==0,0,1))
clone_type_presence<-apply(clone_type_IGHV,1,function(x) ifelse(x==0,0,1))

###Filter by clones that at least are share in 2 samples
clone_type_filter<-clone_type_presence[which(rowSums(clone_type_presence)>1),] #
clone_type_filter_IGKV<-clone_type_filter[,which(colSums(clone_type_filter)>0)] #21436 180
clone_type_filter_IGLV<-clone_type_filter[,which(colSums(clone_type_filter)>0)] #1031 180
clone_type_filter_IGHV<-clone_type_filter[,which(colSums(clone_type_filter)>0)] #431 141

##Heatmp
tiff(paste0("Results/heatmap_common_clones_IGH.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGHV),border_color=F,show_rownames = F, show_colnames = F,
         color = c("white","#BEAED4"),breaks = c(0,0.9,1))
dev.off()
tiff(paste0("Results/heatmap_common_clones_IGK.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGKV),show_rownames = F, show_colnames = F,border_color=F,
         color = c("white","#7FC97F"),breaks = c(0,0.9,1))
dev.off()
tiff(paste0("Results/heatmap_common_clones_IGL.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGLV),border_color=F,show_rownames = F, show_colnames = F,
         color = c("white","#FDC086"),breaks = c(0,0.9,1))
dev.off()

####Study clones ####
##IGL
d <- dist(clone_type_filter_IGLV, method = "euclidean") # distance matrix
fit.k <- hclust(d) 
plot(fit.k) # display dendogram
groups <- cutree(fit.k, k=3) # 2/3 61/86 clones 
id.clones<-match(names(which(groups!=1)),colnames(clone_type_IGLV))
clones_IGL<-clone_type_IGLV[,id.clones]
pheatmap(t(clones_IGL),show_rownames = F, show_colnames = F,border_color=F,
         color = colorRampPalette(c("white", "#FDC086"))(20))
pca<-prcomp(clone_type_IGLV)
tiff("Results/Clone_type_IGL_PCA.tiff",width = 6000, height = 2500, res = 300)
fviz_pca_ind(pca)
dev.off()
which(pca$x[,2]>6000) #d7be9e8f-9c5f-4e1c-940b-7ef67fad520d
which(pca$x[,1]>5000) #5e9e81e2-0e8b-4aca-aced-6ce451fa3262
##
id.clones<-match(names(which(groups!=1)),data_qc_chain_IGLV$V_J_lenghCDR3_CloneId)
data_IGLV<-data_qc_chain_IGLV[id.clones,]
summary(data_IGLV$CDR3_length)
tiff(paste0("Results/barplot_clones_IGL_Vgenes.tiff"),width = 4500, height = 2500, res = 300)
barplot(table(factor(data_IGLV$bestVGene)),las=2)
dev.off()
tiff(paste0("Results/barplot_clones_IGL_jgenes.tiff"),width = 4500, height = 2500, res = 300)
barplot(table(factor(data_IGLV$bestJGene)),las=2)
dev.off()


##IGK
d <- dist(clone_type_filter_IGKV, method = "euclidean") # distance matrix
fit.k <- hclust(d) 
plot(fit.k) # display dendogram
groups <- cutree(fit.k, k=2) #190 clones 
id.clones<-match(names(which(groups==2)),colnames(clone_type_IGKV))
clones_IGK<-clone_type_IGKV[,id.clones]
pheatmap(t(clones_IGK),show_rownames = F, show_colnames = F,border_color=F,
         color = colorRampPalette(c("white", "#7FC97F"))(20000))
pca<-prcomp(clone_type_IGKV)
fviz_pca_ind(pca)
which(pca$x[,2]>15000) #08c4c14a-94a8-4063-bbc3-916579960078
which(pca$x[,1]>15000) #5e9e81e2-0e8b-4aca-aced-6ce451fa3262
##
id.clones<-match(names(which(groups==2)),data_qc_chain_IGKV$V_J_lenghCDR3_CloneId)
data_IGKV<-data_qc_chain_IGKV[id.clones,]
summary(data_IGKV$CDR3_length)
tiff(paste0("Results/barplot_clones_IGK_Vgenes.tiff"),width = 4500, height = 2500, res = 300)
barplot(table(factor(data_IGKV$bestVGene)),las=2)
dev.off()
tiff(paste0("Results/barplot_clones_IGK_jgenes.tiff"),width = 4500, height = 2500, res = 300)
barplot(table(factor(data_IGKV$bestJGene)),las=2)
dev.off()




id<-match(rownames(clone_type_filter2),colnames(clone_type))
id.sample<-match(colnames(clone_type_filter2),rownames(clone_type))
clone_type_abundance_filter<-clone_type[id.sample,id]
write.csv(clone_type_abundance_filter, file = "Results/common_clones_TCR.csv")


############
### 2. Normalization by transforming to relative abundance
############
clone_type_relative<-t(apply(clone_type_abundance_filter,1,function (x) x/sum(x))) #180/2898
##Only pancreas tumors
clone_type_relative_tumors<-clone_type_relative[which(PAAD.repertoire.diversity$Tumor_type=="Tumor_pancreas"),] #148

write.csv(clone_type[,id], file = "Results/common_clones_TCR.csv")



############
## Nonmetric MDS
############
library(MASS)
library(vegan)
d <- vegdist(clone_type_relative_tumors, method="bray") # euclidean distances between the rows
fit <- isoMDS(d, k=2) # k is the number of dim
#fit <- cmdscale(d,eig=TRUE, k=2)
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
tiff("Results/NMDS.tiff",width = 2000, height =2000, res=300)
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Nonmetric MDS")
dev.off()
text(x, y, labels = row.names(clone_type_relative_tumors), cex=.7)

############
## Compositioanl
############
library(Compositional)
additive_log_ratio<-alr(clone_type_abundance_filter+1) ##Adding 1 for the zeros


library("DirichletMultinomial")
fit1<-dmn(clone_type_abundance_filter,k=1)
fit2<-dmn(clone_type_abundance_filter,k=2)
fit3<-dmn(clone_type_abundance_filter,k=3)
fit4<-dmn(clone_type_abundance_filter,k=4)
fit5<-dmn(clone_type_abundance_filter,k=5)
fit<-list(fit1,fit2,fit3,fit4,fit5)
lplc <- sapply(fit, laplace)
best <- fit[[which.min(lplc)]]
heatmapdmn(clone_type_abundance_filter, fit[[1]], best, 30)
