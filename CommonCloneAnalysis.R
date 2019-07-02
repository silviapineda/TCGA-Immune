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
library(reshape2)

##Colors: brewer.pal(4, "Accent")
##"#BEAED4" (IGH) "#7FC97F" (IGK) "#FDC086" (IGL)
setwd("~/TCGA-Immune/")

load("Data/PAAD/PAAD_FullData.Rdata")

###################################################
## Common clones across samples
#################################################
###Prepare data to run the python script using the aa
data_merge$V_J_lenghCDR3aa = paste(data_merge$bestVGene, data_merge$bestJGene, nchar(as.character(data_merge$aaSeqCDR3)),sep="_")
data_merge_cdr3_Ig<-data_merge[which(data_merge$chainType=="IGH" |
                                       data_merge$chainType=="IGK" |
                                       data_merge$chainType=="IGL"),]
data_merge_cdr3_TCR<-data_merge[which(data_merge$chainType=="TRA" | 
                                       data_merge$chainType=="TRB" |
                                       data_merge$chainType=="TRD" | 
                                       data_merge$chainType=="TRG"),]

#data_clonesInference_Ig_aa<-data_merge_cdr3_Ig[,c("seqID","sample","aaSeqCDR3","bestVGene","bestJGene","V_J_lenghCDR3aa")]
#write.table(data_clonesInference_Ig_aa,file="Data/PAAD/data_clonesInference_Ig_aa.txt",row.names = F,sep="\t")
#data_clonesInference_TCR_aa<-data_merge_cdr3_TCR[,c("seqID","sample","aaSeqCDR3","bestVGene","bestJGene","V_J_lenghCDR3aa")]
#write.table(data_clonesInference_TCR_aa,file="Data/PAAD/data_clonesInference_TCR_aa.txt",row.names = F,sep="\t")

data_merge_Ig_aa<-read.csv("Data/PAAD/ClonesInfered_Ig_aa.csv")
data_merge_TCR_aa<-read.csv("Data/PAAD/ClonesInfered_TCR_aa.csv")
data_merge_aa<-rbind(data_merge_Ig_aa,data_merge_TCR_aa)
data_merge<-merge(data_merge,data_merge_aa[,c("seqID","CloneId")],by=c("seqID"))
data_merge$V_J_lenghCDR3aa_CloneId<-paste(data_merge$V_J_lenghCDR3aa,data_merge$CloneId.y,sep="_")

##Filter for those that are pancreas or normal-pseudonormal
data_merge_qc<-data_merge[which(is.na(match(data_merge$sample,rownames(PAAD.repertoire.diversity)))==F),]
#data_merge_qc<-data_merge[which(data_merge$sample %in% rownames(PAAD.repertoire.diversity)[which(PAAD.repertoire.diversity$Tumor_type_2categ=="Tumor_pancreas")]),]

#chain=c("IGHV","IGKV","IGLV")
############
## 1. Build the matrix with the clones by samples
data_qc_chain_IGK<-data_merge_qc[which(data_merge_qc$chainType=="IGK"),]
clone_type_IGK<-t(as.data.frame(unclass(table(data_qc_chain_IGK$V_J_lenghCDR3aa_CloneId,factor(data_qc_chain_IGK$sample))))) 
data_qc_chain_IGL<-data_merge_qc[which(data_merge_qc$chainType=="IGL"),]
clone_type_IGL<-t(as.data.frame(unclass(table(data_qc_chain_IGL$V_J_lenghCDR3aa_CloneId,factor(data_qc_chain_IGL$sample))))) 
data_qc_chain_IGH<-data_merge_qc[which(data_merge_qc$chainType=="IGH"),]
clone_type_IGH<-t(as.data.frame(unclass(table(data_qc_chain_IGH$V_J_lenghCDR3aa_CloneId,factor(data_qc_chain_IGH$sample))))) 
data_qc_chain_IG<-data_merge_qc[which(data_merge_qc$chainType=="IGH" |
                                        data_merge_qc$chainType=="IGK" |
                                        data_merge_qc$chainType=="IGL"),]
clone_type_IG<-t(as.data.frame(unclass(table(data_qc_chain_IG$V_J_lenghCDR3aa_CloneId,factor(data_qc_chain_IG$sample))))) 

data_qc_chain_TCR<-data_merge_qc[which(data_merge_qc$chainType=="TRA" | 
                                         data_merge_qc$chainType=="TRB" |
                                         data_merge_qc$chainType=="TRD" | 
                                         data_merge_qc$chainType=="TRG"),]
clone_type_TCR<-t(as.data.frame(unclass(table(data_qc_chain_TCR$V_J_lenghCDR3aa_CloneId,factor(data_qc_chain_TCR$sample))))) 

##Build present vs no present
clone_type_presence_IGK<-apply(clone_type_IGK,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGL<-apply(clone_type_IGL,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGH<-apply(clone_type_IGH,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IG<-apply(clone_type_IG,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TCR<-apply(clone_type_TCR,1,function(x) ifelse(x==0,0,1))

###Filter by clones that at least are share in more than 2 samples
clone_type_filter_IGK<-clone_type_presence_IGK[which(rowSums(clone_type_presence_IGK)>2),] #
clone_type_filter_IGK<-clone_type_filter_IGK[,which(colSums(clone_type_filter_IGK)>0)] #
clone_type_filter_IGL<-clone_type_presence_IGL[which(rowSums(clone_type_presence_IGL)>2),] #
clone_type_filter_IGL<-clone_type_filter_IGL[,which(colSums(clone_type_filter_IGL)>0)] #
clone_type_filter_IGH<-clone_type_presence_IGH[which(rowSums(clone_type_presence_IGH)>2),] #
clone_type_filter_IGH<-clone_type_filter_IGH[,which(colSums(clone_type_filter_IGH)>0)] #
clone_type_filter_IG<-clone_type_presence_IG[which(rowSums(clone_type_presence_IG)>2),] #
clone_type_filter_IG<-clone_type_filter_IG[,which(colSums(clone_type_filter_IG)>0)] #
clone_type_filter_TCR<-clone_type_presence_TCR[which(rowSums(clone_type_presence_TCR)>2),] #
clone_type_filter_TCR<-clone_type_filter_TCR[,which(colSums(clone_type_filter_TCR)>0)] #

##Heatmp
annotation_row = data.frame(PAAD.repertoire.diversity$Tumor_type_2categ)
ann_colors = list (Tumor_type_2categ = c("normal_pseudonormal_pancreas" = brewer.pal(3,"Accent")[1],
                                         "Tumor_pancreas"= brewer.pal(3,"Accent")[2]))
colnames(annotation_row)<-"Tumor_type_2categ"
rownames(annotation_row)<-rownames(PAAD.repertoire.diversity)

cols = brewer.pal(10,name = "Paired")

#IG
tiff(paste0("Results/heatmap_common_clones_IG.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IG),border_color=F,show_rownames = F, show_colnames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols[6]),breaks = c(0,0.9,1))
dev.off()

#IGH
tiff(paste0("Results/heatmap_common_clones_IGH.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGH),border_color=F,show_rownames = F, show_colnames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols[2]),breaks = c(0,0.9,1))
dev.off()

#IGK
tiff(paste0("Results/heatmap_common_clones_IGK.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGK),show_rownames = F, show_colnames = F,border_color=F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols[4]),breaks = c(0,0.9,1))
dev.off()

#IGL
tiff(paste0("Results/heatmap_common_clones_IGL.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGL),border_color=F,show_rownames = F, show_colnames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols[8]),breaks = c(0,0.9,1))
dev.off()

#TCR
tiff(paste0("Results/heatmap_common_clones_TCR.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_TCR),border_color=F,show_rownames = F, show_colnames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols[10]),breaks = c(0,0.9,1))
dev.off()

#####################################
### Principal Components Analysis ###
#####################################
library(factoextra)
pca <- prcomp(t(clone_type_filter_IGL), scale = TRUE)

SPP <- PAAD.repertoire.diversity$Tumor_type_2categ
levels.SPP <- factor(c("normal_pseudonormal_pancreas", "Tumor_pancreas"))
#cols<-brewer.pal(3,name = "Accent")
pc <- c(1,2)
plot(pca$x[,pc[1]], pca$x[,pc[2]], col=cols[SPP],pch=20,xlab="PCA1",ylab="PCA2")
legend("topleft", legend=levels(levels.SPP), col=cols,pch=20,cex=0.8)

#####Fisher test to find differences by clones in normal vs. tumor
p_value=NULL
id<-match(colnames(clone_type_filter_IG),rownames(PAAD.repertoire.diversity))
for(i in 1:dim(clone_type_filter_IG)[1]){
  print(i)
  tab<-table(clone_type_filter_IG[i,],PAAD.repertoire.diversity$Tumor_type_2categ[id])
  p_value[i]=fisher.test(tab)$p.value
}
clone_type_filter_IG_sign<-clone_type_filter_IG[which(p_value<0.05),]
annotation_row = data.frame(PAAD.repertoire.diversity$Tumor_type_2categ[id])
colnames(annotation_row)<-"Tumor_type_2categ"
rownames(annotation_row)<-rownames(PAAD.repertoire.diversity)[id]
ann_colors = list (Tumor_type_2categ = c("normal_pseudonormal_pancreas" = brewer.pal(3,"Accent")[1],
                                         "Tumor_pancreas"= brewer.pal(3,"Accent")[2]))

tiff(paste0("Results/heatmap_common_clones_IG_sign.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IG_sign),border_color=F,show_rownames = F, annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols[6]),breaks = c(0,0.9,1))
dev.off()

##TCR
p_value=NULL
id<-match(colnames(clone_type_filter_TCR),rownames(PAAD.repertoire.diversity))
for(i in 1:dim(clone_type_filter_TCR)[1]){
  print(i)
  tab<-table(clone_type_filter_TCR[i,],PAAD.repertoire.diversity$Tumor_type_2categ[id])
  p_value[i]=fisher.test(tab)$p.value
}
clone_type_filter_TCR_sign<-clone_type_filter_TCR[which(p_value<0.05),]
annotation_row = data.frame(PAAD.repertoire.diversity$Tumor_type_2categ[id])
colnames(annotation_row)<-"Tumor_type_2categ"
rownames(annotation_row)<-rownames(PAAD.repertoire.diversity)[id]
ann_colors = list (Tumor_type_2categ = c("normal_pseudonormal_pancreas" = brewer.pal(3,"Accent")[1],
                                         "Tumor_pancreas"= brewer.pal(3,"Accent")[2]))

tiff(paste0("Results/heatmap_common_clones_TCR_sign.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_TCR_sign),border_color=F,show_rownames = F, annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols[6]),breaks = c(0,0.9,1))
dev.off()

############
### 2. Normalization by transforming to relative abundance
############
id<-match(rownames(clone_type_filter_IG),colnames(clone_type_IG))
clone_type_abundance_IG<-clone_type_IG[,id]
clone_type_relative<-t(apply(clone_type_abundance_IG,1,function (x) x/sum(x))) #180/2898
p_value=NULL
for(i in 1:dim(clone_type_relative)[2]){
  print(i)
  mod<-glm(PAAD.repertoire.diversity$Tumor_type_2categ~clone_type_relative[,i],family = "binomial")
  p_value[i]=coefficients(summary(mod))[2,4]
}

clone_type_relative_sign<-clone_type_relative[,which(p_value<0.05)]

##Plot results
annotation_row = data.frame(PAAD.repertoire.diversity$Tumor_type_2categ)
colnames(annotation_row)<-"Tumor_type_2categ"
rownames(annotation_row)<-rownames(PAAD.repertoire.diversity)
ann_colors = list (Tumor_type_2categ = c("normal_pseudonormal_pancreas" = brewer.pal(3,"Accent")[1],
                                         "Tumor_pancreas"= brewer.pal(3,"Accent")[2]))

tiff("Results/heatmap_relativeAbundance_Ig_sign.tiff",width = 5000, height = 3000, res = 300)
pheatmap(clone_type_relative_sign,scale = "column",border_color=F,show_rownames = F, annotation_row = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))
dev.off()

clone_type_relative_sign_order<-clone_type_relative_sign[order(PAAD.repertoire.diversity$Tumor_type_2categ),]
df_long <- melt(clone_type_relative_sign_order, id.vars = "Sample", variable.name = "Clones")
colnames(df_long)<-c("Sample","Clones","value")
library(randomcoloR)
n <- 29
palette <- distinctColorPalette(n)
tiff("Results/barplot_relativeAbundance_Ig_sign.tiff",width = 5000, height = 3000, res = 300)
ggplot(df_long, aes(x = Sample, y = value, fill = Clones)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=palette)
dev.off()



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
