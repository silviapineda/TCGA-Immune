rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. V gene usage analysis
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: V gene usage analysis
###         
###
### Author: Silvia Pineda
### Date: July, 2019
############################################################################################
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(glmnet)
library(factoextra)

##Colors: brewer.pal(4, "Accent")
##"#BEAED4" (IGH) "#7FC97F" (IGK) "#FDC086" (IGL)
setwd("~/TCGA-Immune/")

vgene_counts<-readRDS("Data/PAAD/vgene_counts_all_cancers.rds")##Kat
vgene_counts<-vgene_counts[which(vgene_counts$type=="Primary_Solid_Tumor"),]
#write.csv(vgene_counts,"Data/PAAD/vgene_counts_all_cancers.csv")
vgene_counts_IGH<-vgene_counts[,which(grepl('^IGH', colnames(vgene_counts))==T)]
vgene_counts_IGK<-vgene_counts[,which(grepl('^IGK', colnames(vgene_counts))==T)]
vgene_counts_IGL<-vgene_counts[,which(grepl('^IGL', colnames(vgene_counts))==T)]
vgene_counts_TRA<-vgene_counts[,which(grepl('^TRA', colnames(vgene_counts))==T)]
vgene_counts_TRB<-vgene_counts[,which(grepl('^TRB', colnames(vgene_counts))==T)]

vgene_counts_IG<-vgene_counts[,which(grepl('^IG', colnames(vgene_counts))==T)]

##Filter by vgenes share at least in 2% of the samples (4569*0.02 = 91) or 2 samples
vgene_counts_IGH_filter<-vgene_counts_IGH[,which(colSums(vgene_counts_IGH!=0)>=2)] #
vgene_counts_IGK_filter<-vgene_counts_IGK[,which(colSums(vgene_counts_IGK!=0)>=2)] #
vgene_counts_IGL_filter<-vgene_counts_IGL[,which(colSums(vgene_counts_IGL!=0)>=2)] #
vgene_counts_TRA_filter<-vgene_counts_TRA[,which(colSums(vgene_counts_TRA!=0)>=2)] #
vgene_counts_TRB_filter<-vgene_counts_TRB[,which(colSums(vgene_counts_TRB!=0)>=2)] #

vgene_counts_IG_filter<-vgene_counts_IG[,which(colSums(vgene_counts_IG!=0)>=2)] #


#### 
vgene_counts_IG_zerosubs<-vgene_counts_IGH_filter+1
z_vgene_IG<-log(vgene_counts_IG_zerosubs)
clrx_vgene_IG <- apply(z_vgene_IG, 2, function(x) x - rowMeans(z_vgene_IG))

annotation_row<-data.frame(vgene_counts$cancer)
rownames(annotation_row)<-rownames(clrx_vgene_IG)

tiff("Pan-cancer/heatmap_clrx_data.tiff",res=300,h=3000,w=3500)
pheatmap(t(clrx_vgene_IG),scale = "column",border_color=F,show_colnames = F, annotation_col = annotation_row)
dev.off()

pca <- prcomp(clrx_vgene_IG, scale = TRUE)
pc <- c(1,2)
groups <- as.factor(vgene_counts$cancer)
# Define shapes
#shapes = c(0:25,0:4)
#shapes <- shapes[as.numeric(groups)]
#library(RColorBrewer)
#n <- 40
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#col_vector<-col_vector[-c(4,27)]
tiff("Pan-cancer/PCA_clrx_vgene_IGH.tiff",res=300,h=2000,w=2500)
fviz_pca_ind(pca,
             geom.ind = "point",# show points only (nbut not "text")
             #pointshape = 21,
             fill.ind = groups,
             col.ind = groups, # color by groups
             #palette = col_vector,
             #addEllipses = TRUE, # Concentration ellipses
             legend.title = "Tumor types",
             repel = TRUE)

dev.off()

id_COAD<-which(vgene_counts$cancer=="COAD")
id_UCEC<-which(vgene_counts$cancer=="UCEC")
clrx_vgene_IG_COAD<-clrx_vgene_IG[id_COAD,]
clrx_vgene_IG_UCEC<-clrx_vgene_IG[id_UCEC,]
clrx_vgene_IG_Others<-clrx_vgene_IG[-c(id_COAD,id_UCEC),]
data<-rbind(clrx_vgene_IG_COAD,clrx_vgene_IG_UCEC,clrx_vgene_IG_Others)
group<-c(rep("COAD",nrow(clrx_vgene_IG_COAD)),rep("UCEC",nrow(clrx_vgene_IG_UCEC)),rep("Other",nrow(clrx_vgene_IG_Others)))
annotation_row<-data.frame(group)
rownames(annotation_row)<-rownames(data)

tiff("Pan-cancer/heatmap_clrx_data_COAD_UCEC.tiff",res=300,h=3000,w=3500)
pheatmap(t(data),scale="column",border_color=F,show_colnames = F, annotation_col = annotation_row)
dev.off()

pca <- prcomp(data, scale = TRUE)
pc <- c(1,2)
groups <- as.factor(group)
tiff("Pan-cancer/PCA_clrx_vgene_COAD_UCEC_IGH.tiff",res=300,h=2000,w=2500)
fviz_pca_ind(pca,
             geom.ind = "point",# show points only (nbut not "text")
             #pointshape = 21,
             fill.ind = groups,
             col.ind = groups, # color by groups
             #palette = col_vector,
             #addEllipses = TRUE, # Concentration ellipses
             legend.title = "Tumor types",
             repel = TRUE)
dev.off()



#####################
#### This is the data from Kat:
#############
vgenes_usage<-readRDS("Data/PAAD/vgene_usage_all_cancers.rds")
#vgenes_usage<-readRDS("Data/PAAD/vgene_usage_all_cancers_silvia.rds")
vgenes_usage_IG<-vgenes_usage[,which(grepl('^IGH', colnames(vgenes_usage))==T)]

##Filter by vgenes share at least in 20% of the samples (6234*0.2 = 91) or 2 samples
vgenes_usage_filter<-vgenes_usage_IG[,which(colSums(vgenes_usage_IG!=0)>=(dim(vgenes_usage_IG)[1]*0.2))] #67
vgenes_usage_filter_2<-apply(as.matrix.noquote(vgenes_usage_filter),2,as.numeric)
rownames(vgenes_usage_filter_2)<-rownames(vgenes_usage_filter)
pca <- prcomp(vgenes_usage_filter_2, scale = TRUE)
pc <- c(1,2)
groups <- as.factor(vgenes_usage$cancer)
tiff("Pan-cancer/PCA_Katdata.tiff",res=300,h=2000,w=2500)
fviz_pca_ind(pca,
             geom.ind = "point", # show points only (nbut not "text")
             #pointshape = 21,
             fill.ind = groups,
             col.ind = groups, # color by groups
             #palette = c("#7FC97F","#BEAED4","#FDC086","#B3CDE3"),
             #addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups",
             repel = TRUE)
dev.off()

annotation_row<-data.frame(vgenes_usage$cancer)
rownames(annotation_row)<-rownames(vgenes_usage_filter)
###Heatmap with alltumor types
tiff("Pan-cancer/heatmap_Katdata_all.tiff",res=300,h=3000,w=3500)
pheatmap(t(vgenes_usage_filter_2),scale="column",border_color=F,show_colnames = F, annotation_col = annotation_row)
dev.off()

vgenes_usage_filter_COAD<-vgenes_usage_filter[which(vgenes_usage$cancer=="COAD"),]
vgenes_usage_filter_UCEC<-vgenes_usage_filter[which(vgenes_usage$cancer=="UCEC"),]
vgenes_usage_filter_Others<-vgenes_usage_filter[which(vgenes_usage$cancer!="UCEC" | vgenes_usage$cancer=="COAD"),]
data<-rbind(vgenes_usage_filter_COAD,vgenes_usage_filter_UCEC,vgenes_usage_filter_Others)
group<-c(rep("COAD",nrow(vgenes_usage_filter_COAD)),rep("UCEC",nrow(vgenes_usage_filter_UCEC)),rep("Other",nrow(vgenes_usage_filter_Others)))
annotation_row<-data.frame(group)
rownames(annotation_row)<-rownames(data)

tiff("Pan-cancer/heatmap_Katdata_COAD_UCEC.tiff",res=300,h=3000,w=3500)
pheatmap(t(vgenes_usage_filter),scale="column",border_color=F,show_colnames = F, annotation_col = annotation_row)
dev.off()


###########
####Number of clones  that uses Vgenes
##########
vgene_counts<-readRDS("Data/PAAD/vgene_num_clones_all_cancers.rds")
vgene_counts<-vgene_counts[which(vgene_counts$tumor_type=="Primary_Solid_Tumor"),]
vgene_counts<-vgene_counts[which(vgene_counts$cancer %in% c("BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PRAD","THCA","UCEC")),]
#write.csv(vgene_counts,"Data/PAAD/vgene_counts_all_cancers.csv")
vgene_counts_IGH<-vgene_counts[,which(grepl('^IGH', colnames(vgene_counts))==T)]
vgene_counts_IGK<-vgene_counts[,which(grepl('^IGK', colnames(vgene_counts))==T)]
vgene_counts_IGL<-vgene_counts[,which(grepl('^IGL', colnames(vgene_counts))==T)]
vgene_counts_TRA<-vgene_counts[,which(grepl('^TRA', colnames(vgene_counts))==T)]
vgene_counts_TRB<-vgene_counts[,which(grepl('^TRB', colnames(vgene_counts))==T)]

vgene_counts_IG<-vgene_counts[,which(grepl('^IG', colnames(vgene_counts))==T)]

##Filter by vgenes share at least in 2% of the samples (4569*0.02 = 91) or 2 samples
vgene_counts_IGH_filter<-vgene_counts_IGH[,which(colSums(vgene_counts_IGH!=0)>=2)] #
vgene_counts_IGK_filter<-vgene_counts_IGK[,which(colSums(vgene_counts_IGK!=0)>=2)] #
vgene_counts_IGL_filter<-vgene_counts_IGL[,which(colSums(vgene_counts_IGL!=0)>=2)] #
vgene_counts_TRA_filter<-vgene_counts_TRA[,which(colSums(vgene_counts_TRA!=0)>=2)] #
vgene_counts_TRB_filter<-vgene_counts_TRB[,which(colSums(vgene_counts_TRB!=0)>=2)] #

vgene_counts_IG_filter<-vgene_counts_IG[,which(colSums(vgene_counts_IG!=0)>=2)] #


#### 
vgene_counts_IG_zerosubs<-vgene_counts_IGH_filter+1
z_vgene_IG<-log(vgene_counts_IG_zerosubs)
clrx_vgene_IG <- apply(z_vgene_IG, 2, function(x) x - rowMeans(z_vgene_IG))

annotation_row<-data.frame(vgene_counts$cancer)
rownames(annotation_row)<-rownames(clrx_vgene_IG)

tiff("Pan-cancer/heatmap_clrx_Kat_data.tiff",res=300,h=3000,w=3500)
pheatmap(t(clrx_vgene_IG),scale = "column",border_color=F,show_colnames = F, annotation_col = annotation_row)
dev.off()

pca <- prcomp(clrx_vgene_IG, scale = TRUE)
pc <- c(1,2)
groups <- as.factor(vgene_counts$cancer)
# Define shapes
#shapes = c(0:25,0:4)
#shapes <- shapes[as.numeric(groups)]
#library(RColorBrewer)
#n <- 40
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#col_vector<-col_vector[-c(4,27)]
tiff("Pan-cancer/PCA_clrx_Kat_vgene_IGH.tiff",res=300,h=2000,w=2500)
fviz_pca_ind(pca,
             geom.ind = "point",# show points only (nbut not "text")
             #pointshape = 21,
             fill.ind = groups,
             col.ind = groups, # color by groups
             #palette = col_vector,
             #addEllipses = TRUE, # Concentration ellipses
             legend.title = "Tumor types",
             repel = TRUE)

dev.off()

id_COAD<-which(vgene_counts$cancer=="COAD")
id_UCEC<-which(vgene_counts$cancer=="UCEC")
clrx_vgene_IG_COAD<-clrx_vgene_IG[id_COAD,]
clrx_vgene_IG_UCEC<-clrx_vgene_IG[id_UCEC,]
clrx_vgene_IG_Others<-clrx_vgene_IG[-c(id_COAD,id_UCEC),]
data<-rbind(clrx_vgene_IG_COAD,clrx_vgene_IG_UCEC,clrx_vgene_IG_Others)
group<-c(rep("COAD",nrow(clrx_vgene_IG_COAD)),rep("UCEC",nrow(clrx_vgene_IG_UCEC)),rep("Other",nrow(clrx_vgene_IG_Others)))
annotation_row<-data.frame(group)
rownames(annotation_row)<-rownames(data)

tiff("Pan-cancer/heatmap_clrx_Kat_data_COAD_UCEC.tiff",res=300,h=3000,w=3500)
pheatmap(t(data),scale="column",border_color=F,show_colnames = F, annotation_col = annotation_row)
dev.off()

pca <- prcomp(data, scale = TRUE)
pc <- c(1,2)
groups <- as.factor(group)
tiff("Pan-cancer/PCA_clrx_Katdata_vgene_COAD_UCEC_IGH.tiff",res=300,h=2000,w=2500)
fviz_pca_ind(pca,
             geom.ind = "point",# show points only (nbut not "text")
             #pointshape = 21,
             fill.ind = groups,
             col.ind = groups, # color by groups
             #palette = col_vector,
             #addEllipses = TRUE, # Concentration ellipses
             legend.title = "Tumor types",
             repel = TRUE)
dev.off()
