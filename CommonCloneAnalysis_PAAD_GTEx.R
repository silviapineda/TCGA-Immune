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

load("Data/PAAD_GTEx_ValTumor_ValNormal/PAAD_GTEx_ValTumor_ValNormal_FullData.Rdata")

###################################################
## Common clones across ALL samples
#################################################
Repertoire.Diversity<-PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity[which(PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome=="normal_pancreas (GTEx)"|
                                                                      PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome== "PDAC (TCGA)" |
                                                                        PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome== "PDAC (Val)" |  
                                                                PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome== "normal_pancreas (Val)"),]
Repertoire.Diversity$outcome<-factor(Repertoire.Diversity$outcome)
# Rename all levels, by name
levels(Repertoire.Diversity$outcome) <- list("GTEX-Normal"="normal_pancreas (GTEx)", "TCGA-PDAC"= "PDAC (TCGA)", "Validation-PDAC"="PDAC (Val)","Validation-Normal"="normal_pancreas (Val)")
id<-match(data_merge$sample,rownames(Repertoire.Diversity))
data_merge_filter<-data_merge[which(is.na(id)==F),]

#chain=c("IGHV","IGKV","IGLV")
############
## 1. Build the matrix with the clones by samples
data_chain_IGK<-data_merge_filter[which(data_merge_filter$chainType=="IGK"),]
clone_type_IGK<-t(as.data.frame(unclass(table(data_chain_IGK$V_J_lenghCDR3_CloneId,factor(data_chain_IGK$sample))))) 
data_chain_IGL<-data_merge_filter[which(data_merge_filter$chainType=="IGL"),]
clone_type_IGL<-t(as.data.frame(unclass(table(data_chain_IGL$V_J_lenghCDR3_CloneId,factor(data_chain_IGL$sample))))) 
data_chain_IGH<-data_merge_filter[which(data_merge_filter$chainType=="IGH"),]
clone_type_IGH<-t(as.data.frame(unclass(table(data_chain_IGH$V_J_lenghCDR3_CloneId,factor(data_chain_IGH$sample))))) 

data_chain_TRA<-data_merge_filter[which(data_merge_filter$chainType=="TRA"),]
clone_type_TRA<-t(as.data.frame(unclass(table(data_chain_TRA$V_J_lenghCDR3_CloneId,factor(data_chain_TRA$sample))))) 
data_chain_TRB<-data_merge_filter[which(data_merge_filter$chainType=="TRB"),]
clone_type_TRB<-t(as.data.frame(unclass(table(data_chain_TRB$V_J_lenghCDR3_CloneId,factor(data_chain_TRB$sample))))) 
data_chain_TRD<-data_merge_filter[which(data_merge_filter$chainType=="TRD"),]
clone_type_TRD<-t(as.data.frame(unclass(table(data_chain_TRD$V_J_lenghCDR3_CloneId,factor(data_chain_TRD$sample))))) 
data_chain_TRG<-data_merge_filter[which(data_merge_filter$chainType=="TRG"),]
clone_type_TRG<-t(as.data.frame(unclass(table(data_chain_TRG$V_J_lenghCDR3_CloneId,factor(data_chain_TRG$sample))))) 

##Build present vs no present
clone_type_presence_IGK<-apply(clone_type_IGK,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGL<-apply(clone_type_IGL,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGH<-apply(clone_type_IGH,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRA<-apply(clone_type_TRA,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRB<-apply(clone_type_TRB,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRD<-apply(clone_type_TRD,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRG<-apply(clone_type_TRG,1,function(x) ifelse(x==0,0,1))

###Filter by clones share at least by 2 samples
clone_type_filter_IGK<-clone_type_presence_IGK[which(rowSums(clone_type_presence_IGK)>2),] #
clone_type_filter_IGL<-clone_type_presence_IGL[which(rowSums(clone_type_presence_IGL)>2),] #
clone_type_filter_IGH<-clone_type_presence_IGH[which(rowSums(clone_type_presence_IGH)>2),] #
clone_type_filter_TRA<-clone_type_presence_TRA[which(rowSums(clone_type_presence_TRA)>2),] #
clone_type_filter_TRB<-clone_type_presence_TRB[which(rowSums(clone_type_presence_TRB)>2),] #
clone_type_filter_TRD<-clone_type_presence_TRD[which(rowSums(clone_type_presence_TRD)>2),] #
clone_type_filter_TRG<-clone_type_presence_TRG[which(rowSums(clone_type_presence_TRG)>2),] #

##Heatmp

#IGK
brewer.pal(4,name = "Accent")
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")
id<-match(colnames(clone_type_filter_IGK),Repertoire.Diversity$sample)
annotation_col = data.frame(Repertoire.Diversity$outcome[id])
table(Repertoire.Diversity$outcome)
ann_colors = list (outcome = c("GTEX-Normal" = cols[1],
                               "TCGA-PDAC" = cols[2],
                               "Validation-Normal"= cols[3],
                               "Validation-PDAC" = cols[4]))
colnames(annotation_col)<-"outcome"
rownames(annotation_col)<-colnames(clone_type_filter_IGK)

cols2 = brewer.pal(12,name = "Paired")

tiff(paste0("Results/CommonClones//heatmap_common_clones_IGK.tiff"),width = 2000, height = 1500, res = 300)
pheatmap(clone_type_filter_IGK,border_color=F,show_rownames = F, show_colnames = F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = c("white",cols2[2]),breaks = c(0,0.9,1))
dev.off()


#IGH
brewer.pal(4,name = "Accent")
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")
id<-match(colnames(clone_type_filter_IGH),Repertoire.Diversity$sample)
annotation_col = data.frame(Repertoire.Diversity$outcome[id])
table(Repertoire.Diversity$outcome)
ann_colors = list (outcome = c("GTEX-Normal" = cols[1],
                               "TCGA-PDAC" = cols[2],
                               "Validation-Normal"= cols[3],
                               "Validation-PDAC" = cols[4]))
colnames(annotation_col)<-"outcome"
rownames(annotation_col)<-colnames(clone_type_filter_IGH)

cols2 = brewer.pal(12,name = "Paired")

tiff(paste0("Results/CommonClones//heatmap_common_clones_IGH.tiff"),width = 2000, height = 1500, res = 300)
pheatmap(clone_type_filter_IGH,border_color=F,show_rownames = F, show_colnames = F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = c("white",cols2[4]),breaks = c(0,0.9,1))
dev.off()

#IGL
brewer.pal(4,name = "Accent")
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")
id<-match(colnames(clone_type_filter_IGL),Repertoire.Diversity$sample)
annotation_col = data.frame(Repertoire.Diversity$outcome[id])
table(Repertoire.Diversity$outcome)
ann_colors = list (outcome = c("GTEX-Normal" = cols[1],
                               "TCGA-PDAC" = cols[2],
                               "Validation-Normal"= cols[3],
                               "Validation-PDAC" = cols[4]))
colnames(annotation_col)<-"outcome"
rownames(annotation_col)<-colnames(clone_type_filter_IGL)

cols2 = brewer.pal(12,name = "Paired")

tiff(paste0("Results/CommonClones//heatmap_common_clones_IGL.tiff"),width = 2000, height = 1500, res = 300)
pheatmap(clone_type_filter_IGL,border_color=F,show_rownames = F, show_colnames = F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = c("white",cols2[6]),breaks = c(0,0.9,1))
dev.off()

#TRA
brewer.pal(4,name = "Accent")
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")
id<-match(colnames(clone_type_filter_TRA),Repertoire.Diversity$sample)
annotation_col = data.frame(Repertoire.Diversity$outcome[id])
table(Repertoire.Diversity$outcome)
ann_colors = list (outcome = c("GTEX-Normal" = cols[1],
                               "TCGA-PDAC" = cols[2],
                               "Validation-Normal"= cols[3],
                               "Validation-PDAC" = cols[4]))
colnames(annotation_col)<-"outcome"
rownames(annotation_col)<-colnames(clone_type_filter_TRA)

cols2 = brewer.pal(12,name = "Paired")

tiff(paste0("Results/CommonClones//heatmap_common_clones_TRA.tiff"),width = 2000, height = 1500, res = 300)
pheatmap(clone_type_filter_TRA,border_color=F,show_rownames = F, show_colnames = F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = c("white",cols2[8]),breaks = c(0,0.9,1))
dev.off()


#TRB
brewer.pal(4,name = "Accent")
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")
id<-match(colnames(clone_type_filter_TRB),Repertoire.Diversity$sample)
annotation_col = data.frame(Repertoire.Diversity$outcome[id])
table(Repertoire.Diversity$outcome)
ann_colors = list (outcome = c("GTEX-Normal" = cols[1],
                               "TCGA-PDAC" = cols[2],
                               "Validation-Normal"= cols[3],
                               "Validation-PDAC" = cols[4]))
colnames(annotation_col)<-"outcome"
rownames(annotation_col)<-colnames(clone_type_filter_TRB)

cols2 = brewer.pal(12,name = "Paired")

tiff(paste0("Results/CommonClones//heatmap_common_clones_TRB.tiff"),width = 2000, height = 1500, res = 300)
pheatmap(clone_type_filter_TRB,border_color=F,show_rownames = F, show_colnames = F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = c("white",cols2[10]),breaks = c(0,0.9,1))
dev.off()



########################################
## GTEx vs. TCGA Analysis ##############
########################################
#AAD.GTEx.repertoire.diversity.tumor.normmal<-PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity[which(PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome=="normal_pancreas (GTEx)"|
#                                                                                                        PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome== "PDAC (TCGA)"),]
#PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome<-factor(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome)



##########################
#####Fisher test to find differences by clones in normal vs. tumor
##########################
#IGH
p_value=NULL
for(i in 1:dim(clone_type_filter_IGH)[1]){
  print(i)
  id<-match(colnames(clone_type_filter_IGH),rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal))
  tab<-table(clone_type_filter_IGH[i,],PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
  p_value[i]=fisher.test(tab)$p.value
}
p.adj<-p.adjust(p_value,"fdr")
clone_type_filter_IGH_sign<-clone_type_filter_IGH[which(p.adj<0.05),] 
annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal)[id]
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F", "#BEAED4")
ann_colors = list (outcome = c("normal_pancreas (GTEx)" = cols[1],
                               "PDAC (TCGA)" = cols[2]))

tiff(paste0("Results/CommonClones//heatmap_common_clones_IGH_sign.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGH_sign),border_color=F,show_rownames = F, annotation_row = annotation_row,
         annotation_colors = ann_colors,color =c("white",cols2[2]),breaks = c(0,0.9,1))
dev.off()

#IGK
p_value=NULL
for(i in 1:dim(clone_type_filter_IGK)[1]){
  print(i)
  id<-match(colnames(clone_type_filter_IGK),rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal))
  tab<-table(clone_type_filter_IGK[i,],PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
  p_value[i]=fisher.test(tab)$p.value
}
p.adj<-p.adjust(p_value,"fdr")
clone_type_filter_IGK_sign<-clone_type_filter_IGK[which(p.adj<0.05),]
annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal)[id]
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F", "#BEAED4")
ann_colors = list (outcome = c("normal_pancreas (GTEx)" = cols[1],
                                         "PDAC (TCGA)" = cols[2]))

tiff(paste0("Results/CommonClones//heatmap_common_clones_IGK_sign.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGK_sign),border_color=F,show_rownames = F, annotation_row = annotation_row,show_colnames = F,
         annotation_colors = ann_colors,color =c("white",cols2[4]),breaks = c(0,0.9,1))
dev.off()

#IGL
p_value=NULL
for(i in 1:dim(clone_type_filter_IGL)[1]){
  print(i)
  id<-match(colnames(clone_type_filter_IGL),rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal))
  tab<-table(clone_type_filter_IGL[i,],PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
  p_value[i]=fisher.test(tab)$p.value
}
p.adj<-p.adjust(p_value,"fdr")
clone_type_filter_IGL_sign<-clone_type_filter_IGL[which(p.adj<0.05),] 
annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal)[id]
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F", "#BEAED4")
ann_colors = list (outcome = c("normal_pancreas (GTEx)" = cols[1],
                               "PDAC (TCGA)" = cols[2]))

tiff(paste0("Results/CommonClones//heatmap_common_clones_IGL_sign.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGL_sign),border_color=F,show_rownames = F, annotation_row = annotation_row,show_colnames = F,
         annotation_colors = ann_colors,color =c("white",cols2[6]),breaks = c(0,0.9,1))
dev.off()

#TRA
p_value=NULL
for(i in 1:dim(clone_type_filter_TRA)[1]){
  print(i)
  id<-match(colnames(clone_type_filter_TRA),rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal))
  tab<-table(clone_type_filter_TRA[i,],PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
  p_value[i]=fisher.test(tab)$p.value
}
p.adj<-p.adjust(p_value,"fdr")
clone_type_filter_TRA_sign<-clone_type_filter_TRA[which(p.adj<0.05),] 
rownames(clone_type_filter_TRA)[which(p.adj<0.05)]

#TRB
p_value=NULL
for(i in 1:dim(clone_type_filter_TRB)[1]){
  print(i)
  id<-match(colnames(clone_type_filter_TRB),rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal))
  tab<-table(clone_type_filter_TRB[i,],PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
  p_value[i]=fisher.test(tab)$p.value
}
p.adj<-p.adjust(p_value,"fdr")
clone_type_filter_TRB_sign<-clone_type_filter_TRB[which(p.adj<0.05),]
annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal)[id]
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F", "#BEAED4")
ann_colors = list (outcome = c("normal_pancreas (GTEx)" = cols[1],
                               "PDAC (TCGA)" = cols[2]))

tiff(paste0("Results/CommonClones//heatmap_common_clones_TRB_sign.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_TRB_sign),border_color=F,show_rownames = F, annotation_row = annotation_row,
         annotation_colors = ann_colors,color =c("white",cols2[10]),breaks = c(0,0.9,1))
dev.off()

#TRD
p_value=NULL
for(i in 1:dim(clone_type_filter_TRD)[1]){
  print(i)
  id<-match(colnames(clone_type_filter_TRD),rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal))
  tab<-table(clone_type_filter_TRD[i,],PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
  p_value[i]=fisher.test(tab)$p.value
}
p.adj<-p.adjust(p_value,"fdr")
clone_type_filter_TRD_sign<-clone_type_filter_TRD[which(p.adj<0.05),] 
rownames(clone_type_filter_TRD)[which(p.adj<0.05)]#"TRDV2_TRDJ1_18_1" (present: 3 GTEx and 0 TCGA 130 ; No-present: 4 GTEx and 27 TCGA)
annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
table(clone_type_filter_TRD_sign,PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])

#TRG
p_value=NULL
for(i in 1:dim(clone_type_filter_TRG)[1]){
  print(i)
  id<-match(colnames(clone_type_filter_TRG),rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal))
  tab<-table(clone_type_filter_TRG[i,],PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
  p_value[i]=fisher.test(tab)$p.value
}
p.adj<-p.adjust(p_value,"fdr")
clone_type_filter_TRG_sign<-clone_type_filter_TRG[which(p.adj<0.05),] 
rownames(clone_type_filter_TRG)[which(p.adj<0.05)]


########################################
## All Common Analysis ##############
########################################
PAAD.Val.repertoire.diversity.tumor.normmal<-PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity[which(PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome=="PDAC (Val)"|
                                                                                                        PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome== "PDAC (TCGA)" |
                                                                                                       PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome=="normal_pancreas (GTEx)"|
                                                                                                       PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome== "PDAC (TCGA)" |
                                                                                                 PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome=="normal_pancreas (Val)"),]
PAAD.Val.repertoire.diversity.tumor.normmal$outcome<-factor(PAAD.Val.repertoire.diversity.tumor.normmal$outcome)


## Common clones across samples

id<-match(data_merge$sample,rownames(PAAD.Val.repertoire.diversity.tumor.normmal))
data_merge<-data_merge[which(is.na(id)==F),]

#chain=c("IGHV","IGKV","IGLV")

## 1. Build the matrix with the clones by samples
data_chain_IGK<-data_merge[which(data_merge$chainType=="IGK"),]
clone_type_IGK<-t(as.data.frame(unclass(table(data_chain_IGK$V_J_lenghCDR3_CloneId,factor(data_chain_IGK$sample))))) 
data_chain_IGL<-data_merge[which(data_merge$chainType=="IGL"),]
clone_type_IGL<-t(as.data.frame(unclass(table(data_chain_IGL$V_J_lenghCDR3_CloneId,factor(data_chain_IGL$sample))))) 
data_chain_IGH<-data_merge[which(data_merge$chainType=="IGH"),]
clone_type_IGH<-t(as.data.frame(unclass(table(data_chain_IGH$V_J_lenghCDR3_CloneId,factor(data_chain_IGH$sample))))) 

data_chain_TRA<-data_merge[which(data_merge$chainType=="TRA"),]
clone_type_TRA<-t(as.data.frame(unclass(table(data_chain_TRA$V_J_lenghCDR3_CloneId,factor(data_chain_TRA$sample))))) 
data_chain_TRB<-data_merge[which(data_merge$chainType=="TRB"),]
clone_type_TRB<-t(as.data.frame(unclass(table(data_chain_TRB$V_J_lenghCDR3_CloneId,factor(data_chain_TRB$sample))))) 
data_chain_TRD<-data_merge[which(data_merge$chainType=="TRD"),]
clone_type_TRD<-t(as.data.frame(unclass(table(data_chain_TRD$V_J_lenghCDR3_CloneId,factor(data_chain_TRD$sample))))) 
data_chain_TRG<-data_merge[which(data_merge$chainType=="TRG"),]
clone_type_TRG<-t(as.data.frame(unclass(table(data_chain_TRG$V_J_lenghCDR3_CloneId,factor(data_chain_TRG$sample))))) 

##Build present vs no present
clone_type_presence_IGK<-apply(clone_type_IGK,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGL<-apply(clone_type_IGL,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGH<-apply(clone_type_IGH,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRA<-apply(clone_type_TRA,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRB<-apply(clone_type_TRB,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRD<-apply(clone_type_TRD,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRG<-apply(clone_type_TRG,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRDG<-apply(clone_type_TRDG,1,function(x) ifelse(x==0,0,1))

###Filter by clones share at least in 5% of the samples (311*0.05= 15) or 2 samples
clone_type_filter_IGK<-clone_type_presence_IGK[which(rowSums(clone_type_presence_IGK)>2),] #
clone_type_filter_IGL<-clone_type_presence_IGL[which(rowSums(clone_type_presence_IGL)>2),] #
clone_type_filter_IGH<-clone_type_presence_IGH[which(rowSums(clone_type_presence_IGH)>2),] #
clone_type_filter_TRA<-clone_type_presence_TRA[which(rowSums(clone_type_presence_TRA)>2),] #
clone_type_filter_TRB<-clone_type_presence_TRB[which(rowSums(clone_type_presence_TRB)>2),] #
clone_type_filter_TRD<-clone_type_presence_TRD[which(rowSums(clone_type_presence_TRD)>2),] #
clone_type_filter_TRG<-clone_type_presence_TRG[which(rowSums(clone_type_presence_TRG)>2),] #

##Heatmp
brewer.pal(4,name = "Accent")
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")
id<-match(colnames(clone_type_filter_IGH),PAAD.Val.repertoire.diversity.tumor.normmal$sample)
annotation_row = data.frame(PAAD.Val.repertoire.diversity.tumor.normmal$outcome[id])
ann_colors = list (outcome = c("normal_pancreas (GTEx)" = cols[1],
                               "PDAC (TCGA)" = cols[2],
                               "normal_pancreas (Val)"= cols[3],
                               "PDAC (Val)" = cols[4]))
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(PAAD.Val.repertoire.diversity.tumor.normmal)[id]

cols2 = brewer.pal(12,name = "Paired")

#IGH
tiff(paste0("Results/CommonClones//heatmap_common_clones_IGH_All.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGH),border_color=F,show_rownames = F, show_colnames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols2[2]),breaks = c(0,0.9,1))
dev.off()

#"IGHV3-21_IGHJ4_12_656"  "IGHV3-23_IGHJ4_11_242"  "IGHV3-23_IGHJ4_14_1055" "IGHV3-7_IGHJ4_13_265"   "IGHV3-9_IGHJ6_13_34"  
id2<-match(rownames(clone_type_filter_IGH_sign),rownames(clone_type_filter_IGH))
table(clone_type_filter_IGH[id2[2],],PAAD.Val.repertoire.diversity.tumor.normmal$outcome[id])


#IGK
tiff(paste0("Results/CommonClones/heatmap_common_clones_IGK_All.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGK),show_rownames = F, show_colnames = F,border_color=F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols2[4]),breaks = c(0,0.9,1))
dev.off()

id2<-match(rownames(clone_type_filter_IGK_sign),rownames(clone_type_filter_IGK))
table(clone_type_filter_IGK[id2,],PAAD.Val.repertoire.diversity.tumor.normmal$outcome[id])


#IGL
tiff(paste0("Results/CommonClones/heatmap_common_clones_IGL_All.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGL),border_color=F,show_rownames = F, show_colnames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols2[6]),breaks = c(0,0.9,1))
dev.off()

#TRA
tiff(paste0("Results/CommonClones/heatmap_common_clones_TRA_All.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_TRA),border_color=F,show_rownames = F, show_colnames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols2[8]),breaks = c(0,0.9,1))
dev.off()

#TRB
tiff(paste0("Results/CommonClones/heatmap_common_clones_TRB_All.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_TRB),border_color=F,show_rownames = F,annotation_row = annotation_row,show_colnames = F,
         annotation_colors = ann_colors,color = c("white",cols2[10]),breaks = c(0,0.9,1))

dev.off()

#TRD
tiff(paste0("Results/CommonClones/heatmap_common_clones_TRD_All.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_TRD),border_color=F,show_rownames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols2[12]),breaks = c(0,0.9,1))

dev.off()
id2<-match(names(clone_type_filter_TRD_sign),colnames(clone_type_filter_TRD))
id<-match(colnames(clone_type_filter_TRD)[id2],PAAD.Val.repertoire.diversity.tumor.normmal$sample)
table(clone_type_filter_TRD[id2],PAAD.Val.repertoire.diversity.tumor.normmal$outcome[id])

#TRG
tiff(paste0("Results/CommonClones/heatmap_common_clones_TRG_All.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_TRG),border_color=F,show_rownames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols2[2]),breaks = c(0,0.9,1))

dev.off()


############
### 2. Normalization by transforming to relative abundance
###########

###Ig
clone_type_IG_relative<-100*t(apply(clone_type_IG,1,function (x) x/sum(x))) ##Each cell divided by the total number of counts per sample
id<-match(rownames(clone_type_filter_IG),colnames(clone_type_IG_relative)) ##Filter by clones being in at least two samples
clone_type_IG_relative_filter<-clone_type_IG_relative[,id] ##1252

#### Logistoc regression
p_value=NULL
for(i in 1:dim(clone_type_IG_relative_filter)[2]){
  print(i)
  mod<-glm(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome~clone_type_IG_relative_filter[,i]+clone_type_filter_IG[i,],family = "binomial")
  p_value[i]=coefficients(summary(mod))[2,4]
}

p.adj<-p.adjust(p_value,"fdr")
clone_type_relative_sign<-clone_type_IG_relative_filter[,which(p.adj<0.05)] #215

#ENET
alphalist<-seq(0.1,0.9,by=0.01)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(clone_type_IG_relative_filter,PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome,family="binomial"
                                                        ,standardize=TRUE,alpha=a,nfolds=5))})
xx<-rep(NA,length(alphalist))
yy<-rep(NA,length(alphalist))
for (j in 1:length(alphalist)) {
  #print(j)
  if(class(elasticnet[[j]]) != "try-error"){
    xx[j]<-elasticnet[[j]]$lambda.min
    id.cv.opt<-grep(elasticnet[[j]]$lambda.min,elasticnet[[j]]$lambda,fixed=TRUE)
    yy[j]<-elasticnet[[j]]$cvm[id.cv.opt]
  }
}
id.min<-which(yy==min(yy,na.rm=TRUE))
lambda<-xx[id.min]
alpha<-alphalist[id.min]

enet<-glmnet(clone_type_IG_relative_filter,PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome,family="binomial",standardize=TRUE,alpha=alpha,lambda=lambda)
clones<-rownames(enet$beta)[which(enet$beta!=0)]

clone_type_relative_sign<-clone_type_IG_relative_filter[,match(clones,colnames(clone_type_IG_relative_filter))] #23

##Plot results
annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome)
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal)
ann_colors = list (outcome =  c("normal-pancreas (GTEx)" = cols[1],
                                "tumor-pancreas (TCGA)" = cols[2]))

tiff("Results/heatmap_relativeAbundance_Ig_sign.tiff",width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_relative_sign),border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(9,name="Reds"))(100))
dev.off()

clone_type_relative_sign_order<-clone_type_relative_sign[order(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome),]
df_long <- melt(clone_type_relative_sign_order, id.vars = "Sample", variable.name = "Clones")
colnames(df_long)<-c("Sample","Clones","value")
library(randomcoloR)
n <- 215
palette <- distinctColorPalette(n)
tiff("Results/barplot_relativeAbundance_Ig_sign.tiff",width = 5000, height = 3000, res = 300)
ggplot(df_long, aes(x = Sample, y = value, fill = Clones)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=palette)
dev.off()
