rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. Compositional Analysis
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP:
###         
###
### Author: Silvia Pineda
### Date: September, 2019
############################################################################################
library(zCompositions)
library(easyCODA)
library(ellipse)
library(factoextra)
library(RColorBrewer)
library(sva)
library(pheatmap)
library(vegan)
library(MASS)

##Colors: brewer.pal(4, "Accent")
##"#BEAED4" (IGH) "#7FC97F" (IGK) "#FDC086" (IGL)
setwd("~/TCGA-Immune/")

load("Data/PAAD/PAAD_FullData.Rdata")

##Filter for those that are pancreas or normal
PAAD.repertoire.diversity.tumor.normal<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$outcome=="normal-pancreas" |
                                        PAAD.repertoire.diversity$outcome=="tumor-pancreas (TCGA)"),]
data_merge_qc<-data_merge[which(is.na(match(data_merge$sample,rownames(PAAD.GTEx.repertoire.diversity.tumor.normal)))==F),]

############
## 1. Build the matrix with the clones by samples
###########
## Pancreas ans pseudonormal
data_qc_chain_IGK<-data_merge_qc[which(data_merge_qc$chainType=="IGK"),]
clone_type_IGK<-t(as.data.frame(unclass(table(data_qc_chain_IGK$V_J_lenghCDR3_CloneId,factor(data_qc_chain_IGK$sample))))) 
data_qc_chain_IGL<-data_merge_qc[which(data_merge_qc$chainType=="IGL"),]
clone_type_IGL<-t(as.data.frame(unclass(table(data_qc_chain_IGL$V_J_lenghCDR3_CloneId,factor(data_qc_chain_IGL$sample))))) 
data_qc_chain_IGH<-data_merge_qc[which(data_merge_qc$chainType=="IGH"),]
clone_type_IGH<-t(as.data.frame(unclass(table(data_qc_chain_IGH$V_J_lenghCDR3_CloneId,factor(data_qc_chain_IGH$sample))))) 
data_qc_chain_IG<-data_merge_qc[which(data_merge_qc$chainType=="IGH" |
                                        data_merge_qc$chainType=="IGK" |
                                        data_merge_qc$chainType=="IGL"),]
clone_type_IG<-t(as.data.frame(unclass(table(data_qc_chain_IG$V_J_lenghCDR3_CloneId,factor(data_qc_chain_IG$sample))))) 

data_qc_chain_TCR<-data_merge_qc[which(data_merge_qc$chainType=="TRA" | 
                                         data_merge_qc$chainType=="TRB" |
                                         data_merge_qc$chainType=="TRD" | 
                                         data_merge_qc$chainType=="TRG"),]
clone_type_TCR<-t(as.data.frame(unclass(table(data_qc_chain_TCR$V_J_lenghCDR3_CloneId,factor(data_qc_chain_TCR$sample))))) 

#We filter out non-expressed genes, by requiring more than 5 reads in at least two samples for each gene.
filter <- apply(clone_type_IG,2, function(x) length(x[x>5])>=2)
clone_type_IG_filter<-clone_type_IG[,filter] #635
### Calculate relative abundances 
clone_type_IG_relative <- clone_type_IG_filter / colSums(clone_type_IG_filter)

#logtransform (need to add a constant for the zeros which is half the min relative abudance)
min<-min(clone_type_IG_relative[clone_type_IG_relative!=0])/2
clone_type_IG_relative <- clone_type_IG_relative+min
clone_type_IG_relative_log<-log2(clone_type_IG_relative)

# ##Build present vs no present
# clone_type_presence_IGK<-apply(clone_type_IGK,1,function(x) ifelse(x==0,0,1))
# clone_type_presence_IGL<-apply(clone_type_IGL,1,function(x) ifelse(x==0,0,1))
# clone_type_presence_IGH<-apply(clone_type_IGH,1,function(x) ifelse(x==0,0,1))
 clone_type_presence_IG<-apply(clone_type_IG_filter,1,function(x) ifelse(x==0,0,1))
# clone_type_presence_TCR<-apply(clone_type_TCR,1,function(x) ifelse(x==0,0,1))
# 
# #Filter by at least having two common samples
# clone_type_filter_IG<-clone_type_presence_IG[which(rowSums(clone_type_presence_IG)>2),] #
# clone_type_filter_IGK<-clone_type_presence_IGK[which(rowSums(clone_type_presence_IGK)>2),] #
# clone_type_filter_IGL<-clone_type_presence_IGL[which(rowSums(clone_type_presence_IGL)>2),] #
# clone_type_filter_IGH<-clone_type_presence_IGH[which(rowSums(clone_type_presence_IGH)>2),] #
# clone_type_filter_IG<-clone_type_presence_IG[which(rowSums(clone_type_presence_IG)>2),] #
# clone_type_filter_TCR<-clone_type_presence_TCR[which(rowSums(clone_type_presence_TCR)>2),] #


### Correct by batch effect using Combat
clone_type_IG_relative_log_combat = ComBat(dat=t(clone_type_IG_relative_log), 
                                                  batch=PAAD.GTEx.repertoire.diversity.tumor.normal$outcome)


####Write the matrix to run the percentile normalization
##Filter by at least having two samples the clones

write.table(clone_type_IG_filter,"Data/PAAD_GTEx/clone_type_IG.txt",quote=F,sep="\t")
write.table(PAAD.GTEx.repertoire.diversity$sample[which(PAAD.GTEx.repertoire.diversity$outcome=="normal-pancreas (GTEx)")],
            "Data/PAAD_GTEx/control_samples.txt",row.names = F,
            quote = F)
write.table(PAAD.GTEx.repertoire.diversity$sample[which(PAAD.GTEx.repertoire.diversity$outcome=="tumor-pancreas (TCGA)")],
            "Data/PAAD_GTEx/cases_samples.txt",row.names = F, 
            quote = F)

##Read after applying percentile normalization from python script
clone_type_IG_relative_abundances_perc_norm<-read.table("Data/PAAD_GTEx/out_percentile_norm.txt")
id<-match(rownames(PAAD.GTEx.repertoire.diversity.tumor.normal),rownames(clone_type_IG_relative_abundances_perc_norm))
clone_type_IG_relative_abundances_perc_norm<-clone_type_IG_relative_abundances_perc_norm[id,]

####SVA
count_matrix<-clone_type_IG_filter
coldata<-matrix(NA,nrow(count_matrix),2)
coldata[,2]<-PAAD.GTEx.repertoire.diversity.tumor.normal$outcome
coldata[,1]<-ifelse(coldata[,2]=="normal-pancreas (GTEx)","GTEx","TCGA")
rownames(coldata)<-rownames(PAAD.GTEx.repertoire.diversity.tumor.normal)
colnames(coldata)<-c("batch","outcome")

dds <- DESeqDataSetFromMatrix(t(count_matrix), coldata, ~ batch) 
normalized_vst <- varianceStabilizingTransformation(dds)
norm_data_vst<-assay(normalized_vst) 

################################################################################
### Principal Components Analysis and boxplot to check for the normalization###
##############################################################################
SPP <- factor(PAAD.GTEx.repertoire.diversity.tumor.normal$outcome)
levels.SPP <- factor(c("normal-pancreas (GTEx)", "tumor-pancreas (TCGA)"))

brewer.pal(4,name = "Accent")
cols=c( "#7FC97F","#BEAED4")

pc <- c(1,2)
##Present clones
tiff(paste0("Results/PCA_clones_present_IG.tiff"),width = 2000, height = 2000, res = 300)
pca <- prcomp(t(clone_type_presence_IG)) #nxp
plot(pca$x[,pc[1]], pca$x[,pc[2]], col=cols[SPP],pch=20,xlab="PCA1",ylab="PCA2",main="clones defined as present vs. not present")
legend("bottomright", legend=levels(levels.SPP), col=cols,pch=20,cex=0.8)
dev.off()

##Relative abundances
tiff(paste0("Results/PCA_clones_relative_abundance_IG.tiff"),width = 2000, height = 2000, res = 300)
pca <- prcomp(clone_type_IG_relative_log) #nxp
plot(pca$x[,pc[1]], pca$x[,pc[2]], col=cols[SPP],pch=20,xlab="PCA1",ylab="PCA2",main="clones defined as relative abundance")
legend("bottomright", legend=levels(levels.SPP), col=cols,pch=20,cex=0.8)
dev.off()

##Relative abundance corrected by combat
tiff(paste0("Results/PCA_clones_relative_abundance_IG_combat.tiff"),width = 2000, height = 2000, res = 300)
pca <- prcomp(t(clone_type_IG_relative_log_combat))
plot(pca$x[,pc[1]], pca$x[,pc[2]], col=cols[SPP],pch=20,xlab="PCA1",ylab="PCA2",main="clones defined as relative abundance corrected by ComBat")
legend("bottomright", legend=levels(levels.SPP), col=cols,pch=20,cex=0.8)
dev.off()

##Relative abundance corrected by percentile
tiff(paste0("Results/PCA_clones_relative_abundance_IG_perc_norm.tiff"),width = 2000, height = 2000, res = 300)
pca <- prcomp(clone_type_IG_relative_abundances_perc_norm)
plot(pca$x[,pc[1]], pca$x[,pc[2]], col=cols[SPP],pch=20,xlab="PCA1",ylab="PCA2",main="clones defined as relative abundance corrected by percNorm")
legend("bottomright", legend=levels(levels.SPP), col=cols,pch=20,cex=0.8)
dev.off()

################
##Heatmaps#####
##############
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F","#BEAED4")

annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normal$outcome)
ann_colors = list (outcome = c("normal-pancreas (GTEx)" = cols[1],"tumor-pancreas (TCGA)" = cols[2]))
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(PAAD.GTEx.repertoire.diversity.tumor.normal)

cols2 = brewer.pal(10,name = "Paired")

#IG
tiff(paste0("Results/heatmap_common_clones_IG.tiff"),width = 5000, height = 3000, res = 300)
  pheatmap(t(clone_type_IG_relative_log_combat),scale="row",border_color=F,show_rownames = F, show_colnames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(9,name="PuOr"))(10))
dev.off()
