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

load("Data/PAAD_GTEx/PAAD_GTEx_FullData.Rdata")

##Restrict only to Normal vs Tumor 
PAAD.GTEx.repertoire.diversity.tumor.normmal<-PAAD.GTEx.repertoire.diversity[which(PAAD.GTEx.repertoire.diversity$outcome=="normal-pancreas (GTEx)"|
                                                                                     PAAD.GTEx.repertoire.diversity$outcome== "tumor-pancreas (TCGA)"),]
PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome<-factor(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome)

###################################################
## Common clones across samples
#################################################
id<-match(data_merge$sample,rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal))
data_merge<-data_merge[which(is.na(id)==F),]

#chain=c("IGHV","IGKV","IGLV")
############
## 1. Build the matrix with the clones by samples
data_chain_IGK<-data_merge[which(data_merge$chainType=="IGK"),]
clone_type_IGK<-t(as.data.frame(unclass(table(data_chain_IGK$V_J_lenghCDR3_CloneId,factor(data_chain_IGK$sample))))) 
data_chain_IGL<-data_merge[which(data_merge$chainType=="IGL"),]
clone_type_IGL<-t(as.data.frame(unclass(table(data_chain_IGL$V_J_lenghCDR3_CloneId,factor(data_chain_IGL$sample))))) 
data_chain_IGH<-data_merge[which(data_merge$chainType=="IGH"),]
clone_type_IGH<-t(as.data.frame(unclass(table(data_chain_IGH$V_J_lenghCDR3_CloneId,factor(data_chain_IGH$sample))))) 
data_chain_IG<-data_merge[which(data_merge$chainType=="IGH" |
                                        data_merge$chainType=="IGK" |
                                        data_merge$chainType=="IGL"),]
clone_type_IG<-t(as.data.frame(unclass(table(data_chain_IG$V_J_lenghCDR3_CloneId,factor(data_chain_IG$sample))))) 

data_chain_TCR<-data_merge[which(data_merge$chainType=="TRA" | 
                                         data_merge$chainType=="TRB" |
                                         data_merge$chainType=="TRD" | 
                                         data_merge$chainType=="TRG"),]
clone_type_TCR<-t(as.data.frame(unclass(table(data_chain_TCR$V_J_lenghCDR3_CloneId,factor(data_chain_TCR$sample))))) 

##Build present vs no present
clone_type_presence_IGK<-apply(clone_type_IGK,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGL<-apply(clone_type_IGL,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGH<-apply(clone_type_IGH,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IG<-apply(clone_type_IG,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TCR<-apply(clone_type_TCR,1,function(x) ifelse(x==0,0,1))

###Filter by clones share at least in 20% of the samples
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
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F","#BEAED4")

annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome)
ann_colors = list (outcome = c("normal-pancreas (GTEx)" = cols[1],"tumor-pancreas (TCGA)" = cols[2]))
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal)

cols2 = brewer.pal(10,name = "Paired")

#IG
tiff(paste0("Results/heatmap_common_clones_IG.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IG),border_color=F,show_rownames = F, show_colnames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols2[6]),breaks = c(0,0.9,1))
dev.off()

#IGH
tiff(paste0("Results/heatmap_common_clones_IGH.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGH),border_color=F,show_rownames = F, show_colnames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols2[2]),breaks = c(0,0.9,1))
dev.off()

#IGK
tiff(paste0("Results/heatmap_common_clones_IGK.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGK),show_rownames = F, show_colnames = F,border_color=F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols2[4]),breaks = c(0,0.9,1))
dev.off()

#IGL
tiff(paste0("Results/heatmap_common_clones_IGL.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IGL),border_color=F,show_rownames = F, show_colnames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols2[8]),breaks = c(0,0.9,1))
dev.off()

#TCR
tiff(paste0("Results/heatmap_common_clones_TCR.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_TCR),border_color=F,show_rownames = F, show_colnames = F,annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white",cols2[10]),breaks = c(0,0.9,1))
dev.off()

#####################################
### Principal Components Analysis ###
#####################################
library(factoextra)
pca <- prcomp(t(clone_type_filter_IG), scale = TRUE)

SPP <- PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome
levels.SPP <- factor(c("normal-pancreas (GTEx)", "tumor-pancreas (TCGA)"))

pc <- c(1,2)
tiff(paste0("Results/PCA_IG.tiff"),width = 2000, height = 2000, res = 300)
plot(pca$x[,pc[1]], pca$x[,pc[2]], col=cols[SPP],pch=20,xlab="PCA1",ylab="PCA2")
legend("bottomright", legend=levels(levels.SPP), col=cols,pch=20,cex=0.8)
dev.off()

#####Fisher test to find differences by clones in normal vs. tumor
p_value=NULL
for(i in 1:dim(clone_type_filter_IG)[1]){
  print(i)
  tab<-table(clone_type_filter_IG[i,],PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome)
  p_value[i]=fisher.test(tab)$p.value
}
p.adj<-p.adjust(p_value,"fdr")
clone_type_filter_IG_sign<-clone_type_filter_IG[which(p.adj<0.05),] #808 out of 1252 are significant
annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome)
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal)
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F", "#BEAED4")
ann_colors = list (outcome = c("normal-pancreas (GTEx)" = cols[1],
                                         "tumor-pancreas (TCGA)" = cols[2]))

tiff(paste0("Results/heatmap_common_clones_IG_sign.tiff"),width = 5000, height = 3000, res = 300)
pheatmap(t(clone_type_filter_IG_sign),border_color=F,show_rownames = F, annotation_row = annotation_row,
         annotation_colors = ann_colors,color = c("white","red"),breaks = c(0,0.9,1))
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

enet<-glmnet(xcell.data.tumor.filter_Igreads_mat,log10(PAAD.repertoire.diversity_Igreads$IG_expression),family="gaussian",standardize=TRUE,alpha=alpha,lambda=lambda)
cells<-rownames(enet$beta)[which(enet$beta!=0)]
coef<-enet$beta[which(enet$beta!=0)]

significant_cells<-xcell.data.tumor.filter_Igreads_mat[,match(cells,colnames(xcell.data.tumor.filter_Igreads_mat))] #15

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
