rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. Clustering Analysis 
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Clustering
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
library(ggpubr)
library(viridis)
setwd("~/TCGA-Immune/")

load("Data/PAAD/PAAD_FullData.Rdata")

####Filter by p-value
xCell.pvalue.filter<-t(apply(xCell.pvalue.PAAD,1,function (x) replace(x,x>=0.2,1)))
xcell.data.filter <-  t(xCell.data.PAAD[rowSums(xCell.pvalue.filter==1) <= 160*0.8,]) ##(160*0.8)  45 cells


####Filter by variabilty (mean +- 2*SD)
filter.out<-NULL
for(i in 1:ncol(xcell.data.filter)){
  pos<-mean(xcell.data.filter[,i])+2*sd(xcell.data.filter[,i])
  neg<-mean(xcell.data.filter[,i])-2*sd(xcell.data.filter[,i])
  
  if(length(which(xcell.data.filter[,i]> pos)) == 0 & 
     length(which(xcell.data.filter[,i]< neg)) == 0){
    
    filter.out[i]<-i
  } else {
    filter.out[i]<-NA
  }
}

#Heatmp
annotation_col = data.frame(PAAD.repertoire.diversity$Tumor_type_4categ)
cols=c( "#7FC97F", "#FBB4AE","#BEAED4", "#FDC086")
ann_colors = list (Tumor_type_4categ = c("normal_pancreas" = cols[1],
                                         "PAC-Other" = cols[2],
                                         "PDAC"= cols[3],
                                         "pseudonormal_pancreas" = cols[4]))
colnames(annotation_col)<-"Tumor_type_4categ"
rownames(annotation_col)<-PAAD.repertoire.diversity$TCGA_sample

cols = colorRampPalette(brewer.pal(6,name="PuOr"))(12)
tiff("Results/xCell_heatmap.tiff",res=300,w=3500,h=3000)
pheatmap(t(xcell.data.filter),scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12),
         annotation_col = annotation_col,annotation_colors = ann_colors,cluster_cols=F)
dev.off()


# #Delete the scores
# mat<-t(xcell.data.filter[ , -which(colnames(xcell.data.filter) %in% c("ImmuneScore","StromaScore","MicroenvironmentScore"))]) ##Delete the scores
# tiff("Results/xCell_CompareSamples.tiff",res=300,w=3500,h=3000)
# pheatmap(mat,scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))
# dev.off()


###################################
##### Only tumors #################
###################################
PAAD.repertoire.tumor<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type_4categ=="PDAC"),]
PAAD.repertoire.tumor$TCGA_sample<-substr(PAAD.repertoire.tumor$TCGA_sample,1,15)
xCell.data.tumor<-t(xCell.data.PAAD[,match(PAAD.repertoire.tumor$TCGA_sample,colnames(xCell.data.PAAD))])
xCell.pvalue.tumor<-t(xCell.pvalue.PAAD[,match(PAAD.repertoire.tumor$TCGA_sample,colnames(xCell.pvalue.PAAD))])

####Filter by p-value
xCell.pvalue.tumor.filter<-t(apply(xCell.pvalue.tumor,1,function (x) replace(x,x>=0.2,1)))
xcell.data.tumor.filter <-  xCell.data.tumor[,colSums(xCell.pvalue.tumor.filter==1) <= 131*0.8] ##(131*0.8)  45 cells

####Filter by variabilty (mean +- 2*SD)
filter.out<-NULL
for(i in 1:ncol(xcell.data.tumor.filter)){
  pos<-mean(xcell.data.tumor.filter[,i])+2*sd(xcell.data.tumor.filter[,i])
  neg<-mean(xcell.data.tumor.filter[,i])-2*sd(xcell.data.tumor.filter[,i])
  
  if(length(which(xcell.data.tumor.filter[,i]> pos)) == 0 & 
     length(which(xcell.data.tumor.filter[,i]< neg)) == 0){
    
    filter.out[i]<-i
  } else {
    filter.out[i]<-NA
  }
}

#Delete the scores
#mat<-t(xcell.data.tumor.filter[ , -which(colnames(xcell.data.tumor.filter) %in% c("ImmuneScore","StromaScore","MicroenvironmentScore"))]) ##Delete the scores
tiff("Results/xCell_Tumor.tiff",res=300,w=3500,h=3000)
pheatmap(t(xcell.data.tumor.filter),scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))
dev.off()

####################
#### Clustering ###
###################
mat<-t(xCell.data.tumor)
# Prepare Data
mydata <- scale(t(mat)) # standardize variables
rownames(mydata)<-substr(rownames(mydata),1,12)
# Determine number of clusters
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

###############
# 1.  K-Means Cluster Analysis
################
fit.k <- kmeans(mydata, 3) # 2 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit.k$cluster),FUN=mean)
library(fpc)
tiff("Results/xCell_cluster_kmeans.tiff",res=300,w=2000,h=2000)
plotcluster(mydata, fit.k$cluster)
dev.off()

library(cluster) 
tiff("Results/xCell_cluster_kmeans2.tiff",res=300,w=2000,h=2000)
clusplot(mydata, fit.k$cluster, color=TRUE, shade=TRUE, lines=0)
dev.off()

library(factoextra)
tiff("Results/xCell_fviz_cluster_kmeans3.tiff",res=300,w=2000,h=2000)
fviz_cluster(object = fit.k, data = mydata, show.clust.cent = TRUE,labelsize = 8,
             star.plot = TRUE, repel = TRUE) +
  labs(title = "Results clustering K-means") +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

# library(ComplexHeatmap)
# library(viridis)
# library("circlize")
# colors <- colorRampPalette(brewer.pal(6,name="PuOr"))(5)
# cols = colorRamp2(c(-6,-3,0,3,6),c("#B35806" ,"#F4B25D", "#EBDDD0", "#A8A1CD" ,"#542788"))
# Heatmap(matrix = mat, name = "K-means", col = cols,
#         column_title = "samples",
#         row_names_gp = gpar(fontsize = 7),
#         clustering_distance_columns = "euclidean",
#         clustering_distance_rows = "euclidean",
#         clustering_method_columns = "average",
#         clustering_method_rows = "average",
#         km = 3)
# 

cluster.order<-fit.k$cluster[order(fit.k$cluster)]
annotation_col = data.frame(cluster = factor(cluster.order))
colnames(mat)<-substr(colnames(mat),1,12)
mat.k<-mat[,match(rownames(annotation_col),colnames(mat))]
ann_colors = list (cluster = c("1" = brewer.pal(3,"Set1")[1], "2"= brewer.pal(3,"Set1")[3], "3" = brewer.pal(3,"Set1")[2]))
tiff("Results/xCell_heatmap_kmeans3.tiff",res=300,w=3500,h=3000)
pheatmap(mat.k,scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12),
         annotation_col = annotation_col,annotation_colors = ann_colors,cluster_cols=F)
dev.off()
xcell.cluster<-fit.k$cluster
id<-match(substr(rownames(xcell.data.tumor.filter),1,12),names(xcell.cluster))
xcell.data.tumor.filter<-data.frame(xcell.data.tumor.filter)
xcell.data.tumor.filter$cluster<-xcell.cluster[id]

#################
# 2. Ward Hierarchical Clustering
#################
d <- dist(mydata, method = "euclidean") # distance matrix
fit.h <- hclust(d, method="ward.D2") 
plot(fit.h) # display dendogram
groups <- cutree(fit.h, k=4) # cut tree into 2 clusters
# draw dendogram with red borders around the 2 clusters 
rect.hclust(fit.h, k=2, border="red")

tiff("Results/xCell_cluster_hierarchical.tiff",res=300,w=2000,h=2000)
plotcluster(mydata, groups)
dev.off()

tiff("Results/xCell_cluster_hierarchical2.tiff",res=300,w=2000,h=2000)
clusplot(mydata, groups, color=TRUE, shade=TRUE, lines=0)
dev.off()

annotation_col = data.frame(cluster = factor(groups))
ann_colors = list (cluster = c("1" = brewer.pal(3,"Set2")[1], "2"= brewer.pal(3,"Set2")[2], "3" = brewer.pal(3,"Set2")[3],"4" = brewer.pal(4,"Set2")[4]))
pheatmap(mat,scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12),
         annotation_col = annotation_col,annotation_colors = ann_colors,clustering_method = "ward.D2")

#############
# 3. Model Based Clustering
############
library(mclust)
fit <- Mclust(mydata)
plot(fit) # plot results 
summary(fit) # display the best model

##Plotting kmeans with 2 clusters
tiff("Results/xCell_cluster_mclust.tiff",res=300,w=2000,h=2000)
plotcluster(mydata, fit$classification)
dev.off()

tiff("Results/xCell_cluster_mclust2.tiff",res=300,w=2000,h=2000)
clusplot(mydata, fit$classification, color=TRUE, shade=TRUE, lines=0)
dev.off()



##############################################
### Survival Analysis  ####
##############################################
library(survival)
library(survminer)
library(survMisc)

####Xcell
clinical.follow.up.xcell<-clinical.folow_up[match(substr(rownames(xcell.data.tumor.filter),1,12),clinical.folow_up$bcr_patient_barcode),]
clinical.patient.follow.up.xcell<-merge(clinical.follow.up.xcell,clinical.patient,by="bcr_patient_barcode")
#Pathologic_stage
clinical.patient.follow.up.xcell$pathologic_stage<-factor(ifelse(clinical.patient.follow.up.xcell$stage_event_pathologic_stage=="Stage IA" | 
                                                                   clinical.patient.follow.up.xcell$stage_event_pathologic_stage=="Stage IB", "Stage I",
                                                       ifelse(clinical.patient.follow.up.xcell$stage_event_pathologic_stage == "Stage IIA" |
                                                                clinical.patient.follow.up.xcell$stage_event_pathologic_stage=="Stage IIB","Stage II",
                                                              ifelse(clinical.patient.follow.up.xcell$stage_event_pathologic_stage=="Stage III","Stage III",
                                                                     ifelse(clinical.patient.follow.up.xcell$stage_event_pathologic_stage=="Stage IV", "Stage IV",NA)))))

##OS
surv_object <- Surv(time = clinical.patient.follow.up.xcell$OS.time, event = clinical.patient.follow.up.xcell$OS)
xcell.data.tumor.filter$cluster_2cat<-ifelse(xcell.data.tumor.filter$cluster=="1" | xcell.data.tumor.filter$cluster=="2", "Cluster1+2","Cluster3")
res.cox <- coxph(surv_object~xcell.data.tumor.filter[,"cluster_2cat"]+clinical.patient.follow.up.xcell$pathologic_stage)
summary(res.cox)
##Categorical
fit1 <- survfit(surv_object ~  factor(xcell.data.tumor.filter[,"cluster_2cat"]))
fit1
tiff("Results/xCell_cluster_OS.tiff",res=300,w=2000,h=2000)
ggsurvplot(fit1, data = as.data.frame(xcell.data.tumor.filter),pval = "p = 0.03 \n(Gehan-Breslow)" ,conf.int = T,risk.table = TRUE,
           legend.labs=c("Cluster 1+2","Cluster 3"))
dev.off()
comp(ten(fit1))$tests$lrTests

#DSS
clinical.follow.up.xcell.DSS<-clinical.follow.up.xcell[which(clinical.follow.up.xcell$DSS!="#N/A"),]
clinical.follow.up.xcell.DSS$DSS<-as.integer(as.character(clinical.follow.up.xcell.DSS$DSS))
surv_object <- Surv(time = clinical.follow.up.xcell.DSS$DSS.time, event = clinical.follow.up.xcell.DSS$DSS)
res.cox <- coxph(surv_object~xcell.data.tumor.filter[,"cluster_2cat"][which(clinical.follow.up.xcell.DSS$DSS!="#N/A")])
summary(res.cox)
#Categorical
fit1 <- survfit(surv_object ~ xcell.data.tumor.filter[,"cluster_2cat"][which(clinical.follow.up.xcell.DSS$DSS!="#N/A")])
fit1
ggsurvplot(fit1, data = xcell.data.tumor.filter[which(clinical.follow.up.xcell.DSS$DSS!="#N/A"),])
comp(ten(fit1))$tests$lrTests


##PFI
surv_object <- Surv(time = clinical.follow.up.xcell$PFI.time, event = clinical.follow.up.xcell$PFI)
res.cox <- coxph(surv_object~xcell.data.tumor.filter[,"cluster_2cat"])
summary(res.cox)
fit1 <- survfit(surv_object ~  factor(xcell.data.tumor.filter[,"cluster_2cat"]))
fit1
comp(ten(fit1))$tests$lrTests


