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
PAAD.repertoire.diversity$Tumor_type_2categ<-ifelse(PAAD.repertoire.diversity$Tumor_type=="Tumor_pancreas","Tumor_pancres",
                                                    ifelse(PAAD.repertoire.diversity$Tumor_type=="Solid_tissue_normal","Adjacent_normal_pancreas",
                                                           ifelse(PAAD.repertoire.diversity$Tumor_type=="Adjacent_normal_pancreas","Adjacent_normal_pancreas",NA)))
PAAD.repertoire.diversity$Tumor_type_2categ<-as.factor(PAAD.repertoire.diversity$Tumor_type_2categ)

###################################
##### Only tumors #################
###################################
PAAD.repertoire.tumor<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type=="Tumor_pancreas"),] #14
PAAD.repertoire.tumor$TCGA_sample<-substr(PAAD.repertoire.tumor$TCGA_sample,1,15)
xCell.data.tumor<-t(xCell.data.PAAD[,match(PAAD.repertoire.tumor$TCGA_sample,colnames(xCell.data.PAAD))])
xCell.pvalue.tumor<-t(xCell.pvalue.PAAD[,match(PAAD.repertoire.tumor$TCGA_sample,colnames(xCell.pvalue.PAAD))])

####Filter by p-value
xCell.pvalue.tumor.filter<-t(apply(xCell.pvalue.tumor,1,function (x) replace(x,x>=0.2,1)))
xcell.data.tumor.filter <-  xCell.data.tumor[,colSums(xCell.pvalue.tumor.filter==1) <= 148*0.8] ##(148*0.8)  45 cells

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
mat<-t(xcell.data.tumor.filter[ , -which(colnames(xcell.data.tumor.filter) %in% c("ImmuneScore","StromaScore","MicroenvironmentScore"))]) ##Delete the scores
tiff("Results/xCell_CompareSamples.tiff",res=300,w=3500,h=3000)
pheatmap(mat,scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))
dev.off()

####################
#### Clustering ###
###################

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

library(ComplexHeatmap)
library(viridis)
colors <- colorRampPalette(brewer.pal(6,name="PuOr"))(12)
Heatmap(matrix = mat, name = "K-means", col = colors,
        column_title = "samples",
        row_names_gp = gpar(fontsize = 7),
        clustering_distance_columns = "euclidean",
        clustering_distance_rows = "euclidean",
        clustering_method_columns = "average",
        clustering_method_rows = "average",
        km = 3)


annotation_col = data.frame(cluster = factor(fit.k$cluster))
ann_colors = list (cluster = c("1" = brewer.pal(3,"Set2")[1], "2"= brewer.pal(3,"Set2")[2], "3" = brewer.pal(3,"Set2")[3],"4" = brewer.pal(4,"Set2")[4]))
pheatmap(mat,scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12),
         annotation_col = annotation_col,annotation_colors = ann_colors,clustering_method = "ward.D2")

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

generate_breaks = function(x, n, center = F){
  if(center){
    m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    res = seq(-m, m, length.out = n + 1)
  }
  else{
    res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
  }
  
  return(res)
}

scaled_data <- scale_rows(mat)

color =  colorRampPalette(brewer.pal(6,name="PuOr"))(12)
breaks = generate_breaks(scaled_data, length(color), center = T)
pheatmap(scaled_data, breaks = breaks, color = color)

xcell.cluster<-fit.k$cluster


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
