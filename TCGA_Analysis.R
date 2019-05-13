rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. Descriptive analysis 
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
library(ggpubr)
library(viridis)
setwd("~/TCGA-Immune/")

load("Data/PAAD/PAAD_FullData.Rdata")
PAAD.repertoire.diversity$Tumor_type_2categ<-ifelse(PAAD.repertoire.diversity$Tumor_type=="Tumor_pancreas","Tumor_pancres",
                                             ifelse(PAAD.repertoire.diversity$Tumor_type=="Solid_tissue_normal","Adjacent_normal_pancreas",
                                             ifelse(PAAD.repertoire.diversity$Tumor_type=="Adjacent_normal_pancreas","Adjacent_normal_pancreas",NA)))
PAAD.repertoire.diversity$Tumor_type_2categ<-as.factor(PAAD.repertoire.diversity$Tumor_type_2categ)
##################
####Descriptive analysis to see if there are differences by tumor and adjacent_normal
##################
###T markers
PAAD.repertoire.diversity_treads<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$T_Reads>100),]
T_markers<-c("TRA_expression","TRB_expression","TRD_expression","TRG_expression","Alpha_Beta_ratio_expression",
             "clones_recon_TRA","clones_recon_TRB","clones_recon_TRD","clones_recon_TRG","entropy_recon_TRA",
             "entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG")
p.T_markers=NULL
for(i in 1:length(T_markers)){
  p.T_markers[i]<-coef(summary(glm(PAAD.repertoire.diversity_treads[,T_markers[i]]~PAAD.repertoire.diversity_treads$Tumor_type_2categ)))[2,4]
  PAAD.repertoire.diversity_treads$Marker<-PAAD.repertoire.diversity_treads[,T_markers[i]]
  tiff(paste0("Results/boxplot_",T_markers[i],".tiff"),res=300,h=2000,w=2000)
  print(ggplot(PAAD.repertoire.diversity_treads[which(is.na(PAAD.repertoire.diversity_treads$Tumor_type_2categ)==F),], 
         aes(y=Marker, fill=Tumor_type_2categ, x=Tumor_type_2categ)) + geom_boxplot() + scale_y_continuous(name=T_markers[i])
        + stat_compare_means())
  dev.off()
  
}
T_markers[which(p.T_markers<0.05)] ## "TRD_expression" "TRG_expression"

##Ig markers
PAAD.repertoire.diversity_Igreads<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$IG_Reads>100),]
Ig_markers<-c("IGH_expression","IGK_expression","IGL_expression","KappaLambda_ratio_expression",
             "clones_recon_IGH","clones_recon_IGK","clones_recon_IGL","entropy_recon_IGH",
             "entropy_recon_IGK","entropy_recon_IGL")
p.Ig_markers=NULL
for(i in 1:length(Ig_markers)){
  p.Ig_markers[i]<-coef(summary(glm(PAAD.repertoire.diversity_treads[,Ig_markers[i]]~PAAD.repertoire.diversity_treads$Tumor_type_2categ)))[2,4]
  PAAD.repertoire.diversity_treads$Marker<-PAAD.repertoire.diversity_treads[,Ig_markers[i]]
  tiff(paste0("Results/boxplot_",Ig_markers[i],".tiff"),res=300,h=2000,w=2000)
  print(ggplot(PAAD.repertoire.diversity_treads[which(is.na(PAAD.repertoire.diversity_treads$Tumor_type_2categ)==F),], 
               aes(y=Marker, fill=Tumor_type_2categ, x=Tumor_type_2categ)) + geom_boxplot() + scale_y_continuous(name=Ig_markers[i])
        + stat_compare_means())
  dev.off()
}
Ig_markers[which(p.Ig_markers<0.05)] ## 0

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
# Determine number of clusters
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
                                     centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

###############
# 1.  K-Means Cluster Analysis
################
fit.k <- kmeans(mydata, 2) # 2 cluster solution
# get cluster means 
aggregate(mydata,by=list(fit.k$cluster),FUN=mean)
# append cluster assignment
#mydata <- data.frame(mydata, fit.k$cluster)
library(fpc)
tiff("Results/xCell_cluster_kmeans.tiff",res=300,w=2000,h=2000)
plotcluster(mydata, fit.k$cluster)
dev.off()

library(cluster) 
tiff("Results/xCell_cluster_kmeans2.tiff",res=300,w=2000,h=2000)
clusplot(mydata, fit.k$cluster, color=TRUE, shade=TRUE, lines=0)
dev.off()

#################
# 2. Ward Hierarchical Clustering
#################
d <- dist(mydata, method = "euclidean") # distance matrix
fit.h <- hclust(d, method="ward.D2") 
plot(fit.h) # display dendogram
groups <- cutree(fit.h, k=2) # cut tree into 2 clusters
# draw dendogram with red borders around the 2 clusters 
rect.hclust(fit.h, k=2, border="red")

tiff("Results/xCell_cluster_hierarchical.tiff",res=300,w=2000,h=2000)
plotcluster(mydata, groups)
dev.off()

tiff("Results/xCell_cluster_hierarchical2.tiff",res=300,w=2000,h=2000)
clusplot(mydata, groups, color=TRUE, shade=TRUE, lines=0)
dev.off()

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

annotation_col = data.frame(cluster = groups)
ann_colors = list (cluster = c("1" = brewer.pal(3,"Set2")[1], "2"= brewer.pal(3,"Set2")[2]))
pheatmap(mat,scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12),
         annotation_col = annotation_col,annotation_colors = ann_colors,clustering_method = "ward.D2")

xcell.cluster<-fit.k$cluster
  
###Applied ENET to find cells associated with T or B expression
##############
##T_reads
##############
load("~/Downloads/PAAD_repertoire_xCell_clinical.Rdata")
PAAD.repertoire.diversity_treads<-PAAD.repertoire.tumor[which(PAAD.repertoire.tumor$T_Reads>100),] #139
xcell.data.tumor.filter_treads<-xcell.data.tumor.filter[which(PAAD.repertoire.tumor$T_Reads>100),] #139
xcell.data.tumor.filter_treads_mat<-xcell.data.tumor.filter_treads[ , -which(colnames(xcell.data.tumor.filter_treads) %in% c("ImmuneScore","StromaScore","MicroenvironmentScore"))]

library(dplyr)
PAAD.repertoire.diversity_treads$T_expression_tertiles<-ntile(PAAD.repertoire.diversity_treads$T_expression,3)

alphalist<-seq(0.01,0.99,by=0.01)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(xcell.data.tumor.filter_treads_mat,PAAD.repertoire.diversity_treads$TRA_expression,family="gaussian"
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

enet<-glmnet(xcell.data.tumor.filter_treads_mat,PAAD.repertoire.diversity_treads$TRA_expression,family="gaussian",standardize=TRUE,alpha=alpha,lambda=lambda)
cells<-rownames(enet$beta)[which(enet$beta!=0)]
coef<-enet$beta[which(enet$beta!=0)]

significant_cells<-xcell.data.tumor.filter_treads_mat[,match(cells,colnames(xcell.data.tumor.filter_treads_mat))] #12

## Heatmap ####
annotation_col = data.frame(
  T_expression = PAAD.repertoire.diversity_treads$TRA_expression)

rownames(annotation_col)<-rownames(significant_cells)
ann_colors = list (T_expression = brewer.pal(9,"Reds"))
tiff("Results/T_expression.tiff",res=300,w=3500,h=3000)
pheatmap(t(significant_cells),scale="row",annotation_col = annotation_col,annotation_colors = ann_colors,show_colnames = F,border_color=F,
         color = colorRampPalette(brewer.pal(9,name="PuOr"))(12))
dev.off()

##############
##IG_reads
#############
PAAD.repertoire.diversity_Igreads<-PAAD.repertoire.tumor[which(PAAD.repertoire.tumor$IGH_Reads>100),] #147
xcell.data.tumor.filter_Igreads<-xcell.data.tumor.filter[which(PAAD.repertoire.tumor$IGH_Reads>100),] #147
xcell.data.tumor.filter_Igreads_mat<-xcell.data.tumor.filter_Igreads[ , -which(colnames(xcell.data.tumor.filter_Igreads) %in% c("ImmuneScore","StromaScore","MicroenvironmentScore"))]

alphalist<-seq(0.01,0.99,by=0.01)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(xcell.data.tumor.filter_Igreads_mat,PAAD.repertoire.diversity_Igreads$IG_expression,family="gaussian"
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

enet<-glmnet(xcell.data.tumor.filter_Igreads_mat,PAAD.repertoire.diversity_Igreads$IG_expression,family="gaussian",standardize=TRUE,alpha=alpha,lambda=lambda)
cells<-rownames(enet$beta)[which(enet$beta!=0)]
coef<-enet$beta[which(enet$beta!=0)]

significant_cells<-xcell.data.tumor.filter_Igreads_mat[,match(cells,colnames(xcell.data.tumor.filter_Igreads_mat))] #15

## Heatmap ####
annotation_col = data.frame(
  IG_expression = PAAD.repertoire.diversity_Igreads$IG_expression)

rownames(annotation_col)<-rownames(significant_cells)
ann_colors = list (IG_expression = brewer.pal(6,"Greens"))
tiff("Results/IG_expression.tiff",res=300,w=3500,h=3000)
pheatmap(t(significant_cells),scale="row",annotation_col = annotation_col,annotation_colors = ann_colors,show_colnames = F,border_color=F,
         color = colorRampPalette(brewer.pal(9,name="PuOr"))(12))
dev.off()



###############################
## Merge with Clinical data ###
###############################
clinical.patient.tumor<-clinical.patient[match(substr(PAAD.repertoire.tumor$TCGA_sample,1,12),clinical.patient$bcr_patient_barcode),]
biospecimen.slide.tumor<-biospecimen.slide[match(substr(PAAD.repertoire.tumor$TCGA_sample,1,12),biospecimen.slide$bcr_patient_barcode),]

markers<-c("IG_expression","IGH_expression","IGK_expression","IGL_expression","T_expression","TRA_expression","TRB_expression","TRD_expression",
  "TRG_expression","Alpha_Beta_ratio_expression", "KappaLambda_ratio_expression","clones_recon_IGH","clones_recon_IGK","clones_recon_IGL",
  "clones_recon_TRA","clones_recon_TRB","clones_recon_TRD","clones_recon_TRG","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL",           
  "entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG")

##Histological type
clinical.patient.tumor$histological_type_2cat<-factor(ifelse(clinical.patient.tumor$histological_type=="Pancreas-Adenocarcinoma Ductal Type","PDAC","Other"))
clinical.patient.tumor$histological_type_2cat<-factor(clinical.patient.tumor$histological_type_2cat)
association.test.immuneRep(clinical.patient.tumor,"histological_type_2cat",PAAD.repertoire.tumor,markers)

##anatomic_neoplasm_subdivision
clinical.patient.tumor$anatomic_neoplasm_subdivision<-factor(clinical.patient.tumor$anatomic_neoplasm_subdivision)
association.test.immuneRep(clinical.patient.tumor,"anatomic_neoplasm_subdivision",PAAD.repertoire.tumor,markers)

##gender
clinical.patient.tumor$gender<-factor(clinical.patient.tumor$gender)
association.test.immuneRep(clinical.patient.tumor,"gender",PAAD.repertoire.tumor,markers)

##race_list
clinical.patient.tumor$race_list<-factor(clinical.patient.tumor$race_list)
association.test.immuneRep(clinical.patient.tumor,"race_list",PAAD.repertoire.tumor,markers)

##History of Prior Malignancy
clinical.patient.tumor$other_dx<-factor(clinical.patient.tumor$other_dx)
association.test.immuneRep(clinical.patient.tumor,"other_dx",PAAD.repertoire.tumor,markers)

##lymph_node_examined_count
association.test.immuneRep(clinical.patient.tumor,"lymph_node_examined_count",PAAD.repertoire.tumor,markers)

##neoplasm_histologic_grade
clinical.patient.tumor$neoplasm_histologic_grade_3cat<-factor(ifelse(clinical.patient.tumor$neoplasm_histologic_grade=="G1","G1",
                                                                     ifelse(clinical.patient.tumor$neoplasm_histologic_grade=="G2","G2",
                                                                            ifelse(clinical.patient.tumor$neoplasm_histologic_grade=="G3","G3",""))))
association.test.immuneRep(clinical.patient.tumor,"neoplasm_histologic_grade_3cat",PAAD.repertoire.tumor,markers)

##Age 
association.test.immuneRep(clinical.patient.tumor,"age_at_initial_pathologic_diagnosis",PAAD.repertoire.tumor,markers)

##vital_status
clinical.patient.tumor$vital_status<-factor(clinical.patient.tumor$vital_status)
association.test.immuneRep(clinical.patient.tumor,"vital_status",PAAD.repertoire.tumor,markers)

##Smoking
clinical.patient.tumor$smoking<-factor(ifelse(clinical.patient.tumor$tobacco_smoking_history_master=="Current smoker (includes daily smokers and non-daily smokers or occasional smokers)","Current",
                                       ifelse(clinical.patient.tumor$tobacco_smoking_history_master=="Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime)","Non-smoker","Former")))
association.test.immuneRep(clinical.patient.tumor,"smoking",PAAD.repertoire.tumor,markers)

##number_pack_years_smoked
association.test.immuneRep(clinical.patient.tumor,"number_pack_years_smoked",PAAD.repertoire.tumor,markers)

##Alcohol
clinical.patient.tumor$alcohol_history_documented<-factor(clinical.patient.tumor$alcohol_history_documented)
association.test.immuneRep(clinical.patient.tumor,"alcohol_history_documented",PAAD.repertoire.tumor,markers)

##Alcohol category
clinical.patient.tumor$alcoholic_exposure_category2<-ifelse(clinical.patient.tumor$alcohol_history_documented=="NO","No-drinker",
                                                     ifelse(clinical.patient.tumor$alcohol_history_documented=="YES" & clinical.patient.tumor$alcoholic_exposure_category=="",NA,
                                                     ifelse(clinical.patient.tumor$alcoholic_exposure_category=="None","None-Drinker",
                                                     ifelse(clinical.patient.tumor$alcoholic_exposure_category=="Occasional Drinker","Occasional-Drinker",
                                                     ifelse(clinical.patient.tumor$alcoholic_exposure_category=="Daily Drinker","Daily-Drinker",
                                                     ifelse(clinical.patient.tumor$alcoholic_exposure_category=="Social Drinker","Social-Drinker",
                                                     ifelse(clinical.patient.tumor$alcoholic_exposure_category=="Weekly Drinker","Weekly-Drinker",NA)))))))
clinical.patient.tumor$alcoholic_exposure_category2<-factor(clinical.patient.tumor$alcoholic_exposure_category2)
association.test.immuneRep(clinical.patient.tumor,"alcoholic_exposure_category2",PAAD.repertoire.tumor,markers)

##family history
clinical.patient.tumor$family_history_of_cancer<-factor(clinical.patient.tumor$family_history_of_cancer)
association.test.immuneRep(clinical.patient.tumor,"family_history_of_cancer",PAAD.repertoire.tumor,markers)

##radiation_therapy
clinical.patient.tumor$radiation_therapy<-factor(clinical.patient.tumor$radiation_therapy)
association.test.immuneRep(clinical.patient.tumor,"radiation_therapy",PAAD.repertoire.tumor,markers)

##primary_therapy_outcome_success
clinical.patient.tumor$primary_therapy_outcome_success<-factor(clinical.patient.tumor$primary_therapy_outcome_success)
association.test.immuneRep(clinical.patient.tumor,"primary_therapy_outcome_success",PAAD.repertoire.tumor,markers)

##history_chronic_pancreatitis
clinical.patient.tumor$history_of_chronic_pancreatitis<-factor(clinical.patient.tumor$history_of_chronic_pancreatitis)
association.test.immuneRep(clinical.patient.tumor,"history_of_chronic_pancreatitis",PAAD.repertoire.tumor,markers)

##Function to run the association between clinical outcome and BCR/TCR
association.test.immuneRep<- function (clinical.patient.tumor,clinical.var,PAAD.repertoire.tumor,markers){
  if (class(clinical.patient.tumor[,clinical.var])=="factor"){
    p.markers<-NULL
    for(i in 1:length(markers)){
      p.markers[i]<-kruskal.test(PAAD.repertoire.tumor[which(clinical.patient.tumor[,clinical.var]!=""),markers[i]]~
                                   clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var])$p.value
    }
    dfplot <- data.frame(PAAD.repertoire.tumor[which(clinical.patient.tumor[,clinical.var]!=""),markers[which(p.markers<0.05)]],
                         clinical.var = clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var])
    colnames(dfplot)[ncol(dfplot)]<- clinical.var
    
    if (dim(dfplot)[2]>1){
      for(i in 1:length(markers[which(p.markers<0.05)])){
        print(i)
        dfplot$marker<-PAAD.repertoire.tumor[which(clinical.patient.tumor[,clinical.var]!=""),markers[which(p.markers<0.05)][i]]
        tiff(paste0("Results/boxplot_",clinical.var,"_",markers[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
        print(ggplot(dfplot,aes(y=marker, fill=clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var], 
              x=clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var])) + geom_boxplot() 
              + scale_y_continuous(name= markers[which(p.markers<0.05)][i]) + stat_compare_means()  + scale_x_discrete(name=clinical.var)
              + scale_fill_discrete(name=clinical.var))
        dev.off()
      }
    }
  } else {
    p.markers<-NULL
    for(i in 1:length(markers)){
      p.markers[i]<-coef(summary(glm(PAAD.repertoire.tumor[which(clinical.patient.tumor[,clinical.var]!=""),markers[i]]~
                                       clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var])))[2,4]
    }
    dfplot <- data.frame(marker =PAAD.repertoire.tumor[which(clinical.patient.tumor[,clinical.var]!=""),markers[which(p.markers<0.05)]],
                         clinical.var = clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var])
    colnames(dfplot)[ncol(dfplot)]<- clinical.var
    
    if (dim(dfplot)[2]>1){
      for(i in 1:length(markers[which(p.markers<0.05)])){
        print(i)
        dfplot$marker<-PAAD.repertoire.tumor[which(clinical.patient.tumor[,clinical.var]!=""),markers[which(p.markers<0.05)][i]]
        tiff(paste0("Results/boxplot_",clinical.var,"_",markers[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
        print(ggplot(dfplot,aes(y=marker, fill=clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var], 
                                x=clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var])) + geom_point() + geom_smooth(method='lm')
              + scale_y_continuous(name= markers[which(p.markers<0.05)][i]) + stat_cor(method = "pearson") + scale_x_continuous(name=clinical.var)
              + scale_fill_continuous(name=clinical.var))
        dev.off()
      }
    }
  }
  return("Done")
}


#############
##clinical outcome with xcell
############
clinical.patient.xcell<-clinical.patient[match(substr(rownames(xcell.data.tumor.filter),1,12),clinical.patient$bcr_patient_barcode),]

##Function to run the association between clinical outcome and xcell
association.test.xcell<- function (clinical.patient.xcell,clinical.var,xcell.data.tumor.filter){
  if (class(clinical.patient.xcell[,clinical.var])=="factor"){
    p.markers<-NULL
    for(i in 1: ncol(xcell.data.tumor.filter)){
      p.markers[i]<-kruskal.test(xcell.data.tumor.filter[which(clinical.patient.xcell[,clinical.var]!=""),i]~
                                   clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var])$p.value
    }
    dfplot <- data.frame(xcell.data.tumor.filter[which(clinical.patient.xcell[,clinical.var]!=""),which(p.markers<0.05)],
                         clinical.var = clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var])
    colnames(dfplot)[ncol(dfplot)]<- clinical.var
    
    if (dim(dfplot)[2]>1){
      for(i in 1:length(which(p.markers<0.05))){
        print(i)
        dfplot$marker<-xcell.data.tumor.filter[which(clinical.patient.xcell[,clinical.var]!=""),which(p.markers<0.05)[i]]
        tiff(paste0("Results/boxplot_",clinical.var,"_",colnames(xcell.data.tumor.filter)[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
        print(ggplot(dfplot,aes(y=marker, fill=clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var], 
                                x=clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var])) + geom_boxplot() 
              + scale_y_continuous(name=colnames(xcell.data.tumor.filter)[which(p.markers<0.05)][i]) + stat_compare_means()  + scale_x_discrete(name=clinical.var)
              + scale_fill_discrete(name=clinical.var))
        dev.off()
      }
    }
  } else {
    p.markers<-NULL
    for(i in 1:ncol(xcell.data.tumor.filter)){
      p.markers[i]<-coef(summary(glm(xcell.data.tumor.filter[which(clinical.patient.xcell[,clinical.var]!=""),i]~
                                       clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var])))[2,4]
    }
    dfplot <- data.frame(marker =xcell.data.tumor.filter[which(clinical.patient.xcell[,clinical.var]!=""),which(p.markers<0.05)],
                         clinical.var = clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var])
    colnames(dfplot)[ncol(dfplot)]<- clinical.var
    
    if (dim(dfplot)[2]>1){
      for(i in 1:length(which(p.markers<0.05))){
        print(i)
        dfplot$marker<-xcell.data.tumor.filter[which(clinical.patient.xcell[,clinical.var]!=""),which(p.markers<0.05)[i]]
        tiff(paste0("Results/boxplot_",clinical.var,"_",colnames(xcell.data.tumor.filter)[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
        print(ggplot(dfplot,aes(y=marker, fill=clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var], 
                                x=clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var])) + geom_point() + geom_smooth(method='lm')
              + scale_y_continuous(name= colnames(xcell.data.tumor.filter)[which(p.markers<0.05)][i]) + stat_cor(method = "pearson") + scale_x_continuous(name=clinical.var)
              + scale_fill_continuous(name=clinical.var))
        dev.off()
      }
    }
  }
  return("Done")
}

##Histological type
clinical.patient.xcell$histological_type_2cat<-factor(ifelse(clinical.patient.xcell$histological_type=="Pancreas-Adenocarcinoma Ductal Type","PDAC","Other"))
association.test.xcell(clinical.patient.xcell,"histological_type_2cat",xcell.data.tumor.filter)


##anatomic_neoplasm_subdivision
clinical.patient.xcell$anatomic_neoplasm_subdivision<-factor(clinical.patient.xcell$anatomic_neoplasm_subdivision)
association.test.xcell(clinical.patient.xcell,"anatomic_neoplasm_subdivision",xcell.data.tumor.filter)

##gender
clinical.patient.xcell$gender<-factor(clinical.patient.xcell$gender)
association.test.xcell(clinical.patient.xcell,"gender",xcell.data.tumor.filter)

##race_list
clinical.patient.xcell$race_list<-factor(clinical.patient.xcell$race_list)
association.test.xcell(clinical.patient.xcell,"race_list",xcell.data.tumor.filter)

##History of Prior Malignancy
clinical.patient.xcell$other_dx<-factor(clinical.patient.xcell$other_dx)
association.test.xcell(clinical.patient.xcell,"other_dx",xcell.data.tumor.filter)

##lymph_node_examined_count
association.test.xcell(clinical.patient.xcell,"lymph_node_examined_count",xcell.data.tumor.filter)

##neoplasm_histologic_grade
clinical.patient.xcell$neoplasm_histologic_grade_3cat<-factor(ifelse(clinical.patient.xcell$neoplasm_histologic_grade=="G1","G1",
                                                                     ifelse(clinical.patient.xcell$neoplasm_histologic_grade=="G2","G2",
                                                                            ifelse(clinical.patient.xcell$neoplasm_histologic_grade=="G3","G3",""))))
association.test.xcell(clinical.patient.xcell,"neoplasm_histologic_grade_3cat",xcell.data.tumor.filter)

##Age 
association.test.xcell(clinical.patient.xcell,"age_at_initial_pathologic_diagnosis",xcell.data.tumor.filter)

##vital_status
clinical.patient.xcell$vital_status<-factor(clinical.patient.xcell$vital_status)
association.test.xcell(clinical.patient.xcell,"vital_status",xcell.data.tumor.filter)

## Heatmap ####
annotation_col = data.frame(
  vital_status = factor(clinical.patient.xcell$vital_status))
rownames(annotation_col)<-rownames(xcell.data.tumor.filter)
xcell.data.tumor.filter_vital_status<-xcell.data.tumor.filter[,c("aDC","Adipocytes","CD4+ naive T-cells","CD8+ T-cells","cDC",
                                                                 "DC","Epithelial cells","Fibroblasts","HSC","ImmuneScore",
                                                                 "Keratinocytes","Melanocytes","MicroenvironmentScore","Monocytes",
                                                                 "NKT","StromaScore")]

ann_colors = list (vital_status = c("Alive" = brewer.pal(3,"Set2")[1], "Dead"= brewer.pal(3,"Set2")[2]))
pheatmap(t(xcell.data.tumor.filter_vital_status),scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12),
         annotation_col = annotation_col,annotation_colors = ann_colors)



##Smoking
clinical.patient.xcell$smoking<-factor(ifelse(clinical.patient.xcell$tobacco_smoking_history_master=="Current smoker (includes daily smokers and non-daily smokers or occasional smokers)","Current",
                                              ifelse(clinical.patient.xcell$tobacco_smoking_history_master=="Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime)","Non-smoker","Former")))
association.test.xcell(clinical.patient.xcell,"smoking",xcell.data.tumor.filter)

##number_pack_years_smoked
association.test.xcell(clinical.patient.xcell,"number_pack_years_smoked",xcell.data.tumor.filter)

##Alcohol
clinical.patient.xcell$alcohol_history_documented<-factor(clinical.patient.xcell$alcohol_history_documented)
association.test.xcell(clinical.patient.xcell,"alcohol_history_documented",xcell.data.tumor.filter)

##Alcohol category
clinical.patient.xcell$alcoholic_exposure_category2<-ifelse(clinical.patient.xcell$alcohol_history_documented=="NO","No-drinker",
                                                            ifelse(clinical.patient.xcell$alcohol_history_documented=="YES" & clinical.patient.xcell$alcoholic_exposure_category=="",NA,
                                                                   ifelse(clinical.patient.xcell$alcoholic_exposure_category=="None","None-Drinker",
                                                                          ifelse(clinical.patient.xcell$alcoholic_exposure_category=="Occasional Drinker","Occasional-Drinker",
                                                                                 ifelse(clinical.patient.xcell$alcoholic_exposure_category=="Daily Drinker","Daily-Drinker",
                                                                                        ifelse(clinical.patient.xcell$alcoholic_exposure_category=="Social Drinker","Social-Drinker",
                                                                                               ifelse(clinical.patient.xcell$alcoholic_exposure_category=="Weekly Drinker","Weekly-Drinker",NA)))))))
clinical.patient.xcell$alcoholic_exposure_category2<-factor(clinical.patient.xcell$alcoholic_exposure_category2)
association.test.xcell(clinical.patient.xcell,"alcoholic_exposure_category2",xcell.data.tumor.filter)

##family history
clinical.patient.xcell$family_history_of_cancer<-factor(clinical.patient.xcell$family_history_of_cancer)
association.test.xcell(clinical.patient.xcell,"family_history_of_cancer",xcell.data.tumor.filter)

##radiation_therapy
clinical.patient.xcell$radiation_therapy<-factor(clinical.patient.xcell$radiation_therapy)
association.test.xcell(clinical.patient.xcell,"radiation_therapy",xcell.data.tumor.filter)

##primary_therapy_outcome_success
clinical.patient.xcell$primary_therapy_outcome_success<-factor(clinical.patient.xcell$primary_therapy_outcome_success)
association.test.xcell(clinical.patient.xcell,"primary_therapy_outcome_success",xcell.data.tumor.filter)

##history_chronic_pancreatitis
clinical.patient.xcell$history_of_chronic_pancreatitis<-factor(clinical.patient.xcell$history_of_chronic_pancreatitis)
association.test.xcell(clinical.patient.xcell,"history_of_chronic_pancreatitis",xcell.data.tumor.filter)
