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

##Colors: brewer.pal(4, "Accent")
##"#BEAED4" (IGH) "#7FC97F" (IGK) "#FDC086" (IGL)
setwd("~/TCGA-Immune/")

load("Data/PAAD/PAAD_FullData.Rdata")

##Filter for those that are pancreas or normal-pseudonormal
data_merge_qc<-data_merge[which(is.na(match(data_merge$sample,rownames(PAAD.repertoire.diversity)))==F),]


###########################################
##### Analysis with V gene usage #########
##########################################

###Matrix with the vgenes
vgenes<-as.data.frame(unclass(table(data_merge_qc$sample,data_merge_qc$bestVGene)))
vgenes<-vgenes[,-1]
###Genes 
#Obtain the V usage
v_usage<-matrix(NA,nrow(vgenes),ncol(vgenes))
for (i in 1:ncol(vgenes)){
  if(substr(colnames(vgenes)[i],1,3)=="IGH"){
    v_usage[,i]<-vgenes[,i]/PAAD.repertoire.diversity$clones_IGH
  } 
  if(substr(colnames(vgenes)[i],1,3)=="IGK"){
    v_usage[,i]<-vgenes[,i]/PAAD.repertoire.diversity$clones_IGK
  } 
  if(substr(colnames(vgenes)[i],1,3)=="IGL"){
    v_usage[,i]<-vgenes[,i]/PAAD.repertoire.diversity$clones_IGL
  } 
  if(substr(colnames(vgenes)[i],1,3)=="TRA"){
    v_usage[,i]<-vgenes[,i]/PAAD.repertoire.diversity$clones_TRA
  } 
  if(substr(colnames(vgenes)[i],1,3)=="TRB"){
    v_usage[,i]<-vgenes[,i]/PAAD.repertoire.diversity$clones_TRB
  } 
  if(substr(colnames(vgenes)[i],1,3)=="TRD"){
    v_usage[,i]<-vgenes[,i]/PAAD.repertoire.diversity$clones_TRD
  } 
  if(substr(colnames(vgenes)[i],1,3)=="TRG"){
    v_usage[,i]<-vgenes[,i]/PAAD.repertoire.diversity$clones_TRG
  } 
  
}
colnames(v_usage)<-colnames(vgenes)
rownames(v_usage)<-rownames(vgenes)
v_usage<-t(apply(v_usage,1,function(x) replace(x,x=="NaN",0)))

###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
#xx<-replace(v_usage_filter,v_usage_filter<0.05,0)
###Those who are in lesss than 10%
v_usage_filter<-v_usage[,which(apply(v_usage,2,function(x) sum(x==0))<157)] #259

###Applied ENET to find genes associated with the IRF9_expression
alphalist<-seq(0.01,0.99,by=0.01)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(v_usage_filter,PAAD.repertoire.diversity$Tumor_type_2categ,family="binomial"
                                                        ,standardize=TRUE,alpha=a,nfolds=5))})
xx<-rep(NA,length(alphalist))
yy<-rep(NA,length(alphalist))
for (j in 1:length(alphalist)) {
  print(j)
  if(class(elasticnet[[j]]) != "try-error"){
    xx[j]<-elasticnet[[j]]$lambda.min
    id.cv.opt<-grep(elasticnet[[j]]$lambda.min,elasticnet[[j]]$lambda,fixed=TRUE)
    yy[j]<-elasticnet[[j]]$cvm[id.cv.opt]
  }
}
id.min<-which(yy==min(yy,na.rm=TRUE))
lambda<-xx[id.min]
alpha<-alphalist[id.min]
enet<-glmnet(v_usage_filter,PAAD.repertoire.diversity$Tumor_type_2categ,family="binomial",standardize=TRUE,alpha=alpha,lambda=lambda) #86

p_value=NULL
for(i in 1:dim(v_usage_filter)[2]){
  print(i)
  mod<-glm(PAAD.repertoire.diversity$Tumor_type_2categ~v_usage_filter[,i],family = "binomial")
  p_value[i]=coefficients(summary(mod))[2,4]
}
v_usage_filter_sign<-v_expression_filter[,which(p_value<0.05)]

##Plot results
annotation_row = data.frame(PAAD.repertoire.diversity$Tumor_type_2categ)
colnames(annotation_row)<-"Tumor_type_2categ"
rownames(annotation_row)<-rownames(PAAD.repertoire.diversity)
ann_colors = list (Tumor_type_2categ = c("normal_pseudonormal_pancreas" = brewer.pal(3,"Accent")[1],
                                         "Tumor_pancreas"= brewer.pal(3,"Accent")[2]))

tiff("Results/heatmap_VgeneUsage_sign.tiff",width = 5000, height = 3000, res = 300)
pheatmap(t(v_usage_filter_sign),border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(4,name="Reds"))(1000))
dev.off()


################################################
##### Analysis with V gene expression #########
###############################################
#Read the total reads to extract the gene expression
total_reads<-read.table("Data/PAAD/MIXCR_PAAD/total_reads.txt",sep=";")
id.reads<-match(rownames(vgenes),total_reads$V1)

#Obtain the V expression
v_expression<-matrix(NA,nrow(vgenes),ncol(vgenes))
for (i in 1:ncol(vgenes)){
  v_expression[,i]<-100*(vgenes[,i]/total_reads$V2[id.reads])
}

colnames(v_expression)<-colnames(vgenes)
rownames(v_expression)<-rownames(vgenes)

###FIlTERING
###Those who are at least in 2 samples
v_expression_filter<-v_expression[,which(apply(v_expression,2,function(x) sum(x==0))<157)] #259

####GLM 
p_value=NULL
for(i in 1:dim(v_expression_filter)[2]){
  print(i)
  mod<-glm(PAAD.repertoire.diversity$Tumor_type_2categ~v_expression_filter[,i],family = "binomial")
  p_value[i]=coefficients(summary(mod))[2,4]
}
v_expression_filter_sign<-v_expression_filter[,which(p_value<0.05)]

##Plot results
annotation_row = data.frame(PAAD.repertoire.diversity$Tumor_type_2categ)
colnames(annotation_row)<-"Tumor_type_2categ"
rownames(annotation_row)<-rownames(PAAD.repertoire.diversity)
ann_colors = list (Tumor_type_2categ = c("normal_pseudonormal_pancreas" = brewer.pal(3,"Accent")[1],
                                         "Tumor_pancreas"= brewer.pal(3,"Accent")[2]))

tiff("Results/heatmap_VgeneExpression_sign.tiff",width = 5000, height = 3000, res = 300)
pheatmap(t(v_expression_filter_sign),border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(4,name="Reds"))(1000))
dev.off()

