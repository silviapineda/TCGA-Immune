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

load("Data/PAAD_GTEx/PAAD_GTEx_FullData.Rdata")

##Restrict only to Normal vs Tumor 
PAAD.GTEx.repertoire.diversity.tumor.normmal<-PAAD.GTEx.repertoire.diversity[which(PAAD.GTEx.repertoire.diversity$outcome=="normal-pancreas (GTEx)"|
                                                                                     PAAD.GTEx.repertoire.diversity$outcome== "tumor-pancreas (TCGA)"),]
PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome<-factor(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome)

id<-match(data_merge$sample,rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal))
data_merge_qc<-data_merge[which(is.na(id)==F),]

###########################################
##### Analysis with V gene usage #########
##########################################

###Matrix with the vgenes
vgenes<-as.data.frame(unclass(table(data_merge_qc$sample,data_merge_qc$bestVGene)))
vgenes<-vgenes[,-1]
vgenes_Ig
###Genes 
#Obtain the V usage
v_usage_Ig<-matrix(NA,nrow(vgenes),ncol(vgenes))
v_usage_TCR<-matrix(NA,nrow(vgenes),ncol(vgenes))
id.Ig<-NULL
id.TR<-NULL
for (i in 1:ncol(vgenes)){
  if(substr(colnames(vgenes)[i],1,3)=="IGH"){
    v_usage_Ig[,i]<-vgenes[,i]/PAAD.GTEx.repertoire.diversity.tumor.normmal$clones_IGH
    id.Ig<-c(id.Ig,i)
  } 
  if(substr(colnames(vgenes)[i],1,3)=="IGK"){
    v_usage_Ig[,i]<-vgenes[,i]/PAAD.GTEx.repertoire.diversity.tumor.normmal$clones_IGK
    id.Ig<-c(id.Ig,i)
  } 
  if(substr(colnames(vgenes)[i],1,3)=="IGL"){
    v_usage_Ig[,i]<-vgenes[,i]/PAAD.GTEx.repertoire.diversity.tumor.normmal$clones_IGL
    id.Ig<-c(id.Ig,i)
  } 
  if(substr(colnames(vgenes)[i],1,3)=="TRA"){
    v_usage_TCR[,i]<-vgenes[,i]/PAAD.GTEx.repertoire.diversity.tumor.normmal$clones_TRA
    id.TR<-c(id.TR,i)
  } 
  if(substr(colnames(vgenes)[i],1,3)=="TRB"){
    v_usage_TCR[,i]<-vgenes[,i]/PAAD.GTEx.repertoire.diversity.tumor.normmal$clones_TRB
    id.TR<-c(id.TR,i)
  } 
  if(substr(colnames(vgenes)[i],1,3)=="TRD"){
    v_usage_TCR[,i]<-vgenes[,i]/PAAD.GTEx.repertoire.diversity.tumor.normmal$clones_TRD
    id.TR<-c(id.TR,i)
  } 
  if(substr(colnames(vgenes)[i],1,3)=="TRG"){
    v_usage_TCR[,i]<-vgenes[,i]/PAAD.GTEx.repertoire.diversity.tumor.normmal$clones_TRG
    id.TR<-c(id.TR,i)
  } 
  
}
colnames(v_usage_Ig)<-colnames(vgenes)
rownames(v_usage_Ig)<-rownames(vgenes)
colnames(v_usage_TCR)<-colnames(vgenes)
rownames(v_usage_TCR)<-rownames(vgenes)

v_usage_TCR<-v_usage_TCR[,id.TR]
v_usage_Ig<-v_usage_Ig[,id.Ig]

v_usage_Ig<-t(apply(v_usage_Ig,1,function(x) replace(x,x=="NaN",0)))
v_usage_TCR<-t(apply(v_usage_TCR,1,function(x) replace(x,x=="NaN",0)))
###FIlTERING
###Convert into 0 all those who has a low expression (<0.002)
#xx<-replace(v_usage_filter,v_usage_filter<0.05,0)
###Those who are in lesss than 10%
v_usage_Ig_filter<-v_usage_Ig[,which(apply(v_usage_Ig,2,function(x) sum(x==0))<292)] #211
v_usage_TCR_filter<-v_usage_TCR[,which(apply(v_usage_TCR,2,function(x) sum(x==0))<292)] #211

##IG
###Applied ENET 
alphalist<-seq(0.1,0.9,by=0.1)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(v_usage_Ig_filter,PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome,family="binomial"
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
enet<-glmnet(v_usage_Ig_filter,PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome,family="binomial",standardize=TRUE,alpha=alpha,lambda=lambda) #86
genes<-rownames(enet$beta)[which(enet$beta!=0)]

v_usage_Ig_filter_sign<-v_usage_Ig_filter[,match(genes,colnames(v_usage_Ig_filter))] #51

##Plot results
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F","#BEAED4")

annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome)
ann_colors = list (outcome = c("normal-pancreas (GTEx)" = cols[1],"tumor-pancreas (TCGA)" = cols[2]))
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(v_usage_filter_sign)

tiff("Results/heatmap_VgeneUsage_Ig_sign.tiff",width = 5000, height = 3000, res = 300)
pheatmap(t(v_usage_Ig_filter_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(9,name="PuOr"))(10))
dev.off()


##TCR
###Applied ENET
alphalist<-seq(0.1,0.9,by=0.1)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(v_usage_TCR_filter,PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome,family="binomial"
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
enet<-glmnet(v_usage_TCR_filter,PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome,family="binomial",standardize=TRUE,alpha=alpha,lambda=lambda) #86
genes<-rownames(enet$beta)[which(enet$beta!=0)]

v_usage_TCR_filter_sign<-v_usage_TCR_filter[,match(genes,colnames(v_usage_TCR_filter))] #51

##Plot results
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F","#BEAED4")

annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome)
ann_colors = list (outcome = c("normal-pancreas (GTEx)" = cols[1],"tumor-pancreas (TCGA)" = cols[2]))
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(v_usage_TCR_filter_sign)

tiff("Results/heatmap_VgeneUsage_TCR_sign.tiff",width = 5000, height = 3000, res = 300)
pheatmap(t(v_usage_TCR_filter_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(9,name="PuOr"))(10))
dev.off()



################################################
##### Analysis with V gene expression #########
###############################################
#Read the total reads to extract the gene expression
#Obtain the V expression
v_expression<-matrix(NA,nrow(vgenes),ncol(vgenes))
for (i in 1:ncol(vgenes)){
  v_expression[,i]<-100*(vgenes[,i]/PAAD.GTEx.repertoire.diversity.tumor.normmal$totalReads)
}

colnames(v_expression)<-colnames(vgenes)
rownames(v_expression)<-rownames(vgenes)

v_expression_Ig<-v_expression[,which(substr(colnames(v_expression),1,2)=="IG")]
v_expression_TR<-v_expression[,which(substr(colnames(v_expression),1,2)=="TR")]

###FIlTERING
###Those who are at least in 2 samples
v_expression_filter<-v_expression_Ig[,which(apply(v_expression_Ig,2,function(x) sum(x==0))<292)] 

##TCR
###Applied ENET
alphalist<-seq(0.1,0.9,by=0.1)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(v_expression_filter,PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome,family="binomial"
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
enet<-glmnet(v_expression_filter,PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome,family="binomial",standardize=TRUE,alpha=alpha,lambda=lambda) #86
genes<-rownames(enet$beta)[which(enet$beta!=0)]

v_expression_filter_sign<-v_expression_filter[,match(genes,colnames(v_expression_filter))] 

##Plot results
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F","#BEAED4")

annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome)
ann_colors = list (outcome = c("normal-pancreas (GTEx)" = cols[1],"tumor-pancreas (TCGA)" = cols[2]))
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(v_expression_filter_sign)

tiff("Results/heatmap_Vexpression_Ig_sign.tiff",width = 5000, height = 3000, res = 300)
pheatmap(t(v_expression_filter_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(9,name="PuOr"))(10))
dev.off()

###Overlap with IG v gene usage
intersect(colnames(v_expression_filter_sign),colnames(v_usage_Ig_filter_sign)) #48 (from 57 of v gene usage)
# "IGHV1-58"  "IGHV1-69"  "IGHV3-11"  "IGHV3-15"  "IGHV3-35"  "IGHV3-38"  "IGHV3-43"  "IGHV3-52"  "IGHV3-53"  "IGHV3-60" 
# "IGHV3-72"  "IGHV3-74"  "IGHV4-28"  "IGHV4-61"  "IGKV1-12"  "IGKV1-13"  "IGKV1-37"  "IGKV1-8"   "IGKV2-24"  "IGKV2-28" 
# "IGKV2-29"  "IGKV3-11"  "IGKV3-20"  "IGKV4-1"   "IGLV1-40"  "IGLV1-41"  "IGLV1-44"  "IGLV10-67" "IGLV2-11"  "IGLV2-23" 
# "IGLV2-34"  "IGLV2-5"   "IGLV3-1"   "IGLV3-12"  "IGLV3-16"  "IGLV3-25"  "IGLV4-69"  "IGLV5-52"  "IGLV6-57"  "IGLV7-46" 
# "IGLV9-49"  "IGHV3-47"  "IGKV1-35"  "IGHV3-25"  "IGHV7-4-1" "IGKV2-18"  "IGLV11-55" "IGKV7-3"  
 

##Overlap with TR v gene usage
intersect(colnames(v_expression_filter_sign),colnames(v_usage_TCR_filter_sign)) #13 (All of them detected by v usage)
#"TRAV12-2"  "TRBV23-1"  "TRBV6-2"   "TRBV7-3"   "TRBV5-6"   "TRGV2"     "TRAV23DV6" "TRGV10"    "TRAV29DV5" "TRAV9-2"  
# "TRBV9"     "TRAV2"     "TRAV20" 
