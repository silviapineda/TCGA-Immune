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
PAAD.GTEx.repertoire.diversity.tumor.normmal<-PAAD.GTEx.repertoire.diversity[which(PAAD.GTEx.repertoire.diversity$outcome=="normal_pancreas (GTEx)"|
                                                                                     PAAD.GTEx.repertoire.diversity$outcome== "PDAC (TCGA)"),]
PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome<-factor(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome)

###########################################
##### Analysis with V gene usage #########
##########################################
Vgene_differences<-function(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge){
  ####Filter the genes that has low clonanlity (clones<100)
  assign(paste0("PAAD.GTEx.repertoire.diversity.tumor.normmal_",chain),
         PAAD.GTEx.repertoire.diversity.tumor.normmal[which(PAAD.GTEx.repertoire.diversity.tumor.normmal$clones_IGH>10),])
  id<-match(data_merge$sample,rownames(get(paste0("PAAD.GTEx.repertoire.diversity.tumor.normmal_",chain))))
  data_merge_qc<-data_merge[which(is.na(id)==F),]

  ###Matrix with the vgenes
  vgenes<-as.data.frame(unclass(table(data_merge_qc$sample,data_merge_qc$bestVGene)))
  vgenes<-vgenes[,-1]
  ###Genes 
  #Obtain the V usage
  v_usage<-matrix(NA,nrow(vgenes),ncol(vgenes))
  id<-NULL
  for (i in 1:ncol(vgenes)){
    if(substr(colnames(vgenes)[i],1,3)==chain) {
      v_usage[,i]<-vgenes[,i]/get(paste0("PAAD.GTEx.repertoire.diversity.tumor.normmal_",chain))$clones_IGH
      id<-c(id,i)
    }
  }
  colnames(v_usage)<-colnames(vgenes)
  rownames(v_usage)<-rownames(vgenes)
  v_usage<-v_usage[,id]

  ###FIlTERING
  ###Convert into 0 all those who has a low expression (<0.002)
  #xx<-replace(v_usage_filter,v_usage_filter<0.05,0)
  ###Those who are in lesss than 10%
  v_usage_filter<-v_usage[,which(apply(v_usage,2,function(x) sum(x==0))<dim(v_usage)[1]*0.9)]
 
  ###Applied ENET 
  alphalist<-seq(0.1,0.9,by=0.1)
  set.seed(54)
  elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(v_usage_filter,get(paste0("PAAD.GTEx.repertoire.diversity.tumor.normmal_",chain))$outcome,family="binomial"
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
  enet<-glmnet(v_usage_filter,get(paste0("PAAD.GTEx.repertoire.diversity.tumor.normmal_",chain))$outcome,family="binomial",standardize=TRUE,alpha=alpha,lambda=lambda) #86
  genes<-rownames(enet$beta)[which( as.numeric(enet$beta)!=0)]

  v_usage_sign<-v_usage_filter[,match(genes,colnames(v_usage_filter))]

  ##Plot results
  brewer.pal(4,name = "Accent")
  cols=c( "#7FC97F","#BEAED4")

  annotation_row = data.frame(get(paste0("PAAD.GTEx.repertoire.diversity.tumor.normmal_",chain))$outcome)
  ann_colors = list (outcome = c("normal_pancreas (GTEx)" = cols[1],"PDAC (TCGA)" = cols[2]))
  colnames(annotation_row)<-"outcome"
  rownames(annotation_row)<-rownames(v_usage_sign)

  tiff(paste0("Results/Vgene/heatmap_VgeneUsage_",chain,".tiff"),width = 5000, height = 3000, res = 300)
  pheatmap(t(v_usage_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(9,name="PuOr"))(10))
  dev.off()
  return(v_usage_sign)
}

####Chain
#IGH
chain="IGH"
v_usage_sign_IGH<-Vgene_differences(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)

#IGK
chain="IGK"
v_usage_sign_IGK<-Vgene_differences(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)

#IGL
chain="IGL"
v_usage_sign_IGL<-Vgene_differences(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)

#TRA
chain="TRA"
v_usage_sign_TRA<-Vgene_differences(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)

#TRB
chain="TRB"
v_usage_sign_TRB<-Vgene_differences(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)

#TRD
chain="TRD"
v_usage_sign_TRD<-Vgene_differences(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)

#TRB
chain="TRG"
v_usage_sign_TRG<-Vgene_differences(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)


################################################
##### Analysis with V gene expression #########
###############################################
#Read the total reads to extract the gene expression
#Obtain the V expression
Vgene_differences_expression<-function(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge){
  ####Filter the genes that has low clonanlity (clones<100)
  assign(paste0("PAAD.GTEx.repertoire.diversity.tumor.normmal_",chain),
         PAAD.GTEx.repertoire.diversity.tumor.normmal[which(PAAD.GTEx.repertoire.diversity.tumor.normmal$clones_IGH>10),])
  id<-match(data_merge$sample,rownames(get(paste0("PAAD.GTEx.repertoire.diversity.tumor.normmal_",chain))))
  data_merge_qc<-data_merge[which(is.na(id)==F),]
  
  ###Matrix with the vgenes
  vgenes<-as.data.frame(unclass(table(data_merge_qc$sample,data_merge_qc$bestVGene)))
  vgenes<-vgenes[,-1]
  ###Genes 
  #Obtain the V usage
  v_expression<-matrix(NA,nrow(vgenes),ncol(vgenes))
  id<-NULL
  for (i in 1:ncol(vgenes)){
    if(substr(colnames(vgenes)[i],1,3)==chain) {
         v_expression[,i]<-100*(vgenes[,i]/get(paste0("PAAD.GTEx.repertoire.diversity.tumor.normmal_",chain))$totalReads)
        id<-c(id,i)
             
    } 
  }
  colnames(v_expression)<-colnames(vgenes)
  rownames(v_expression)<-rownames(vgenes)
  v_expression<-v_expression[,id]
  
  ###FIlTERING
  ###Convert into 0 all those who has a low expression (<0.002)
  #xx<-replace(v_usage_filter,v_usage_filter<0.05,0)
  ###Those who are in lesss than 10%
  v_expression_filter<-v_expression[,which(apply(v_expression,2,function(x) sum(x==0))<dim(v_expression)[1]*0.9)]
  
  ###Applied ENET 
  alphalist<-seq(0.1,0.9,by=0.1)
  set.seed(54)
  elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(v_expression_filter,get(paste0("PAAD.GTEx.repertoire.diversity.tumor.normmal_",chain))$outcome,family="binomial"
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
  enet<-glmnet(v_expression_filter,get(paste0("PAAD.GTEx.repertoire.diversity.tumor.normmal_",chain))$outcome,family="binomial",standardize=TRUE,alpha=alpha,lambda=lambda) #86
  genes<-rownames(enet$beta)[which( as.numeric(enet$beta)!=0)]
  
  v_expression_sign<-v_expression_filter[,match(genes,colnames(v_expression_filter))]
  
  ##Plot results
  brewer.pal(4,name = "Accent")
  cols=c( "#7FC97F","#BEAED4")
  
  annotation_row = data.frame(get(paste0("PAAD.GTEx.repertoire.diversity.tumor.normmal_",chain))$outcome)
  ann_colors = list (outcome = c("normal_pancreas (GTEx)" = cols[1],"PDAC (TCGA)" = cols[2]))
  colnames(annotation_row)<-"outcome"
  rownames(annotation_row)<-rownames(v_expression_sign)
  
  tiff(paste0("Results/Vgene/heatmap_VgeneExpression_",chain,".tiff"),width = 5000, height = 3000, res = 300)
  pheatmap(t(v_expression_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
           annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(9,name="PuOr"))(10))
  dev.off()
  return(v_expression_sign)
}


####Chain
#IGH
chain="IGH"
v_expression_sign_IGH<-Vgene_differences_expression(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)

#IGK
chain="IGK"
v_expression_sign_IGK<-Vgene_differences_expression(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)

#IGL
chain="IGL"
v_expression_sign_IGL<-Vgene_differences_expression(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)

#TRA
chain="TRA"
v_expression_sign_TRA<-Vgene_differences_expression(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)

#TRB
chain="TRB"
v_expression_sign_TRB<-Vgene_differences_expression(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)

#TRD
chain="TRD"
v_expression_sign_TRD<-Vgene_differences_expression(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)

#TRG
chain="TRG"
v_expression_sign_TRG<-Vgene_differences_expression(PAAD.GTEx.repertoire.diversity.tumor.normmal,chain,data_merge)
