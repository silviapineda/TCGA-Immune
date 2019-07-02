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

###Matrix with the vgenes
vgenes<-as.data.frame(unclass(table(data_merge$sample,data_merge$bestVGene)))

#Read the total reads to extract the gene expression
total_reads<-read.table("Data/PAAD/MIXCR_PAAD/total_reads.txt",sep=";")
id.reads<-match(rownames(vgenes),total_reads$V1)

#Obtain the V expression
v_expression<-matrix(NA,nrow(vgenes),ncol(vgenes))
for (i in 1:ncol(vgenes)){
  v_expression[,i]<-100000*(vgenes[,i]/total_reads$V2[id.reads])
}
colnames(v_expression)<-colnames(vgenes)
rownames(v_expression)<-rownames(vgenes)

##Only pancreas and normal
v_expression<-v_expression[which(is.na(PAAD.repertoire.diversity$Tumor_type_2categ[id.spec])==F),]
id.spec<-match(rownames(v_expression),rownames(PAAD.repertoire.diversity))

###FIlTERING
###Those who are in lesss than 20%
v_expression_filter<-v_expression[,which(apply(v_expression,2,function(x) sum(x==0))<=160-160*0.2)] ##212 genes
###Applied ENET to find genes associated with the IRF9_expression
alphalist<-seq(0.01,0.99,by=0.01)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(v_expression_filter,PAAD.repertoire.diversity$Tumor_type_2categ[id.spec],family="binomial"
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
enet<-glmnet(v_expression_filter,PAAD.repertoire.diversity$Tumor_type_2categ[id.spec],family="binomial",standardize=TRUE,alpha=alpha,lambda=lambda) #86

# p_value=NULL
# for(i in 1:dim(v_expression_filter)[2]){
#   print(i)
#   mod<-glm(PAAD.repertoire.diversity$Tumor_type_2categ[id.spec]~v_expression_filter[,i],family = "binomial")
#   p_value[i]=coefficients(summary(mod))[2,4]
# }
# v_expression_filter_sign<-v_expression_filter[,which(p_value<0.05)]
##Extract genes that are selected by ENET

genes<-rownames(enet$beta)[which(enet$beta!=0)]
v_expression_filter_sign<-v_expression_filter[,match(genes,colnames(v_expression_filter))]

##Plot results
annotation_row = data.frame(PAAD.repertoire.diversity$Tumor_type_2categ)
colnames(annotation_row)<-"Tumor_type_2categ"
rownames(annotation_row)<-rownames(PAAD.repertoire.diversity)
ann_colors = list (Tumor_type_2categ = c("normal_pseudonormal_pancreas" = brewer.pal(3,"Accent")[1],
                                         "Tumor_pancreas"= brewer.pal(3,"Accent")[2]))

tiff("Results/heatmap_VgeneUsage_sign.tiff",width = 5000, height = 3000, res = 300)
pheatmap(t(v_expression_filter_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))
dev.off()


