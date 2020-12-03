rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. Compositional Analysis for PDAC 
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP:
###         
###
### Author: Silvia Pineda
### Date: July, 2020
############################################################################################
# copy the repository from https://github.com/UVic-omics/CoDA-Penalized-Regression
#system('git clone https://github.com/UVic-omics/CoDA-Penalized-Regression')
# cran.packages <- c('knitr', 'glmnet', 'ggplot2', 'gridExtra',
#                    'UpSetR', 'ggforce')
# install.packages(cran.packages)
# devtools::install_github(repo = 'UVic-omics/selbal')

library(knitr) # rbookdown, kable
library(glmnet) # glmnet
library(selbal) # selbal
library(ggplot2) # draw selbal
library(gridExtra) # grid.arrange
library(UpSetR) # upset
library(ggforce) # selbal-like plot
library(grid) # grid.draw
# source coda-lasso functions
source(file = './CoDA-Penalized-Regression/R/functions_coda_penalized_regression.R')
library(survival)
library(survminer)
library(survMisc)

setwd("~/TCGA-Immune/")
load("Data/PAAD/PAAD_FullData.Rdata")


id<-match(data_merge$sample,rownames(PAAD.repertoire.tumor.filter))
data_merge_filter<-data_merge[which(is.na(id)==F),]


#############
## 1. Build the matrix with the clones by samples
###########
data_chain_IGK<-data_merge_filter[which(data_merge_filter$chainType=="IGK"),]
clone_type_IGK<-t(as.data.frame(unclass(table(data_chain_IGK$V_J_lenghCDR3_CloneId,factor(data_chain_IGK$sample))))) 
data_chain_IGL<-data_merge_filter[which(data_merge_filter$chainType=="IGL"),]
clone_type_IGL<-t(as.data.frame(unclass(table(data_chain_IGL$V_J_lenghCDR3_CloneId,factor(data_chain_IGL$sample))))) 
data_chain_IGH<-data_merge_filter[which(data_merge_filter$chainType=="IGH"),]
clone_type_IGH<-t(as.data.frame(unclass(table(data_chain_IGH$V_J_lenghCDR3_CloneId,factor(data_chain_IGH$sample))))) 

# data_chain_TRA<-data_merge_filter[which(data_merge_filter$chainType=="TRA"),]
# clone_type_TRA<-t(as.data.frame(unclass(table(data_chain_TRA$V_J_lenghCDR3_CloneId,factor(data_chain_TRA$sample))))) 
# data_chain_TRB<-data_merge_filter[which(data_merge_filter$chainType=="TRB"),]
# clone_type_TRB<-t(as.data.frame(unclass(table(data_chain_TRB$V_J_lenghCDR3_CloneId,factor(data_chain_TRB$sample))))) 
# data_chain_TRD<-data_merge_filter[which(data_merge_filter$chainType=="TRD"),]
# clone_type_TRD<-t(as.data.frame(unclass(table(data_chain_TRD$V_J_lenghCDR3_CloneId,factor(data_chain_TRD$sample))))) 
# data_chain_TRG<-data_merge_filter[which(data_merge_filter$chainType=="TRG"),]
# clone_type_TRG<-t(as.data.frame(unclass(table(data_chain_TRG$V_J_lenghCDR3_CloneId,factor(data_chain_TRG$sample))))) 

##Build present vs no present
clone_type_presence_IGK<-apply(clone_type_IGK,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGL<-apply(clone_type_IGL,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGH<-apply(clone_type_IGH,1,function(x) ifelse(x==0,0,1))
# clone_type_presence_TRA<-apply(clone_type_TRA,1,function(x) ifelse(x==0,0,1))
# clone_type_presence_TRB<-apply(clone_type_TRB,1,function(x) ifelse(x==0,0,1))
# clone_type_presence_TRD<-apply(clone_type_TRD,1,function(x) ifelse(x==0,0,1))
# clone_type_presence_TRG<-apply(clone_type_TRG,1,function(x) ifelse(x==0,0,1))

#
###Filter by clones share at least in 2% of the samples (311*0.05 = 15) or 2 samples
clone_type_filter_IGK<-clone_type_IGK[,which(rowSums(clone_type_presence_IGK)>dim(clone_type_IGK)[1]*0.02)] #
clone_type_filter_IGL<-clone_type_IGL[,which(rowSums(clone_type_presence_IGL)>dim(clone_type_IGL)[1]*0.02)] #
clone_type_filter_IGH<-clone_type_IGH[,which(rowSums(clone_type_presence_IGH)>dim(clone_type_IGH)[1]*0.02)] #
# clone_type_filter_TRA<-clone_type_TRA[,which(rowSums(clone_type_presence_TRA)>dim(clone_type_TRA)[1]*0.02)] #
# clone_type_filter_TRB<-clone_type_TRB[,which(rowSums(clone_type_presence_TRB)>dim(clone_type_TRB)[1]*0.02)] #
# clone_type_filter_TRD<-clone_type_TRD[,which(rowSums(clone_type_presence_TRD)>dim(clone_type_TRD)[1]*0.02)] #
# clone_type_filter_TRG<-clone_type_TRG[,which(rowSums(clone_type_presence_TRG)>dim(clone_type_TRG)[1]*0.02)] #

########################
### Applied Geometric Bayesian Multiplicative to substitute the zeros
########################
#clone_type_IGH.GBM <- cmultRepl(clone_type_filter_IGH, output="p-counts")
#clone_type_IGK.GBM <- cmultRepl(clone_type_filter_IGK, output="p-counts")
#clone_type_IGL.GBM <- cmultRepl(clone_type_filter_IGL, output="p-counts")

#####################
### Add an offset of 1 to substitute the zeros
clone_type_filter_IGH_zerosubs<-clone_type_filter_IGH+1
clone_type_filter_IGK_zerosubs<-clone_type_filter_IGK+1
clone_type_filter_IGL_zerosubs<-clone_type_filter_IGL+1
# clone_type_filter_TRA_zerosubs<-clone_type_filter_TRA+1
# clone_type_filter_TRB_zerosubs<-clone_type_filter_TRB+1
# clone_type_filter_TRD_zerosubs<-clone_type_filter_TRD+1
# clone_type_filter_TRG_zerosubs<-clone_type_filter_TRG+1


#########################
### CLR-transformation ##
########################
## IGK
z_clone_type_IGK<-log(clone_type_filter_IGK_zerosubs)
clrx_clone_type_IGK <- apply(z_clone_type_IGK, 2, function(x) x - rowMeans(z_clone_type_IGK))

id<-match(rownames(clone_type_filter_IGK_zerosubs),rownames(PAAD.repertoire.tumor.filter))

##Heatmap for all 
library(pheatmap)
library(RColorBrewer)
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F","#BEAED4")

res <- pheatmap(t(clrx_clone_type_IGK),scale="row",border_color=F,show_colnames = F,
                color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))

####################
#### Clustering ###
###################
mat.clust <- as.data.frame(cbind(clrx_clone_type_IGK, cluster = cutree(res$tree_col, k = 3)))

annotation_row = data.frame(cluster = factor(mat.clust$cluster))
colnames(annotation_row)<-c("cluster")
rownames(annotation_row)<-rownames(clrx_clone_type_IGK)

pheatmap(t(clrx_clone_type_IGK),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
        color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))

PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC<-mat.clust$cluster
save(data_merge,PAAD.repertoire.diversity,PAAD.repertoire.tumor.filter,xCell.data.PAAD,xCell.pvalue.PAAD,paad.subtype,clinical.drug,clinical.patient,clinical.radiation,clinical.new_tumor_event,clinical.folow_up,biospecimen.slide,annotation,
     file="Data/PAAD/PAAD_FullData.Rdata")


##################
### CLR-LASSO ###
################

#############
##1. Cox ###
############
z_clone_type_IGL<-log(clone_type_filter_IGL_zerosubs)
clrx_clone_type_IGL <- apply(z_clone_type_IGL, 2, function(x) x - rowMeans(z_clone_type_IGL))

PAAD.repertoire.tumor.filter.pos<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$OS.time>0),]
surv_object <- Surv(time = PAAD.repertoire.tumor.filter.pos$OS.time, event = PAAD.repertoire.tumor.filter.pos$OS)

id<-match(rownames(PAAD.repertoire.tumor.filter.pos),rownames(clrx_clone_type_IGL))
clrx_clone_type_IGL_OS<-clrx_clone_type_IGL[id,]


x<-model.matrix( ~ clrx_clone_type_IGL_OS+PAAD.repertoire.tumor.filter.pos$gender+ PAAD.repertoire.tumor.filter.pos$age_at_initial_pathologic_diagnosis
+PAAD.repertoire.tumor.filter.pos$pathologic_stage)
id<-match(rownames(clrx_clone_type_IGL_OS[as.numeric(rownames(x)),]),rownames(PAAD.repertoire.tumor.filter.pos))

set.seed(35)
IGL.surv_clrlasso.cv <- cv.glmnet(x = clrx_clone_type_IGL_OS, y = surv_object,family="cox",
                                  nfolds = 5, alpha=1)

IGL.surv_clrlasso <- glmnet(x = clrx_clone_type_IGL_OS, y = surv_object,family="cox", 
                            alpha=1,lambda=IGL.surv_clrlasso.cv$lambda.min)

clones_IGL_clr<-rownames(IGL.surv_clrlasso$beta)[which(as.numeric(IGL.surv_clrlasso$beta)!=0)]
clrx_clone_type_IGL[,match(clones_IGL_clr,colnames(clrx_clone_type_IGL_OS))]
clrx_clone_type_IGL_sign<-clrx_clone_type_IGL_OS[,match(clones_IGL_clr,colnames(clrx_clone_type_IGL_OS))]
xx<-data.frame(clrx_clone_type_IGL_sign)


RiskScore_compositional<-predict(IGL.surv_clrlasso,clrx_clone_type_IGL_OS)
names(RiskScore_compositional)<-rownames(clrx_clone_type_IGL_sign)
write.csv(RiskScore_compositional,"Data/PAAD/RiskScore_compositional.csv")
id<-match(rownames(PAAD.repertoire.tumor.filter),names(RiskScore_compositional))
PAAD.repertoire.tumor.filter$Risk_Score_compositional<-RiskScore_compositional[id]
save(data_merge,PAAD.repertoire.diversity,PAAD.repertoire.tumor.filter,xCell.data.PAAD,xCell.pvalue.PAAD,paad.subtype,clinical.drug,clinical.patient,clinical.radiation,clinical.new_tumor_event,clinical.folow_up,biospecimen.slide,annotation,
     file="Data/PAAD/PAAD_FullData.Rdata")

# tiff("Results/CompositionalAnalysis/IGK_forest_Sign_IGKV2-30_IGKJ4_33_42.tiff",res=300,h=2000,w=2200)
# res.cox <- coxph(surv_object~clrx_clone_type_IGK_sign)
# ggforest(res.cox,data=clrx_clone_type_IGL_sign)
# dev.off()

tiff("Results/CompositionalAnalysis/IGL_forest_ALL.tiff",res=300,h=2000,w=2200)
res.cox <- coxph(surv_object~IGLV1.36_IGLJ2_42_43+IGLV1.44_IGLJ1_39_102+IGLV10.54_IGLJ2_39_93+ IGLV2.11_IGLJ2_33_176+
                   IGLV2.14_IGLJ2_33_208+IGLV2.18_IGLJ3_36_93+IGLV2.8_IGLJ3_42_32+IGLV3.25_IGLJ3_36_113+IGLV3.27_IGLJ2_33_56+
                   IGLV5.37_IGLJ2_36_11+IGLV6.57_IGLJ3_39_130+IGLV7.43_IGLJ2_39_14+IGLV7.46_IGLJ2_36_106+IGLV7.46_IGLJ3_33_130,data=xx)
ggforest(res.cox,data=clrx_clone_type_IGL_sign)
dev.off()




tiff("Results/CompositionalAnalysis/IGL_forest_Sign.tiff",res=300,h=2000,w=2200)
res.cox <- coxph(surv_object~IGLV2.11_IGLJ2_33_176+
                   IGLV2.14_IGLJ2_33_208+
                   IGLV5.37_IGLJ2_36_11+IGLV6.57_IGLJ3_39_130,data=xx)

summary(res.cox)
ggforest(res.cox,data=clrx_clone_type_IGL_sign)
dev.off()

#############
### 2. Other variables
#PAAD.repertoire.tumor.filter$Mutated.Genes<-as.numeric(as.character(PAAD.repertoire.tumor.filter$Mutated.Genes))
#PAAD.repertoire.tumor.filter$KRAS.Mutated..1.or.0.<-factor(PAAD.repertoire.tumor.filter$KRAS.Mutated..1.or.0.)
#PAAD.repertoire.tumor.filter$ZNF521.mutated<-factor(PAAD.repertoire.tumor.filter$ZNF521)

####Tobacco
PAAD.repertoire.tumor.filter.nas<-PAAD.repertoire.tumor.filter[which(is.na(PAAD.repertoire.tumor.filter$number_pack_years_smoked)==F),]
id2<-match(rownames(PAAD.repertoire.tumor.filter.nas),rownames(clrx_clone_type_IGK))
IGK.test_clrlasso.cv <- cv.glmnet(x = clrx_clone_type_IGK[id2,], y = PAAD.repertoire.tumor.filter.nas$number_pack_years_smoked, 
                                  nfolds = 5, alpha=1)

IGK.test_clrlasso <- glmnet(x = clrx_clone_type_IGK[id2,], y = PAAD.repertoire.tumor.filter.nas$number_pack_years_smoked, 
                           alpha=1,lambda=IGK.test_clrlasso.cv$lambda.min)

clones_IGK_clr<-rownames(IGK.test_clrlasso$beta)[which(as.numeric(IGK.test_clrlasso$beta)!=0)]
clrx_clone_type_IGK[,match(clones_IGK_clr,colnames(clrx_clone_type_IGK))]
clrx_clone_type_IGK_sign<-clrx_clone_type_IGK[,match(clones_IGK_clr,colnames(clrx_clone_type_IGK))]

plot(clrx_clone_type_IGK_sign[id2,9],PAAD.repertoire.tumor.filter.nas$number_pack_years_smoked)

annotation_col = data.frame("Num.pack/years" = PAAD.repertoire.tumor.filter.nas$number_pack_years_smoked)
rownames(annotation_col)<-rownames(clrx_clone_type_IGK_sign[id2,])

tiff("Results/CompositionalAnalysis/HeatmapIGK_tobacco.tiff",res=300,h=2000,w=3000)
pheatmap(t(clrx_clone_type_IGK_sign[id2,]),scale="row",border_color=F,show_colnames = F,annotation_col = annotation_col,
         color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()


## IGL
z_clone_type_IGL<-log(clone_type_filter_IGL_zerosubs)
clrx_clone_type_IGL <- apply(z_clone_type_IGL, 2, function(x) x - rowMeans(z_clone_type_IGL))

##################
### CLR-LASSO ###
################
set.seed(35)
id<-match(rownames(clone_type_filter_IGL_zerosubs),rownames(PAAD.repertoire.tumor.filter))

####Tobacco
PAAD.repertoire.tumor.filter.nas<-PAAD.repertoire.tumor.filter[which(is.na(PAAD.repertoire.tumor.filter$number_pack_years_smoked)==F),]
id2<-match(rownames(PAAD.repertoire.tumor.filter.nas),rownames(clrx_clone_type_IGL))
IGL.test_clrlasso.cv <- cv.glmnet(x = clrx_clone_type_IGL[id2,], y = PAAD.repertoire.tumor.filter.nas$number_pack_years_smoked, 
                                  nfolds = 5, alpha=1)

IGL.test_clrlasso <- glmnet(x = clrx_clone_type_IGL[id2,], y = PAAD.repertoire.tumor.filter.nas$number_pack_years_smoked, 
                            alpha=1,lambda=IGL.test_clrlasso.cv$lambda.min)


clones_IGL_clr<-rownames(IGL.test_clrlasso$beta)[which(as.numeric(IGL.test_clrlasso$beta)!=0)]
clrx_clone_type_IGL[,match(clones_IGL_clr,colnames(clrx_clone_type_IGL))]
clrx_clone_type_IGL_sign<-clrx_clone_type_IGL[,match(clones_IGL_clr,colnames(clrx_clone_type_IGL))]

plot(clrx_clone_type_IGL_sign[id2,9],PAAD.repertoire.tumor.filter.nas$number_pack_years_smoked)

annotation_col = data.frame("Num.pack/years" = PAAD.repertoire.tumor.filter.nas$number_pack_years_smoked)
rownames(annotation_col)<-rownames(clrx_clone_type_IGL_sign[id2,])

tiff("Results/CompositionalAnalysis/HeatmapIGL_tobacco.tiff",res=300,h=2000,w=3000)
pheatmap(t(clrx_clone_type_IGL_sign[id2,]),scale="row",border_color=F,show_colnames = F,annotation_col = annotation_col,
         color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()


####
##Are this specific clonotypes better prognosis?
###

library(survival)
library(survminer)
library(survMisc)

surv_object <- Surv(time = PAAD.repertoire.tumor.filter$OS.time, event = PAAD.repertoire.tumor.filter$OS)
##tertiles los que fuman

res.cox <- coxph(Surv(time = OS.time, event = OS)~clrx_clone_type_IGK_sign[,4],data=PAAD.repertoire.tumor.filter)
summary(res.cox)

ggforest(res.cox)


