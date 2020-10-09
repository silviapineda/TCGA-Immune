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


################
### CLR-lasso ##
################
##
z_clone_type_IGK<-log(clone_type_filter_IGK_zerosubs)
clrx_clone_type_IGK <- apply(z_clone_type_IGK, 2, function(x) x - rowMeans(z_clone_type_IGK))

id<-match(rownames(clone_type_filter_IGK_zerosubs),rownames(PAAD.repertoire.tumor.filter))
set.seed(35)
IGK.test_clrlasso.cv <- cv.glmnet(x = clrx_clone_type_IGK, y = PAAD.repertoire.tumor.filter$Gender[id], 
                                  family = 'binomial', nfolds = 5, alpha=1)

IGK.test_clrlasso <- glmnet(x = clrx_clone_type_IGK, y = PAAD.repertoire.tumor.filter$Gender[id], 
                            family = 'binomial', alpha=1,lambda=IGK.test_clrlasso.cv$lambda.min)

clones_IGK_clr<-rownames(IGK.test_clrlasso$beta)[which(as.numeric(IGK.test_clrlasso$beta)!=0)]
clrx_clone_type_IGK[,match(clones_IGK_clr,colnames(clrx_clone_type_IGK))]
clrx_clone_type_IGK_sign<-clrx_clone_type_IGK[,match(clones_IGK_clr,colnames(clrx_clone_type_IGK))]

##Heatmap
library(pheatmap)
library(RColorBrewer)
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F","#BEAED4")

annotation_row = data.frame(PAAD.repertoire.tumor.filter$gender[id])
#ann_colors = list (outcome = c("FEMALE" = cols[1],"MALE" = cols[2]))
#colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(clrx_clone_type_IGK)

tiff("Results/CompositionalAnalysis/IGK_codaLasso.tiff",width = 2500, height = 1500, res = 300)
pheatmap(t(clrx_clone_type_IGK_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()


PAAD.repertoire.tumor.filter$ZNF521.mutated<-factor(PAAD.repertoire.tumor.filter$ZNF521)
PAAD.repertoire.tumor.filter$smoking<-factor(ifelse(PAAD.repertoire.tumor.filter$tobacco_smoking_history_master=="Current smoker (includes daily smokers and non-daily smokers or occasional smokers)","Current",
                                                    ifelse(PAAD.repertoire.tumor.filter$tobacco_smoking_history_master=="Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime)","Non-smoker",
                                                           ifelse(PAAD.repertoire.tumor.filter$tobacco_smoking_history_master=="Current reformed smoker for > 15 years (greater than 15 years)","Former",
                                                                  ifelse(PAAD.repertoire.tumor.filter$tobacco_smoking_history_master=="Current reformed smoker for ≤15 years (less than or equal to 15 years)","Former",
                                                                         ifelse(PAAD.repertoire.tumor.filter$tobacco_smoking_history_master=="Current reformed smoker, duration not specified","Former",NA))))))
PAAD.repertoire.tumor.filter$smoking2<-factor(ifelse(PAAD.repertoire.tumor.filter$smoking=="Current" |
                                                       PAAD.repertoire.tumor.filter$smoking=="Former" ,"Ever-Smoker",
                                                     ifelse(PAAD.repertoire.tumor.filter$smoking=="Non-smoker","Non-smoker",NA)))

set.seed(35)
lambdaRange_codalasso(X = clone_type_filter_IGK_zerosubs, y = PAAD.repertoire.tumor.filter$smoking2[id], lambdaSeq = seq(0.1, 1, 0.01))
codalasso_IGK <- coda_logistic_lasso(X = clone_type_filter_IGK_zerosubs, y = PAAD.repertoire.tumor.filter$gender[id], lambda = 0.1)
clones_IGK_coda<-codalasso_IGK$`name of selected variables`
coda_clone_type_IGK_sign<-clone_type_filter_IGK_zerosubs[,match(clones_IGK_coda,colnames(clone_type_filter_IGK_zerosubs))]



##Heatmap
library(pheatmap)
library(RColorBrewer)
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F","#BEAED4")

annotation_row = data.frame(PAAD.repertoire.tumor.filter$gender[id])
#ann_colors = list (outcome = c("FEMALE" = cols[1],"MALE" = cols[2]))
#colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(clrx_clone_type_IGK)

tiff("Results/CompositionalAnalysis/IGK_codaLasso.tiff",width = 2500, height = 1500, res = 300)
pheatmap(t(coda_clone_type_IGK_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()
