rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. Compositional Analysis
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
system('git clone https://github.com/UVic-omics/CoDA-Penalized-Regression')
# cran.packages <- c('knitr', 'glmnet', 'ggplot2', 'gridExtra',
#                    'UpSetR', 'ggforce')
# install.packages(cran.packages)
devtools::install_github(repo = 'UVic-omics/selbal')

library(knitr) # rbookdown, kable
library(glmnet) # glmnet
library(selbal) # selbal
library(ggplot2) # draw selbal
library(gridExtra) # grid.arrange
library(UpSetR) # upset
library(ggforce) # selbal-like plot
library(grid) # grid.draw
library(pheatmap)
library(RColorBrewer)
# source coda-lasso functions
source(file = './CoDA-Penalized-Regression/R/functions_coda_penalized_regression.R')

setwd("~/TCGA-Immune/")

load("Data/PAAD_GTEx_ValTumor_ValNormal/PAAD_GTEx_ValTumor_ValNormal_FullData.Rdata")
##Filter for those that are pancreas or normal
########################################
### 5e9e81e2-0e8b-4aca-aced-6ce451fa3262
grep("5e9e81e2-0e8b-4aca-aced-6ce451fa3262",rownames(PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity))
PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity<-PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity[-59,]

PAAD.GTEx.repertoire.diversity.tumor.normmal<-PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity[which(PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome=="normal_pancreas (GTEx)"|
                                                                                                        PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome== "PDAC (TCGA)"),]
PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome<-factor(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome)

id<-match(data_merge$sample,rownames(PAAD.GTEx.repertoire.diversity.tumor.normmal))
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

data_chain_TRA<-data_merge_filter[which(data_merge_filter$chainType=="TRA"),]
clone_type_TRA<-t(as.data.frame(unclass(table(data_chain_TRA$V_J_lenghCDR3_CloneId,factor(data_chain_TRA$sample))))) 
data_chain_TRB<-data_merge_filter[which(data_merge_filter$chainType=="TRB"),]
clone_type_TRB<-t(as.data.frame(unclass(table(data_chain_TRB$V_J_lenghCDR3_CloneId,factor(data_chain_TRB$sample))))) 
data_chain_TRD<-data_merge_filter[which(data_merge_filter$chainType=="TRD"),]
clone_type_TRD<-t(as.data.frame(unclass(table(data_chain_TRD$V_J_lenghCDR3_CloneId,factor(data_chain_TRD$sample))))) 
data_chain_TRG<-data_merge_filter[which(data_merge_filter$chainType=="TRG"),]
clone_type_TRG<-t(as.data.frame(unclass(table(data_chain_TRG$V_J_lenghCDR3_CloneId,factor(data_chain_TRG$sample))))) 

##Build present vs no present
clone_type_presence_IGK<-apply(clone_type_IGK,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGL<-apply(clone_type_IGL,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGH<-apply(clone_type_IGH,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRA<-apply(clone_type_TRA,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRB<-apply(clone_type_TRB,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRD<-apply(clone_type_TRD,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRG<-apply(clone_type_TRG,1,function(x) ifelse(x==0,0,1))

#
###Filter by clones share at least in 2% of the samples (311*0.05 = 15) or 2 samples
clone_type_filter_IGK<-clone_type_IGK[,which(rowSums(clone_type_presence_IGK)>dim(clone_type_IGK)[1]*0.02)] #
clone_type_filter_IGL<-clone_type_IGL[,which(rowSums(clone_type_presence_IGL)>dim(clone_type_IGL)[1]*0.02)] #
clone_type_filter_IGH<-clone_type_IGH[,which(rowSums(clone_type_presence_IGH)>dim(clone_type_IGH)[1]*0.02)] #
clone_type_filter_TRA<-clone_type_TRA[,which(rowSums(clone_type_presence_TRA)>dim(clone_type_TRA)[1]*0.02)] #
clone_type_filter_TRB<-clone_type_TRB[,which(rowSums(clone_type_presence_TRB)>dim(clone_type_TRB)[1]*0.02)] #
clone_type_filter_TRD<-clone_type_TRD[,which(rowSums(clone_type_presence_TRD)>dim(clone_type_TRD)[1]*0.02)] #
clone_type_filter_TRG<-clone_type_TRG[,which(rowSums(clone_type_presence_TRG)>dim(clone_type_TRG)[1]*0.02)] #

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
clone_type_filter_TRA_zerosubs<-clone_type_filter_TRA+1
clone_type_filter_TRB_zerosubs<-clone_type_filter_TRB+1
clone_type_filter_TRD_zerosubs<-clone_type_filter_TRD+1
clone_type_filter_TRG_zerosubs<-clone_type_filter_TRG+1


#########################
### Selbal algorithm ###
########################
#IGH
set.seed(35)
id<-match(rownames(clone_type_filter_IGH),PAAD.GTEx.repertoire.diversity.tumor.normmal$sample)
CV.selbal.IGH <- selbal.cv(x = (clone_type_filter_IGH), y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], n.fold = 5, n.iter = 2,
                           logit.acc = "AUC",zero.rep = "one")
save(CV.selbal.IGH,file="Results/CompositionalAnalysis/Selbal_IGH_filter_1.Rdata")

CV.selbal.IGH$accuracy.nvar
CV.selbal.IGH$var.barplot
grid.draw(CV.selbal.IGH$global.plot)
plot.tab(CV.selbal.IGH$cv.tab)

#IGK
id<-match(rownames(clone_type_filter_IGK),PAAD.GTEx.repertoire.diversity.tumor.normmal$sample)
CV.selbal.IGK <- selbal.cv(x = (clone_type_filter_IGK), y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], n.fold = 5, n.iter = 2,
                           logit.acc = "AUC",zero.rep = "one")
save(CV.selbal.IGK,file="Results/CompositionalAnalysis/Selbal_IGK_filter_1.Rdata")

load("Results/CompositionalAnalysis/Selbal_IGK_filter_1.Rdata")
CV.selbal.IGK$accuracy.nvar
CV.selbal.IGH$var.barplot
grid.draw(CV.selbal.IGH$global.plot)
plot.tab(CV.selbal.IGH$cv.tab)

#IGL
id<-match(rownames(clone_type_filter_IGL),PAAD.GTEx.repertoire.diversity.tumor.normmal$sample)
CV.selbal.IGL <- selbal.cv(x = (clone_type_filter_IGL), y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], n.fold = 5, n.iter = 10,
                           logit.acc = "AUC",zero.rep = "one")
save(CV.selbal.IGL,file="Selbal_IGL_filter_1.Rdata")




################
### clr-lasso ##
################
#IGH
z_clone_type_IGH<-log(clone_type_filter_IGH_zerosubs)
clrx_clone_type_IGH <- apply(z_clone_type_IGH, 2, function(x) x - rowMeans(z_clone_type_IGH))

id<-match(rownames(clrx_clone_type_IGH),PAAD.GTEx.repertoire.diversity.tumor.normmal$sample)
set.seed(35)
IGH.test_clrlasso.cv <- cv.glmnet(x = clrx_clone_type_IGH, y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], 
                                  family = 'binomial', nfolds = 5, alpha=1)

IGH.test_clrlasso <- glmnet(x = clrx_clone_type_IGH, y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], 
                            family = 'binomial', alpha=1,lambda=IGH.test_clrlasso.cv$lambda.min)
##2 selected
clones_IGH_clr<-rownames(IGH.test_clrlasso$beta)[which(as.numeric(IGH.test_clrlasso$beta)!=0)]
clrx_clone_type_IGH[,match(clones_IGH_clr,colnames(clrx_clone_type_IGH))]
clrx_clone_type_IGH_sign<-clrx_clone_type_IGH[,match(clones_IGH_clr,colnames(clrx_clone_type_IGH))]
#boxplot(clrx_clone_type_IGK[,match(clones_IGK,colnames(clrx_clone_type_IGK))[1]]~PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])

##Plot results
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F","#BEAED4")

annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
ann_colors = list (outcome = c("normal_pancreas (GTEx)" = cols[1],"PDAC (TCGA)" = cols[2]))
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(clrx_clone_type_IGH)

tiff("Results/CompositionalAnalysis/IGH_clrLasso.tiff",width = 2500, height = 1500, res = 300)
pheatmap(t(clrx_clone_type_IGH_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()

#IGK
z_clone_type_IGK<-log(clone_type_filter_IGK_zerosubs)
clrx_clone_type_IGK <- apply(z_clone_type_IGK, 2, function(x) x - rowMeans(z_clone_type_IGK))

id<-match(rownames(clrx_clone_type_IGK),PAAD.GTEx.repertoire.diversity.tumor.normmal$sample)
set.seed(35)
IGK.test_clrlasso.cv <- cv.glmnet(x = clrx_clone_type_IGK, y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], 
                            family = 'binomial', nfolds = 5, alpha=1)

IGK.test_clrlasso <- glmnet(x = clrx_clone_type_IGK, y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], 
                              family = 'binomial', alpha=1,lambda=IGK.test_clrlasso.cv$lambda.min)

clones_IGK_clr<-rownames(IGK.test_clrlasso$beta)[which(as.numeric(IGK.test_clrlasso$beta)!=0)]
clrx_clone_type_IGK[,match(clones_IGK_clr,colnames(clrx_clone_type_IGK))]
clrx_clone_type_IGK_sign<-clrx_clone_type_IGK[,match(clones_IGK_clr,colnames(clrx_clone_type_IGK))]
#boxplot(clrx_clone_type_IGK[,match(clones_IGK,colnames(clrx_clone_type_IGK))[1]]~PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])

##Plot results
brewer.pal(4,name = "Accent")
cols=c( "#7FC97F","#BEAED4")

annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
ann_colors = list (outcome = c("normal_pancreas (GTEx)" = cols[1],"PDAC (TCGA)" = cols[2]))
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(clrx_clone_type_IGK)

tiff("Results/CompositionalAnalysis/IGK_clrLasso.tiff",width = 2500, height = 1500, res = 300)
pheatmap(t(clrx_clone_type_IGK_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()

#IGL
z_clone_type_IGL<-log(clone_type_filter_IGL_zerosubs)
clrx_clone_type_IGL <- apply(z_clone_type_IGL, 2, function(x) x - rowMeans(z_clone_type_IGL))

id<-match(rownames(clrx_clone_type_IGL),PAAD.GTEx.repertoire.diversity.tumor.normmal$sample)
set.seed(35)
IGL.test_clrlasso.cv <- cv.glmnet(x = clrx_clone_type_IGL, y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], 
                                  family = 'binomial', nfolds = 5, alpha=1)

IGL.test_clrlasso <- glmnet(x = clrx_clone_type_IGL, y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], 
                            family = 'binomial', alpha=1,lambda=IGL.test_clrlasso.cv$lambda.min)

clones_IGL_clr<-rownames(IGL.test_clrlasso$beta)[which(as.numeric(IGL.test_clrlasso$beta)!=0)]
clrx_clone_type_IGL[,match(clones_IGL_clr,colnames(clrx_clone_type_IGL))]
clrx_clone_type_IGL_sign<-clrx_clone_type_IGL[,match(clones_IGL_clr,colnames(clrx_clone_type_IGL))]

##Heatmap
annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
ann_colors = list (outcome = c("normal_pancreas (GTEx)" = cols[1],"PDAC (TCGA)" = cols[2]))
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(clrx_clone_type_IGL)

tiff("Results/CompositionalAnalysis/IGL_clrLasso.tiff",width = 2500, height = 1500, res = 300)
pheatmap(t(clrx_clone_type_IGL_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()


################
### CODA-lasso ##
################
#IGH
set.seed(35)
id<-match(rownames(clone_type_filter_IGH_zerosubs),PAAD.GTEx.repertoire.diversity.tumor.normmal$sample)
lambdaRange_codalasso(X = clone_type_filter_IGH_zerosubs, y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], lambdaSeq = seq(0.1, 1, 0.01))
codalasso_IGH <- coda_logistic_lasso(X = clone_type_filter_IGH_zerosubs, y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], lambda = 0.13)
clones_IGH_coda<-codalasso_IGH$`name of selected variables`
coda_clone_type_IGH_sign<-clrx_clone_type_IGH[,match(clones_IGK_coda,colnames(clrx_clone_type_IGH))]
##Non results found

#IGK
set.seed(35)
id<-match(rownames(clone_type_filter_IGK_zerosubs),PAAD.GTEx.repertoire.diversity.tumor.normmal$sample)
lambdaRange_codalasso(X = clone_type_filter_IGK_zerosubs, y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], lambdaSeq = seq(0.1, 1, 0.01))
codalasso_IGK <- coda_logistic_lasso(X = clone_type_filter_IGK_zerosubs, y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], lambda = 0.13)
clones_IGK_coda<-codalasso_IGK$`name of selected variables`
coda_clone_type_IGK_sign<-clrx_clone_type_IGK[,match(clones_IGK_coda,colnames(clrx_clone_type_IGK))]

##Heatmap
annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
ann_colors = list (outcome = c("normal_pancreas (GTEx)" = cols[1],"PDAC (TCGA)" = cols[2]))
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(clrx_clone_type_IGK)

tiff("Results/CompositionalAnalysis/IGK_codaLasso.tiff",width = 2500, height = 1500, res = 300)
pheatmap(t(coda_clone_type_IGK_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()

#IGL
set.seed(35)
id<-match(rownames(clone_type_filter_IGL_zerosubs),PAAD.GTEx.repertoire.diversity.tumor.normmal$sample)
lambdaRange_codalasso(X = clone_type_filter_IGL_zerosubs, y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], lambdaSeq = seq(0.1, 1, 0.01))
codalasso_IGL <- coda_logistic_lasso(X = clone_type_filter_IGL_zerosubs, y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], lambda = 0.33)
clones_IGL_coda<-codalasso_IGL$`name of selected variables`
coda_clone_type_IGL_sign<-clrx_clone_type_IGL[,match(clones_IGL_coda,colnames(clrx_clone_type_IGL))]

##Heatmap
annotation_row = data.frame(PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id])
ann_colors = list (outcome = c("normal_pancreas (GTEx)" = cols[1],"PDAC (TCGA)" = cols[2]))
colnames(annotation_row)<-"outcome"
rownames(annotation_row)<-rownames(clrx_clone_type_IGL)

tiff("Results/CompositionalAnalysis/IGL_codaLasso.tiff",width = 2500, height = 1500, res = 300)
pheatmap(t(coda_clone_type_IGL_sign),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         annotation_colors = ann_colors,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()

intersect(colnames(clrx_clone_type_IGK_sign),colnames(coda_clone_type_IGK_sign))
intersect(colnames(clrx_clone_type_IGL_sign),colnames(coda_clone_type_IGL_sign))

###############################
### Are they validated?? ######
###############################
Repertoire.Diversity<-PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity[which(PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome=="normal_pancreas (GTEx)"|
                                                                                PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome== "PDAC (TCGA)" |
                                                                                PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome== "PDAC (Val)" |  
                                                                                PAAD.GTEx.ValTumor.ValNormal.repertoire.diversity$outcome== "normal_pancreas (Val)"),]
Repertoire.Diversity$outcome<-factor(Repertoire.Diversity$outcome)
# Rename all levels, by name
levels(Repertoire.Diversity$outcome) <- list("GTEX-Normal"="normal_pancreas (GTEx)", "TCGA-PDAC"= "PDAC (TCGA)", "Validation-PDAC"="PDAC (Val)","Validation-Normal"="normal_pancreas (Val)")
id<-match(data_merge$sample,rownames(Repertoire.Diversity))
data_merge_filter<-data_merge[which(is.na(id)==F),]

#chain=c("IGHV","IGKV","IGLV")
############
## 1. Build the matrix with the clones by samples
data_chain_IGK<-data_merge_filter[which(data_merge_filter$chainType=="IGK"),]
clone_type_IGK<-t(as.data.frame(unclass(table(data_chain_IGK$V_J_lenghCDR3_CloneId,factor(data_chain_IGK$sample))))) 
data_chain_IGL<-data_merge_filter[which(data_merge_filter$chainType=="IGL"),]
clone_type_IGL<-t(as.data.frame(unclass(table(data_chain_IGL$V_J_lenghCDR3_CloneId,factor(data_chain_IGL$sample))))) 
data_chain_IGH<-data_merge_filter[which(data_merge_filter$chainType=="IGH"),]
clone_type_IGH<-t(as.data.frame(unclass(table(data_chain_IGH$V_J_lenghCDR3_CloneId,factor(data_chain_IGH$sample))))) 

data_chain_TRA<-data_merge_filter[which(data_merge_filter$chainType=="TRA"),]
clone_type_TRA<-t(as.data.frame(unclass(table(data_chain_TRA$V_J_lenghCDR3_CloneId,factor(data_chain_TRA$sample))))) 
data_chain_TRB<-data_merge_filter[which(data_merge_filter$chainType=="TRB"),]
clone_type_TRB<-t(as.data.frame(unclass(table(data_chain_TRB$V_J_lenghCDR3_CloneId,factor(data_chain_TRB$sample))))) 
data_chain_TRD<-data_merge_filter[which(data_merge_filter$chainType=="TRD"),]
clone_type_TRD<-t(as.data.frame(unclass(table(data_chain_TRD$V_J_lenghCDR3_CloneId,factor(data_chain_TRD$sample))))) 
data_chain_TRG<-data_merge_filter[which(data_merge_filter$chainType=="TRG"),]
clone_type_TRG<-t(as.data.frame(unclass(table(data_chain_TRG$V_J_lenghCDR3_CloneId,factor(data_chain_TRG$sample))))) 

##Build present vs no present
clone_type_presence_IGK<-apply(clone_type_IGK,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGL<-apply(clone_type_IGL,1,function(x) ifelse(x==0,0,1))
clone_type_presence_IGH<-apply(clone_type_IGH,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRA<-apply(clone_type_TRA,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRB<-apply(clone_type_TRB,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRD<-apply(clone_type_TRD,1,function(x) ifelse(x==0,0,1))
clone_type_presence_TRG<-apply(clone_type_TRG,1,function(x) ifelse(x==0,0,1))

###Filter by clones share at least by 2 samples
clone_type_filter_IGK<-clone_type_presence_IGK[which(rowSums(clone_type_presence_IGK)>2),] #
clone_type_filter_IGL<-clone_type_presence_IGL[which(rowSums(clone_type_presence_IGL)>2),] #
clone_type_filter_IGH<-clone_type_presence_IGH[which(rowSums(clone_type_presence_IGH)>2),] #
clone_type_filter_TRA<-clone_type_presence_TRA[which(rowSums(clone_type_presence_TRA)>2),] #
clone_type_filter_TRB<-clone_type_presence_TRB[which(rowSums(clone_type_presence_TRB)>2),] #
clone_type_filter_TRD<-clone_type_presence_TRD[which(rowSums(clone_type_presence_TRD)>2),] #
clone_type_filter_TRG<-clone_type_presence_TRG[which(rowSums(clone_type_presence_TRG)>2),] #


#####################
### Add an offset of 1 to substitute the zeros
clone_type_filter_IGH_zerosubs<-clone_type_filter_IGH+1
clone_type_filter_IGK_zerosubs<-clone_type_filter_IGK+1
clone_type_filter_IGL_zerosubs<-clone_type_filter_IGL+1

###IGH
z_clone_type_IGH<-log(clone_type_filter_IGH_zerosubs)
clrx_clone_type_IGH <- apply(z_clone_type_IGH, 2, function(x) x - rowMeans(z_clone_type_IGH))

id.igh<-match(colnames(clrx_clone_type_IGH_sign),rownames(clrx_clone_type_IGH))
clrx_clone_type_IGH_val<-clrx_clone_type_IGH[id.igh,]


#IGH
brewer.pal(4,name = "Accent")
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")
id<-match(colnames(clrx_clone_type_IGH_val),Repertoire.Diversity$sample)
annotation_col = data.frame(Repertoire.Diversity$outcome[id])
table(Repertoire.Diversity$outcome)
ann_colors = list (outcome = c("GTEX-Normal" = cols[1],
                               "TCGA-PDAC" = cols[2],
                               "Validation-Normal"= cols[3],
                               "Validation-PDAC" = cols[4]))
colnames(annotation_col)<-"outcome"
rownames(annotation_col)<-colnames(clrx_clone_type_IGH_val)

cols2 = brewer.pal(12,name = "Paired")

tiff("Results/CompositionalAnalysis/IGH_clrLasso_validation.tiff",width = 2500, height = 1500, res = 300)
pheatmap(clrx_clone_type_IGH_val,scale="row",border_color=F,show_colnames = F, annotation_col = annotation_col,
         annotation_colors = ann_colors,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()


###IGK
z_clone_type_IGK<-log(clone_type_filter_IGK_zerosubs)
clrx_clone_type_IGK <- apply(z_clone_type_IGK, 2, function(x) x - rowMeans(z_clone_type_IGK))

id.IGK<-match(colnames(clrx_clone_type_IGK_sign),rownames(clrx_clone_type_IGK))
clrx_clone_type_IGK_val<-clrx_clone_type_IGK[id.IGK,]

id.IGK<-match(colnames(coda_clone_type_IGK_sign),rownames(clrx_clone_type_IGK))
clrx_clone_type_IGK_val<-clrx_clone_type_IGK[id.IGK,]

#IGK
brewer.pal(4,name = "Accent")
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")
id<-match(colnames(clrx_clone_type_IGK_val),Repertoire.Diversity$sample)
annotation_col = data.frame(Repertoire.Diversity$outcome[id])
table(Repertoire.Diversity$outcome)
ann_colors = list (outcome = c("GTEX-Normal" = cols[1],
                               "TCGA-PDAC" = cols[2],
                               "Validation-Normal"= cols[3],
                               "Validation-PDAC" = cols[4]))
colnames(annotation_col)<-"outcome"
rownames(annotation_col)<-colnames(coda_clone_type_IGK_val)

cols2 = brewer.pal(12,name = "Paired")

tiff("Results/CompositionalAnalysis/IGK_codaLasso_validation.tiff",width = 2500, height = 1500, res = 300)
pheatmap(clrx_clone_type_IGK_val,scale="row",border_color=F,show_colnames = F, annotation_col = annotation_col,
         annotation_colors = ann_colors,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()


###IGL
z_clone_type_IGL<-log(clone_type_filter_IGL_zerosubs)
clrx_clone_type_IGL <- apply(z_clone_type_IGL, 2, function(x) x - rowMeans(z_clone_type_IGL))

id.IGL<-match(colnames(clrx_clone_type_IGL_sign),rownames(clrx_clone_type_IGL))
clrx_clone_type_IGL_val<-clrx_clone_type_IGL[id.IGL,]

id.IGL<-match(colnames(coda_clone_type_IGL_sign),rownames(clrx_clone_type_IGL))
clrx_clone_type_IGL_val<-clrx_clone_type_IGL[id.IGL,]

#IGL
brewer.pal(4,name = "Accent")
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")
id<-match(colnames(clrx_clone_type_IGL_val),Repertoire.Diversity$sample)
annotation_col = data.frame(Repertoire.Diversity$outcome[id])
table(Repertoire.Diversity$outcome)
ann_colors = list (outcome = c("GTEX-Normal" = cols[1],
                               "TCGA-PDAC" = cols[2],
                               "Validation-Normal"= cols[3],
                               "Validation-PDAC" = cols[4]))
colnames(annotation_col)<-"outcome"
rownames(annotation_col)<-colnames(clrx_clone_type_IGL_val)

cols2 = brewer.pal(12,name = "Paired")

tiff("Results/CompositionalAnalysis/IGL_codaLasso_validation.tiff",width = 2500, height = 1500, res = 300)
pheatmap(clrx_clone_type_IGL_val,scale="row",border_color=F,show_colnames = F, annotation_col = annotation_col,
         annotation_colors = ann_colors,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()

################
write.csv(clrx_clone_type_IGH_sign,"Results/CompositionalAnalysis/IGH_clrLasso_clonotypes_sign.csv")
write.csv(clrx_clone_type_IGK_sign,"Results/CompositionalAnalysis/IGK_clrLasso_clonotypes_sign.csv")
write.csv(clrx_clone_type_IGL_sign,"Results/CompositionalAnalysis/IGL_clrLasso_clonotypes_sign.csv")

write.csv(coda_clone_type_IGK_sign,"Results/CompositionalAnalysis/IGK_codaLasso_clonotypes_sign.csv")
write.csv(coda_clone_type_IGL_sign,"Results/CompositionalAnalysis/IGL_codaLasso_clonotypes_sign.csv")



########################
###Study the results ###
########################
load("Data/PAAD/PAAD_FullData.Rdata")

library(survival)
library(survminer)
library(survMisc)

##OS
id<-match(rownames(PAAD.repertoire.tumor.filter),rownames(clrx_clone_type_IGL_sign))
surv_object <- Surv(time = PAAD.repertoire.tumor.filter$OS.time, event = PAAD.repertoire.tumor.filter$OS)
res.cox <- coxph(Surv(time = OS.time, event = OS)~clrx_clone_type_IGL_sign[id,7],data=PAAD.repertoire.tumor.filter)
summary(res.cox)

ggforest(res.cox)

##Categorical
KL_mean<-mean(clrx_clone_type_IGL_sign[id,3],na.rm=T)
PAAD.repertoire.tumor.filter$KL_ratio_2cat<-ifelse(as.numeric(clrx_clone_type_IGL_sign[id,3])<=KL_mean,1,2)
fit1 <- survfit(surv_object ~ PAAD.repertoire.tumor.filter$KL_ratio_2cat)
fit1
tiff("Results/CompositionalAnalysis/IGLV1-40_IGLJ3_14_334_KM.tiff",res=300,h=2000,w=2000)
ggsurvplot(fit1, data = PAAD.repertoire.tumor.filter)
dev.off()
comp(ten(fit1))$tests$lrTests

