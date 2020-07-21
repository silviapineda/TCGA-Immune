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
library(compositions)
library(selbal)


setwd("~/TCGA-Immune/")

load("Data/PAAD_GTEx_ValTumor_ValNormal/PAAD_GTEx_ValTumor_ValNormal_FullData.Rdata")
##Filter for those that are pancreas or normal
########################################
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
###Filter by clones share at least in 5% of the samples (311*0.05 = 15) or 2 samples
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
clone_type_IGK.GBM <- cmultRepl(clone_type_filter_IGK, output="p-counts")
clone_type_IGL.GBM <- cmultRepl(clone_type_filter_IGL, output="p-counts")

#########################
### Selbal algorithm ###
########################
id<-match(rownames(clone_type_IGL.GBM),PAAD.GTEx.repertoire.diversity.tumor.normmal$sample)
CV.selbal.IGL <- selbal.cv(x = (clone_type_IGL.GBM), y = PAAD.GTEx.repertoire.diversity.tumor.normmal$outcome[id], n.fold = 5, n.iter = 10,
                        logit.acc = "AUC")
save(CV.selbal.IGL,file="Results/CompositionalAnalysis/Selbal_IGL.Rdata")

CV.selbal.IGH$accuracy.nvar
CV.selbal.IGH$var.barplot
grid.draw(CV.selbal.IGH$global.plot)
plot.tab(CV.selbal.IGH$cv.tab)
