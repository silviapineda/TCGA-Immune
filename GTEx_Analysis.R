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
### DESCRIP: Validation analysis
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

load("Data/GTEx/Pancreas/GTEx_FullData.Rdata")


###############################
## Merge with Clinical data ###
###############################
Pancreas.repertoire.diversity.annotation<-merge(Pancreas.repertoire.diversity,annotation_gtex_pancreas,by="SUBJID")

##Function to run the association between clinical outcome and BCR/TCR
association.test.immuneRep<- function (PAAD.repertoire.tumor.clinical.patient,clinical.var){
  Ig_expr<-melt(PAAD.repertoire.tumor.clinical.patient[,c("sample","IGH_expression","IGK_expression","IGL_expression",clinical.var)])
  Ig_expr$value<-log10(Ig_expr$value)
  Ig_expr<-na.omit(Ig_expr)
  tiff(paste0("Results/boxplot_Ig_expression_TCGA_",clinical.var,".tiff"),res=300,h=2500,w=3500)
  if(length(table(Ig_expr[,clinical.var]))<10){  
    print(ggboxplot(Ig_expr, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
            rotate_x_text() +
            geom_point(aes(x=Ig_expr[,clinical.var], y=value,color=Ig_expr[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
            stat_compare_means(label = "p.format"))
  } else {  
    print(ggplot(Ig_expr,aes(x = as.numeric(Ig_expr[,clinical.var]), y = value)) + geom_point() + 
            facet_grid(.~Ig_expr$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()
  
  Ig_entropy<-melt(PAAD.repertoire.tumor.clinical.patient[,c("sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL",clinical.var)])
  Ig_entropy<-Ig_entropy[which(Ig_entropy$value!=0),]
  Ig_entropy<-na.omit(Ig_entropy)
  tiff(paste0("Results/boxplot_Ig_entropy_TCGA_",clinical.var,".tiff"),res=300,h=2500,w=3500)
  if(length(table(Ig_entropy[,clinical.var]))<10){  
    print(ggboxplot(Ig_entropy, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
            rotate_x_text() +
            geom_point(aes(x=Ig_entropy[,clinical.var], y=value,color=Ig_entropy[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
            stat_compare_means(label = "p.format"))
  } else {  
    print(ggplot(Ig_entropy,aes(x = as.numeric(Ig_entropy[,clinical.var]), y = value)) + geom_point() + 
            facet_grid(.~Ig_entropy$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()
  
  T_expr<-melt(PAAD.repertoire.tumor.clinical.patient[,c("sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression",clinical.var)])
  T_expr$value<-log10(T_expr$value)
  T_expr<-na.omit(T_expr)
  tiff(paste0("Results/boxplot_T_expr_TCGA_",clinical.var,".tiff"),res=300,h=2500,w=3500)
  if(length(table(T_expr[,clinical.var]))<10){  
    print(ggboxplot(T_expr, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
            rotate_x_text() +
            geom_point(aes(x=T_expr[,clinical.var], y=value,color=T_expr[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
            stat_compare_means(label = "p.format"))
  } else {  
    print(ggplot(T_expr,aes(x = as.numeric(T_expr[,clinical.var]), y = value)) + geom_point() + 
            facet_grid(.~T_expr$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()
  
  T_entropy<-melt(PAAD.repertoire.tumor.clinical.patient[,c("sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG",clinical.var)])
  T_entropy<-T_entropy[which(T_entropy$value!=0),]
  T_entropy<-na.omit(T_entropy)
  tiff(paste0("Results/boxplot_T_entropy_TCGA_",clinical.var,".tiff"),res=300,h=2500,w=3500)
  if(length(table(T_entropy[,clinical.var]))<10){  
    print(ggboxplot(T_entropy, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
            rotate_x_text() +
            geom_point(aes(x=T_entropy[,clinical.var], y=value,color=T_entropy[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
            stat_compare_means(label = "p.format"))
  } else {  
    print(ggplot(T_entropy,aes(x = as.numeric(T_entropy[,clinical.var]), y = value)) + geom_point() + 
            facet_grid(.~T_entropy$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()
  
  
  Alpha_Beta_ratio_expression<-melt(PAAD.repertoire.tumor.clinical.patient[,c("sample","Alpha_Beta_ratio_expression",clinical.var)])
  Alpha_Beta_ratio_expression<-na.omit(Alpha_Beta_ratio_expression)
  tiff(paste0("Results/boxplot_Alpha_Beta_ratio_expression_TCGA_",clinical.var,".tiff"),res=300,h=2500,w=3500)
  if(length(table(Alpha_Beta_ratio_expression[,clinical.var]))<10){  
    print(ggboxplot(Alpha_Beta_ratio_expression, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
            rotate_x_text() +
            geom_point(aes(x=Alpha_Beta_ratio_expression[,clinical.var], y=value,color=Alpha_Beta_ratio_expression[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
            stat_compare_means(label = "p.format"))
  } else {  
    print(ggplot(Alpha_Beta_ratio_expression,aes(x = as.numeric(Alpha_Beta_ratio_expression[,clinical.var]), y = value)) + geom_point() + 
            facet_grid(.~Alpha_Beta_ratio_expression$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson") + xlab(clinical.var))
  }
  dev.off()
  
  KappaLambda_ratio_expression<-melt(PAAD.repertoire.tumor.clinical.patient[,c("sample","KappaLambda_ratio_expression",clinical.var)])
  KappaLambda_ratio_expression<-na.omit(KappaLambda_ratio_expression)
  tiff(paste0("Results/boxplot_KappaLambda_ratio_expression_TCGA_",clinical.var,".tiff"),res=300,h=2500,w=3500)
  if(length(table(KappaLambda_ratio_expression[,clinical.var]))<10){  
    print(ggboxplot(KappaLambda_ratio_expression, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
            rotate_x_text() +
            geom_point(aes(x=KappaLambda_ratio_expression[,clinical.var], y=value,color=KappaLambda_ratio_expression[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
            stat_compare_means(label = "p.format"))
  } else {  
    print(ggplot(KappaLambda_ratio_expression,aes(x = as.numeric(KappaLambda_ratio_expression[,clinical.var]), y = value)) + geom_point() + 
            facet_grid(.~KappaLambda_ratio_expression$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()
}

##Gender
Pancreas.repertoire.diversity.annotation$SEX<-factor(Pancreas.repertoire.diversity.annotation$SEX)
association.test.immuneRep(Pancreas.repertoire.diversity.annotation,"SEX")

##Age
Pancreas.repertoire.diversity.annotation$AGE<-factor(Pancreas.repertoire.diversity.annotation$AGE)
association.test.immuneRep(Pancreas.repertoire.diversity.annotation,"AGE")

##Race
Pancreas.repertoire.diversity.annotation$RACE<-ifelse(Pancreas.repertoire.diversity.annotation$RACE=="Indian American",NA,
                                                             as.character(Pancreas.repertoire.diversity.annotation$RACE))
association.test.immuneRep(Pancreas.repertoire.diversity.annotation,"RACE")

##COD
Pancreas.repertoire.diversity.annotation$COD<-factor(Pancreas.repertoire.diversity.annotation$COD)
association.test.immuneRep(Pancreas.repertoire.diversity.annotation,"COD")

##diabete
Pancreas.repertoire.diversity.annotation$Diabetes_type2<-factor(Pancreas.repertoire.diversity.annotation$Diabetes_type2)
association.test.immuneRep(Pancreas.repertoire.diversity.annotation,"Diabetes_type2")


#####
#### Analysis of he validation pancreas
#####
load("Data/Pancreas_Validation/Pancreas_Validation_FullData.Rdata")
##sex
Pancreas.Validation.repertoire.diversity$tumor_stage<-factor(Pancreas.Validation.repertoire.diversity$tumor_stage)
association.test.immuneRep(Pancreas.Validation.repertoire.diversity,"tumor_stage")

