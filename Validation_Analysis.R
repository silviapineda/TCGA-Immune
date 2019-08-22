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

load("Data/Pancreas_Validation/Pancreas_Validation_FullData.Rdata")

Pancreas.Validation.repertoire.diversity.tumor<-Pancreas.Validation.repertoire.diversity[which(Pancreas.Validation.repertoire.diversity$tissue=="pancreas tumor"),]
table(Pancreas.Validation.repertoire.diversity.tumor$sex)
summary(Pancreas.Validation.repertoire.diversity.tumor$age)
table(Pancreas.Validation.repertoire.diversity.tumor$tumor_grade)
table(Pancreas.Validation.repertoire.diversity.tumor$tumor_stage)
Pancreas.Validation.repertoire.diversity.tumor$tumor_stage_group<-ifelse(Pancreas.Validation.repertoire.diversity.tumor$tumor_stage=="IA","I",
                                                                          ifelse(Pancreas.Validation.repertoire.diversity.tumor$tumor_stage=="IB","I",
                                                                                 ifelse(Pancreas.Validation.repertoire.diversity.tumor$tumor_stage=="IIA","II",
                                                                                        ifelse(Pancreas.Validation.repertoire.diversity.tumor$tumor_stage=="IIB","II",NA))))

cols=brewer.pal(3,name = "Set2")

##Sex
Ig_expr<-melt(Pancreas.Validation.repertoire.diversity.tumor[,c("sample","IGH_expression","IGK_expression","IGL_expression","sex")])
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_VAL_sex.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_expr, x = "sex", y = "value",facet.by = "variable",color = "sex",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=sex, y=value,color=sex), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[2]), labels = c("Female","Male")) +
  stat_compare_means(
    comparisons =list(c("Female","Male")))
dev.off()

##Grade
Ig_expr<-melt(Pancreas.Validation.repertoire.diversity.tumor[,c("sample","IGH_expression","IGK_expression","IGL_expression","tumor_stage_group")])
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_VAL_sex.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_expr, x = "tumor_stage_group", y = "value",facet.by = "variable",color = "tumor_stage_group",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=tumor_stage_group, y=value,color=tumor_stage_group), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[2]), labels = c("Female","Male")) +
  stat_compare_means(
    comparisons =list(c("Female","Male")))
dev.off()
