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

table(annotation_gtex_pancreas$SEX)
summary(annotation_gtex_pancreas$AGE)
summary(annotation_gtex_pancreas$COD)


cols=brewer.pal(3,name = "Set3")

id<-match(Pancreas.repertoire.diversity$SUBJID,annotation_gtex_pancreas$SUBJID)
Pancreas.repertoire.diversity$sex<-annotation_gtex_pancreas[id,"SEX"]
Pancreas.repertoire.diversity$cod<-annotation_gtex_pancreas[id,"COD"]

##Sex
Ig_expr<-melt(Pancreas.repertoire.diversity[,c("SUBJID","IGH_expression","IGK_expression","IGL_expression","sex")])
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_GTEx_sex.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_expr, x = "sex", y = "value",facet.by = "variable",color = "sex",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=sex, y=value,color=sex), position = position_jitterdodge(dodge.width = 0.8)) +
  #scale_color_manual(values = c(cols[1],cols[2]), labels = c("Female","Male")) +
  stat_compare_means(
    comparisons =list(c("Female","Male")))
dev.off()

#COD
Ig_expr<-melt(Pancreas.repertoire.diversity[,c("SUBJID","IGH_expression","IGK_expression","IGL_expression","cod")])
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_GTEx_cod.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_expr, x = "cod", y = "value",facet.by = "variable",color = "cod",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=cod, y=value,color=cod), position = position_jitterdodge(dodge.width = 0.8)) +
  #scale_color_manual(values = c(cols[1],cols[2]), labels = c("Female","Male")) +
  stat_compare_means(label = "p.format")
dev.off()

