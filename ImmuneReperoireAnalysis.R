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
### DESCRIP: Summary plots comparing normal&pseudonormal with tumoral
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
library(reshape)
setwd("~/TCGA-Immune/")


##################
####Descriptive analysis to see if there are differences by tumor and adjacent_normal
##################

#################
### 1. TCGA #####
#################
load("Data/PAAD/PAAD_FullData.Rdata")
brewer.pal(3,name = "Accent")
brewer.pal(3,name = "Pastel1")
cols=c( "#7FC97F", "#FBB4AE","#BEAED4", "#FDC086")

#Barplot
tiff("Results/barplot_Treads.tiff",res=300,h=2500,w=4000)
barplot(PAAD.repertoire.diversity$T_Reads,col=cols[PAAD.repertoire.diversity$Tumor_type_3categ],main="Number of T-Reads",xlab = "Samples", ylab = "Reads",las=2)
abline(h=100)
legend("topright", legend=levels(PAAD.repertoire.diversity$Tumor_type_3categ),col=cols,pch=15, cex=0.8)
dev.off()

tiff("Results/barplot_Igreads.tiff",res=300,h=2500,w=4000)
barplot(PAAD.repertoire.diversity$Ig_Reads_filter,col=cols[PAAD.repertoire.diversity$Tumor_type_4categ],main="Number of Ig-Reads",xlab = "Samples", ylab = "Reads",las=2)
abline(h=100)
legend("topright", legend=levels(PAAD.repertoire.diversity$Tumor_type_4categ),col=cols,pch=15, cex=0.8)
dev.off()


tiff("Results/Corr_plot_reads.tiff",res=300,h=2500,w=4000)
plot(PAAD.repertoire.diversity$Total_Reads,PAAD.repertoire.diversity$totalSeqReads,col=cols[PAAD.repertoire.diversity$Tumor_type_3categ],pch=19,
     xlab = "Total B- and T- reads aligned", ylab = "Total sequencing reads",
     main = paste0("rho = ",round(cor(PAAD.repertoire.diversity$Total_Reads,PAAD.repertoire.diversity$totalSeqReads),2)))
legend("bottomright", legend=levels(PAAD.repertoire.diversity$Tumor_type_3categ),col=cols,pch=15, cex=0.8)
dev.off()

summary(PAAD.repertoire.diversity$IG_Reads[which(PAAD.repertoire.diversity$Tumor_type_3categ=="Tumor_pancreas")])
summary(PAAD.repertoire.diversity$T_Reads[which(PAAD.repertoire.diversity$Tumor_type_3categ=="normal_pancreas")])
summary(PAAD.repertoire.diversity$totalSeqReads)

##Filter by number of reads>100 and clones <100
PAAD.repertoire.diversity_treads<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$T_Reads_filter>100),]
PAAD.repertoire.diversity_Igreads<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Ig_Reads_filter>110),]

#Recon
PAAD.repertoire.diversity$T_clones_recon<-PAAD.repertoire.diversity$clones_recon_TRA+PAAD.repertoire.diversity$clones_recon_TRB+
                                    PAAD.repertoire.diversity$clones_recon_TRD+PAAD.repertoire.diversity$clones_recon_TRG
PAAD.repertoire.diversity$Ig_clones_recon<-PAAD.repertoire.diversity$clones_recon_IGH+PAAD.repertoire.diversity$clones_recon_IGK+
                                     PAAD.repertoire.diversity$clones_recon_IGL
PAAD.repertoire.diversity_tclones_recon<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$T_clones_recon>100),]
PAAD.repertoire.diversity_Igclones_recon<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Ig_clones_recon>100),]

#
PAAD.repertoire.diversity$T_clones<-PAAD.repertoire.diversity$clones_TRA+PAAD.repertoire.diversity$clones_TRB+
  PAAD.repertoire.diversity$clones_TRD+PAAD.repertoire.diversity$clones_TRG
PAAD.repertoire.diversity$Ig_clones<-PAAD.repertoire.diversity$clones_IGH+PAAD.repertoire.diversity$clones_IGK+
  PAAD.repertoire.diversity$clones_IGL
PAAD.repertoire.diversity_tclones<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$T_clones>100),]
PAAD.repertoire.diversity_Igclones<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Ig_clones>100),]

#These are outliers that disappear when filter by cdr3
#id<-match(c("TCGA-F2-6879-01A","TCGA-FB-AAQ3-01A","TCGA-F2-6880-01A","TCGA-HZ-A8P1-01A"),PAAD.repertoire.diversity$TCGA_sample)
####Summary plots
Ig_expr<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","IGH_expression_filter","IGK_expression_filter","IGL_expression_filter","Tumor_type_4categ")])
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_filter_PAAD.tiff",res=300,h=2500,w=3000)
ggboxplot(Ig_expr, x = "Tumor_type_4categ", y = "value",facet.by = "variable",color = "Tumor_type_4categ",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=value, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))

dev.off()

TR_expr<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","TRA_expression_filter","TRB_expression_filter","TRD_expression_filter","TRG_expression_filter","Tumor_type_4categ")])
TR_expr<-TR_expr[which(TR_expr$value!=0),]
TR_expr$value<-log10(TR_expr$value)
tiff("Results/boxplot_TR_expression_filter_PAAD.tiff",res=300,h=2500,w=3500)
ggboxplot(TR_expr, x = "Tumor_type_4categ", y = "value",facet.by = "variable",color = "Tumor_type_4categ",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=value, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))
dev.off()

Ig_entropy<-melt(PAAD.repertoire.diversity_Igclones[,c("TCGA_sample","entropy_IGH","entropy_IGK","entropy_IGL","Tumor_type_4categ")])
tiff("Results/boxplot_Ig_entropy_PAAD.tiff",res=300,h=2500,w=3000)
ggboxplot(Ig_entropy, x = "Tumor_type_4categ", y = "value",facet.by = "variable",color = "Tumor_type_4categ",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=value, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))
dev.off()

Ig_entropy<-melt(PAAD.repertoire.diversity_Igclones[,c("TCGA_sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL","Tumor_type_4categ")])
tiff("Results/boxplot_Ig_entropy_recon_PAAD.tiff",res=300,h=2500,w=3000)
ggboxplot(Ig_entropy, x = "Tumor_type_4categ", y = "value",facet.by = "variable",color = "Tumor_type_4categ",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=value, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))
dev.off()

TR_entropy<-melt(PAAD.repertoire.diversity_tclones[,c("TCGA_sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG","Tumor_type_4categ")])
TR_entropy<-TR_entropy[which(TR_entropy$value!=0),]
tiff("Results/boxplot_TR_entropy_recon_PAAD.tiff",res=300,h=2500,w=3500)
ggboxplot(TR_entropy, x = "Tumor_type_4categ", y = "value",facet.by = "variable",color = "Tumor_type_4categ",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=value, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))
dev.off()

TR_entropy<-melt(PAAD.repertoire.diversity_tclones[,c("TCGA_sample","entropy_TRA","entropy_TRB","entropy_TRD","entropy_TRG","Tumor_type_4categ")])
TR_entropy<-TR_entropy[which(TR_entropy$value!=0),]
tiff("Results/boxplot_TR_entropy_PAAD.tiff",res=300,h=2500,w=3500)
ggboxplot(TR_entropy, x = "Tumor_type_4categ", y = "value",facet.by = "variable",color = "Tumor_type_4categ",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=value, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))
dev.off()

kappa_lambda<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","KappaLambda_ratio_expression_filter","Tumor_type_4categ")])
tiff("Results/boxplot_kappa_lambda_PAAD_filter.tiff",res=300,h=2500,w=3000)
ggboxplot(kappa_lambda, x = "Tumor_type_4categ", y = "value",facet.by = "variable",color = "Tumor_type_4categ",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=value, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))
dev.off()

alpha_beta_ratio<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","Alpha_Beta_ratio_expression_filter","Tumor_type_4categ")])
tiff("Results/boxplot_alpha_beta_PAAD_filter.tiff",res=300,h=2500,w=3500)
ggboxplot(alpha_beta_ratio, x = "Tumor_type_4categ", y = "value",facet.by = "variable",color = "Tumor_type_4categ",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=value, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))
dev.off()

Ig_cdr3length<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","cdr3_length_IGH","cdr3_length_IGK","cdr3_length_IGL","Tumor_type_4categ")])
tiff("Results/boxplot_Ig_cdr3_length_filter.tiff",res=300,h=2500,w=3000)
ggboxplot(Ig_cdr3length, x = "Tumor_type_4categ", y = "value",facet.by = "variable",color = "Tumor_type_4categ",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=value, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))
dev.off()

TCR_cdr3length<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","cdr3_length_TRA","cdr3_length_TRB","cdr3_length_TRD","cdr3_length_TRG","Tumor_type_4categ")])
TCR_cdr3length<-TCR_cdr3length[which(TCR_cdr3length$value!=0),]
tiff("Results/boxplot_TCR_cdr3_length_filter.tiff",res=300,h=2500,w=3000)
ggboxplot(TCR_cdr3length, x = "Tumor_type_4categ", y = "value",facet.by = "variable",color = "Tumor_type_4categ",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=value, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))
dev.off()

#########################
### 2. GTEX Pancreas ###
########################
load("Data/GTEx/Pancreas/GTEx_FullData.Rdata")

brewer.pal(3,name = "Accent")
brewer.pal(3,name = "Pastel1")
cols=c( "#7FC97F", "#FDC086", "#BEAED4" ,"#B3CDE3")
#Barplot
tiff("Results/barplot_Treads_GTEX_pancreas.tiff",res=300,h=2500,w=4000)
barplot(Pancreas.repertoire.diversity$T_Reads,col=cols[4],main="Number of T-Reads",xlab = "Samples", ylab = "Reads",las=2)
abline(h=100)
dev.off()

tiff("Results/barplot_Igreads_GTEX_pancreas.tiff",res=300,h=2500,w=4000)
barplot(Pancreas.repertoire.diversity$Ig_Reads,col=cols[4],main="Number of Ig-Reads",xlab = "Samples", ylab = "Reads",las=2)
abline(h=100)
dev.off()
tiff("Results/Corr_plot_reads_GTEX_pancreas.tiff",res=300,h=2500,w=4000)
plot(Pancreas.repertoire.diversity$read_count,Pancreas.repertoire.diversity$totalReads,col=cols[4],pch=19,
     xlab = "Total B- and T- reads aligned", ylab = "Total sequencing reads",
     main = paste0("rho = ",round(cor(Pancreas.repertoire.diversity$read_count,Pancreas.repertoire.diversity$totalReads),2)))
dev.off()

summary(Pancreas.repertoire.diversity$Ig_Reads)
summary(Pancreas.repertoire.diversity$T_Reads)
summary(Pancreas.repertoire.diversity$totalReads)

##Filter by number of reads and clones <100
Pancreas.repertoire.diversity_treads<-Pancreas.repertoire.diversity[which(Pancreas.repertoire.diversity$T_Reads>100),]
Pancreas.repertoire.diversity_Igreads<-Pancreas.repertoire.diversity[which(Pancreas.repertoire.diversity$Ig_Reads>100),]
Pancreas.repertoire.diversity$T_clones<-Pancreas.repertoire.diversity$clones_recon_TRA+Pancreas.repertoire.diversity$clones_recon_TRB+
  Pancreas.repertoire.diversity$clones_recon_TRD+Pancreas.repertoire.diversity$clones_recon_TRG
Pancreas.repertoire.diversity$Ig_clones<-Pancreas.repertoire.diversity$clones_recon_IGH+Pancreas.repertoire.diversity$clones_recon_IGK+
  Pancreas.repertoire.diversity$clones_recon_IGL
Pancreas.repertoire.diversity_tclones<-Pancreas.repertoire.diversity[which(Pancreas.repertoire.diversity$T_clones>100),]
Pancreas.repertoire.diversity_Igclones<-Pancreas.repertoire.diversity[which(Pancreas.repertoire.diversity$Ig_clones>100),]

####Summary plots
#IgExpression
Ig_expr_PAAD<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","IGH_expression","IGK_expression","IGL_expression","Tumor_type_3categ")])
Pancreas.repertoire.diversity_Igreads$TCGA_sample<-rownames(Pancreas.repertoire.diversity_Igreads)
Pancreas.repertoire.diversity_Igreads$Tumor_type_3categ<-rep("Normal_pancreas")
Ig_expr_Pancreas<-melt(Pancreas.repertoire.diversity_Igreads[,c("TCGA_sample","IGH_expression","IGK_expression","IGL_expression","Tumor_type_3categ")])

Ig_expr<-rbind(Ig_expr_PAAD,Ig_expr_Pancreas)
Ig_expr<-Ig_expr[which(Ig_expr$value!=0),]
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_TCGA_GTEX_Pancreas.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_expr, x = "Tumor_type_3categ", y = "value",facet.by = "variable",color = "Tumor_type_3categ",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_3categ, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2],cols[3],cols[4]), labels = c("normal_pancreas (TCGA)", "pseudonormal_pancreas (TCGA)", 
                                                                              "tumor_pancreas (TCGA)","normal_pancreas (GTEx)")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","Normal_pancreas"),c("pseudonormal_pancreas","Normal_pancreas"),c("Tumor_pancreas","Normal_pancreas")))
dev.off()

#Entropy
Ig_entropy_PAAD<-melt(PAAD.repertoire.diversity_Igclones[,c("TCGA_sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL","Tumor_type_3categ")])
Pancreas.repertoire.diversity_Igclones$TCGA_sample<-rownames(Pancreas.repertoire.diversity_Igclones)
Pancreas.repertoire.diversity_Igclones$Tumor_type_3categ<-rep("Normal_pancreas")
Ig_entropy_Pancreas<-melt(Pancreas.repertoire.diversity_Igclones[,c("TCGA_sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL","Tumor_type_3categ")])

Ig_expr<-rbind(Ig_entropy_PAAD,Ig_entropy_Pancreas)
Ig_expr<-Ig_expr[which(Ig_expr$value!=0),]
tiff("Results/boxplot_Ig_entropy_TCGA_GTEX_Pancreas.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_expr, x = "Tumor_type_3categ", y = "value",facet.by = "variable",color = "Tumor_type_3categ",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_3categ, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2],cols[3],cols[4]), labels = c("normal_pancreas (TCGA)", "pseudonormal_pancreas (TCGA)", 
                                                                              "tumor_pancreas (TCGA)","normal_pancreas (GTEx)")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","Normal_pancreas"),c("pseudonormal_pancreas","Normal_pancreas"),c("Tumor_pancreas","Normal_pancreas")))
dev.off()

#Texpression
T_expr_PAAD<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression","Tumor_type_3categ")])
Pancreas.repertoire.diversity_treads$TCGA_sample<-rownames(Pancreas.repertoire.diversity_treads)
Pancreas.repertoire.diversity_treads$Tumor_type_3categ<-rep("Normal_pancreas")
T_expr_Pancreas<-melt(Pancreas.repertoire.diversity_treads[,c("TCGA_sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression","Tumor_type_3categ")])

T_expr<-rbind(T_expr_PAAD,T_expr_Pancreas)
T_expr<-T_expr[which(T_expr$value!=0),]
T_expr$value<-log10(T_expr$value)
tiff("Results/boxplot_T_expression_TCGA_GTEX_Pancreas.tiff",res=300,h=2500,w=3500)
ggboxplot(T_expr, x = "Tumor_type_3categ", y = "value",facet.by = "variable",color = "Tumor_type_3categ",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_3categ, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2],cols[3],cols[4]), labels = c("normal_pancreas (TCGA)", "pseudonormal_pancreas (TCGA)", 
                                                                              "tumor_pancreas (TCGA)","normal_pancreas (GTEx)")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","Normal_pancreas"),c("pseudonormal_pancreas","Normal_pancreas"),c("Tumor_pancreas","Normal_pancreas")))
dev.off()


#Entropy (T)
T_entropy_PAAD<-melt(PAAD.repertoire.diversity[,c("TCGA_sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG","Tumor_type_3categ")])
Pancreas.repertoire.diversity_tclones$TCGA_sample<-rownames(Pancreas.repertoire.diversity_tclones)
Pancreas.repertoire.diversity_tclones$Tumor_type_3categ<-rep("Normal_pancreas")
T_entropy_Pancreas<-melt(Pancreas.repertoire.diversity_tclones[,c("TCGA_sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG","Tumor_type_3categ")])

T_expr<-rbind(T_entropy_PAAD,T_entropy_Pancreas)
T_expr<-T_expr[which(T_expr$value!=0),]
tiff("Results/boxplot_T_entropy_TCGA_GTEX_Pancreas.tiff",res=300,h=2500,w=3500)
ggboxplot(T_expr, x = "Tumor_type_3categ", y = "value",facet.by = "variable",color = "Tumor_type_3categ",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_3categ, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2],cols[3],cols[4]), labels = c("normal_pancreas (TCGA)", "pseudonormal_pancreas (TCGA)", 
                                                                              "tumor_pancreas (TCGA)","normal_pancreas (GTEx)")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","Normal_pancreas"),c("pseudonormal_pancreas","Normal_pancreas"),c("Tumor_pancreas","Normal_pancreas")))
dev.off()

#KappaLambda
kappa_lambda_PAAD<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","KappaLambda_ratio_expression","Tumor_type_3categ")])
kappa_lambda_Pancreas<-melt(Pancreas.repertoire.diversity_Igreads[,c("TCGA_sample","KappaLambda_ratio_expression","Tumor_type_3categ")])

kappa_lambda<-rbind(kappa_lambda_PAAD,kappa_lambda_Pancreas)
kappa_lambda<-kappa_lambda[which(kappa_lambda$value!=0),]
tiff("Results/boxplot_kappa_lambda_TCGA_GTEX_Pancreas.tiff",res=300,h=2500,w=3500)
ggboxplot(kappa_lambda, x = "Tumor_type_3categ", y = "value",facet.by = "variable",color = "Tumor_type_3categ",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_3categ, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2],cols[3],cols[4]), labels = c("normal_pancreas (TCGA)", "pseudonormal_pancreas (TCGA)", 
                                                                              "tumor_pancreas (TCGA)","normal_pancreas (GTEx)")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","Normal_pancreas"),c("pseudonormal_pancreas","Normal_pancreas"),c("Tumor_pancreas","Normal_pancreas")))
dev.off()


#alphabeta
alpha_beta_ratio_PAAD<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","Alpha_Beta_ratio_expression","Tumor_type_3categ")])
alpha_beta_ratio_Pancreas<-melt(Pancreas.repertoire.diversity_treads[,c("TCGA_sample","Alpha_Beta_ratio_expression","Tumor_type_3categ")])

alpha_beta_ratio<-rbind(alpha_beta_ratio_PAAD,alpha_beta_ratio_Pancreas)
alpha_beta_ratio<-alpha_beta_ratio[which(alpha_beta_ratio$value!=1),]
tiff("Results/boxplot_alpha_beta_TCGA_GTEX_Pancreas.tiff",res=300,h=2500,w=3500)
ggboxplot(alpha_beta_ratio, x = "Tumor_type_3categ", y = "value",facet.by = "variable",color = "Tumor_type_3categ",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_3categ, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2],cols[3],cols[4]), labels = c("normal_pancreas (TCGA)", "pseudonormal_pancreas (TCGA)", 
                                                                              "tumor_pancreas (TCGA)","normal_pancreas (GTEx)")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","Normal_pancreas"),c("pseudonormal_pancreas","Normal_pancreas"),c("Tumor_pancreas","Normal_pancreas")))
dev.off()

###############################
### 3. Validation Pancreas ####
###############################
load("Data/Pancreas_Validation/Pancreas_Validation_FullData.Rdata")


brewer.pal(3,name = "Accent")
brewer.pal(3,name = "Pastel1")
cols=c( "#7FC97F", "#FDC086", "#BEAED4" ,"#B3CDE3")

summary(Pancreas.Validation.repertoire.diversity$Ig_Reads[which(Pancreas.Validation.repertoire.diversity$tissue=="pancreas tumor")])
summary(Pancreas.Validation.repertoire.diversity$T_Reads)
summary(Pancreas.Validation.repertoire.diversity$totalReads)

##Filter by number of reads and clones <100
Pancreas.Validation.repertoire.diversity_treads<-Pancreas.Validation.repertoire.diversity[which(Pancreas.Validation.repertoire.diversity$T_Reads>100),]
Pancreas.Validation.repertoire.diversity_Igreads<-Pancreas.Validation.repertoire.diversity[which(Pancreas.Validation.repertoire.diversity$Ig_Reads>100),]
Pancreas.Validation.repertoire.diversity$T_clones<-Pancreas.Validation.repertoire.diversity$clones_recon_TRA+Pancreas.Validation.repertoire.diversity$clones_recon_TRB+
  Pancreas.Validation.repertoire.diversity$clones_recon_TRD+Pancreas.Validation.repertoire.diversity$clones_recon_TRG
Pancreas.Validation.repertoire.diversity$Ig_clones<-Pancreas.Validation.repertoire.diversity$clones_recon_IGH+Pancreas.Validation.repertoire.diversity$clones_recon_IGK+
  Pancreas.Validation.repertoire.diversity$clones_recon_IGL
Pancreas.Validation.repertoire.diversity_tclones<-Pancreas.Validation.repertoire.diversity[which(Pancreas.Validation.repertoire.diversity$T_clones>100),]
Pancreas.Validation.repertoire.diversity_Igclones<-Pancreas.Validation.repertoire.diversity[which(Pancreas.Validation.repertoire.diversity$Ig_clones>100),]

###T markers
#Barplot
Pancreas.Validation.repertoire.diversity$T_Reads<-Pancreas.Validation.repertoire.diversity$TRA+Pancreas.Validation.repertoire.diversity$TRB+
  Pancreas.Validation.repertoire.diversity$TRD+Pancreas.Validation.repertoire.diversity$TRG
tiff("Results/barplot_Treads_PancreasValidation.tiff",res=300,h=2500,w=4000)
barplot(Pancreas.Validation.repertoire.diversity$T_Reads,col=cols[Pancreas.Validation.repertoire.diversity$tissue],main="Number of T-Reads",xlab = "Samples", ylab = "Reads",las=2)
abline(h=100)
legend("topright", legend=levels(Pancreas.Validation.repertoire.diversity$tissue),col=cols,pch=15, cex=0.8)
dev.off()

tiff("Results/Corr_plot_reads_PancreasValidation.tiff",res=300,h=2500,w=4000)
plot(Pancreas.Validation.repertoire.diversity$read_count,Pancreas.Validation.repertoire.diversity$totalReads,col=cols[Pancreas.Validation.repertoire.diversity$tissue],pch=19,
     xlab = "Total B- and T- reads aligned", ylab = "Total sequencing reads",
     main = paste0("rho = ",round(cor(Pancreas.Validation.repertoire.diversity$read_count,Pancreas.Validation.repertoire.diversity$totalReads),2)))
legend("bottomright", legend=levels(Pancreas.Validation.repertoire.diversity$tissue),col=cols,pch=15, cex=0.8)
dev.off()

##Ig markers
#Barplot
tiff("Results/barplot_Igreads_PancreasValidation.tiff",res=300,h=2500,w=4000)
barplot(Pancreas.Validation.repertoire.diversity$IG_Reads,col=cols[Pancreas.Validation.repertoire.diversity$tissue],main="Number of Ig-Reads",xlab = "Samples", ylab = "Reads",las=2)
abline(h=100)
legend("topright", legend=levels(Pancreas.Validation.repertoire.diversity$tissue),col=cols,pch=15, cex=0.8)
dev.off()

####Summary plots
cols=c( "#7FC97F", "#FDC086", "#BEAED4" ,"#FBB4AE")
#IgExpression
Ig_expr<-melt(Pancreas.Validation.repertoire.diversity_Igreads[,c("sample","IGH_expression","IGK_expression","IGL_expression","tissue")])
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_PancreasValidation.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_expr, x = "tissue", y = "value",facet.by = "variable",color = "tissue",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=tissue, y=value,color=tissue), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[3]), labels = c("normal pancreas (val)","tumor pancreas (val)")) +
  stat_compare_means(
    comparisons =list(c("normal pancreas","pancreas tumor")))
dev.off()

#Entropy
Ig_entropy<-melt(Pancreas.Validation.repertoire.diversity_Igclones[,c("sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL","tissue")])
tiff("Results/boxplot_Ig_entropy_PancreasValidation.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_entropy, x = "tissue", y = "value",facet.by = "variable",color = "tissue",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=tissue, y=value,color=tissue), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[3]), labels = c("normal pancreas (val)","tumor pancreas (val)")) +
  stat_compare_means(
    comparisons =list(c("normal pancreas","pancreas tumor")))
dev.off()

#Texpression
T_expr<-melt(Pancreas.Validation.repertoire.diversity_treads[,c("sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression","tissue")])
T_expr<-T_expr[which(T_expr$value!=0),]
T_expr$value<-log10(T_expr$value)
tiff("Results/boxplot_T_expression_PancreasValidation.tiff",res=300,h=2500,w=3500)
ggboxplot(T_expr, x = "tissue", y = "value",facet.by = "variable",color = "tissue",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=tissue, y=value,color=tissue), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[3]), labels = c("normal pancreas (val)","tumor pancreas (val)")) +
  stat_compare_means(
    comparisons =list(c("normal pancreas","pancreas tumor")))
dev.off()


#Entropy (T)
T_entropy<-melt(Pancreas.Validation.repertoire.diversity_tclones[,c("sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG","tissue")])
T_entropy<-T_entropy[which(T_entropy$value!=0),]
tiff("Results/boxplot_T_entropy_PancreasValidation.tiff",res=300,h=2500,w=3500)
ggboxplot(T_entropy, x = "tissue", y = "value",facet.by = "variable",color = "tissue",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=tissue, y=value,color=tissue), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[3]), labels = c("normal pancreas (val)","tumor pancreas (val)")) +
  stat_compare_means(
    comparisons =list(c("normal pancreas","pancreas tumor")))
dev.off()

#KappaLambda
kappa_lambda<-melt(Pancreas.Validation.repertoire.diversity_Igreads[,c("sample","KappaLambda_ratio_expression","tissue")])
tiff("Results/boxplot_kappa_lambda_PancreasValidation.tiff",res=300,h=2500,w=3500)
ggboxplot(kappa_lambda, x = "tissue", y = "value",facet.by = "variable",color = "tissue",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=tissue, y=value,color=tissue), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[3]), labels = c("normal pancreas (val)","tumor pancreas (val)")) +
  stat_compare_means(
    comparisons =list(c("normal pancreas","pancreas tumor")))
dev.off()


#alphabeta
alpha_beta_ratio<-melt(Pancreas.Validation.repertoire.diversity_treads[,c("sample","Alpha_Beta_ratio_expression","tissue")])
tiff("Results/boxplot_alpha_beta_PancreasValidation.tiff",res=300,h=2500,w=3500)
ggboxplot(alpha_beta_ratio, x = "tissue", y = "value",facet.by = "variable",color = "tissue",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=tissue, y=value,color=tissue), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[3]), labels = c("normal pancreas (val)","tumor pancreas (val)")) +
  stat_compare_means(
    comparisons =list(c("normal pancreas","pancreas tumor")))
dev.off()

#########################################################
### 4. TCGA vs. Validation with clones call together ####
#########################################################
load("Data/PAAD_Val/PAAD_VAL_FullData.Rdata")

###Log transform the expression
PAAD.VAL.repertoire.diversity$IGH_expression_filter_log<-log10(PAAD.VAL.repertoire.diversity$IGH_expression_filter)
PAAD.VAL.repertoire.diversity$IGK_expression_filter_log<-log10(PAAD.VAL.repertoire.diversity$IGK_expression_filter)
PAAD.VAL.repertoire.diversity$IGL_expression_filter_log<-log10(PAAD.VAL.repertoire.diversity$IGL_expression_filter)
PAAD.VAL.repertoire.diversity$TRA_expression_filter_log<-log10(PAAD.VAL.repertoire.diversity$TRA_expression_filter)
PAAD.VAL.repertoire.diversity$TRB_expression_filter_log<-log10(PAAD.VAL.repertoire.diversity$TRB_expression_filter)
PAAD.VAL.repertoire.diversity$TRD_expression_filter_log<-log10(PAAD.VAL.repertoire.diversity$TRD_expression_filter)
PAAD.VAL.repertoire.diversity$TRG_expression_filter_log<-log10(PAAD.VAL.repertoire.diversity$TRG_expression_filter)

##Filter by number of reads and clones <100
PAAD.VAL.repertoire.diversity_treads<-PAAD.VAL.repertoire.diversity[which(PAAD.VAL.repertoire.diversity$T_Reads_filter>100),]
PAAD.VAL.repertoire.diversity_Igreads<-PAAD.VAL.repertoire.diversity[which(PAAD.VAL.repertoire.diversity$Ig_Reads_filter>100),]
PAAD.VAL.repertoire.diversity$T_clones<-PAAD.VAL.repertoire.diversity$clones_TRA+PAAD.VAL.repertoire.diversity$clones_TRB+
  PAAD.VAL.repertoire.diversity$clones_TRD+PAAD.VAL.repertoire.diversity$clones_TRG
PAAD.VAL.repertoire.diversity$Ig_clones<-PAAD.VAL.repertoire.diversity$clones_IGH+PAAD.VAL.repertoire.diversity$clones_IGK+
  PAAD.VAL.repertoire.diversity$clones_IGL
PAAD.VAL.repertoire.diversity_tclones<-PAAD.VAL.repertoire.diversity[which(PAAD.VAL.repertoire.diversity$T_clones>100),]
PAAD.VAL.repertoire.diversity_Igclones<-PAAD.VAL.repertoire.diversity[which(PAAD.VAL.repertoire.diversity$Ig_clones>100),]

##Normalization measures
cols=brewer.pal(3,name = "Set1")

tiff("Results/boxplot_norm_expression_TCGA_VAL.tiff",res=300,h=2000,w=3500)
boxplot(t(PAAD.VAL.repertoire.diversity[, c("IGH_expression_filter_log",   "IGK_expression_filter_log",   "IGL_expression_filter_log", 
                                            "TRA_expression_filter_log",  "TRB_expression_filter_log", "TRD_expression_filter_log", 
                                            "TRG_expression_filter_log")]),
        col=c(cols[1],cols[2],cols[1],cols[1],cols[2])[PAAD.VAL.repertoire.diversity$outcome],main="norm_expression")
legend("topright",c("TCGA","Validation"),fill=cols)
dev.off()

tiff("Results/boxplot_norm_expression_vst_TCGA_VAL.tiff",res=300,h=2000,w=3500)
boxplot(t(PAAD.VAL.repertoire.diversity[, c("IGH_expression_vst",   "IGK_expression_vst",   "IGL_expression_vst", 
                                             "TRA_expression_vst",  "TRB_expression_vst", "TRD_expression_vst", "TRG_expression_vst")]),
        col=c(cols[1],cols[2],cols[1],cols[1],cols[2])[PAAD.VAL.repertoire.diversity$outcome],main="norm_expression_vst")
legend("topright",c("TCGA","GTEx"),fill=cols)
dev.off()

tiff("Results/boxplot_norm_entropy_TCGA_VAL.tiff",res=300,h=2000,w=3500)
boxplot(t(PAAD.VAL.repertoire.diversity[, c("entropy_IGH", "entropy_IGK", "entropy_IGL", "entropy_TRA",
                                            "entropy_TRB", "entropy_TRD", "entropy_TRG")]),
        col=c(cols[1],cols[2],cols[1],cols[1],cols[2])[PAAD.VAL.repertoire.diversity$outcome],main="norm_entropy")
legend("topright",c("TCGA","Validation"),fill=cols)
dev.off()

#Barplot
tiff("Results/barplot_Igreads.tiff",res=300,h=2500,w=4000)
barplot(t(as.matrix(PAAD.VAL.repertoire.diversity_Igclones[, c("IGH_filter",   "IGK_filter",   "IGL_filter")]))
          ,col=c(cols[1],cols[2],cols[1],cols[1],cols[2])[PAAD.VAL.repertoire.diversity_Igclones$outcome],main="Number of Ig-Reads",xlab = "Samples", ylab = "Reads",las=2)

barplot(t(as.matrix(PAAD.VAL.repertoire.diversity_Igclones[, c("clones_IGH",   "clones_IGK",   "clones_IGL")]))
        ,col=c(cols[1],cols[2],cols[1],cols[1],cols[2])[PAAD.VAL.repertoire.diversity_Igclones$outcome],main="Number of Ig-Reads",xlab = "Samples", ylab = "Reads",las=2)
dev.off()

tiff("Results/Corr_plot_reads.tiff",res=300,h=2500,w=4000)
plot(PAAD.VAL.repertoire.diversity_Igclones$Ig_clones,PAAD.VAL.repertoire.diversity_Igclones$totalReads,
     col=c(cols[1],cols[2],cols[1],cols[1],cols[2])[PAAD.VAL.repertoire.diversity_Igclones$outcome],pch=19,
     xlab = "Total Ig clones", ylab = "Total sequencing reads",
     main = paste0("rho = ",round(cor(PAAD.VAL.repertoire.diversity_Igclones$Ig_clones,PAAD.VAL.repertoire.diversity_Igclones$totalReads),2)))
dev.off()


##Boxplots
brewer.pal(3,name = "Accent")
cols=c( "#7FC97F", "#FDC086", "#BEAED4")

#Ig expression
Ig_expr<-melt(PAAD.VAL.repertoire.diversity_Igclones[,c("sample","IGH_expression","IGK_expression","IGL_expression","outcome")])
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_expr, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#Ig expression vst
Ig_expr<-melt(PAAD.VAL.repertoire.diversity_Igclones[,c("sample","IGH_expression_vst","IGK_expression_vst","IGL_expression_vst","outcome")])
tiff("Results/boxplot_Ig_expression_vst_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_expr, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#Entropy
Ig_entropy<-melt(PAAD.VAL.repertoire.diversity_Igclones[,c("sample","entropy_IGH","entropy_IGK","entropy_IGL","outcome")])
tiff("Results/boxplot_Ig_entropy_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
      "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                       c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                       c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#Clones
Ig_clones<-melt(PAAD.VAL.repertoire.diversity_Igclones[,c("sample","clones_IGH","clones_IGK","clones_IGL","outcome")])
tiff("Results/boxplot_Ig_clones_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_clones, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#Clones/reads
PAAD.VAL.repertoire.diversity_Igclones$IGH_ratio<-PAAD.VAL.repertoire.diversity_Igclones$clones_IGH/PAAD.VAL.repertoire.diversity_Igclones$IGH_filter
PAAD.VAL.repertoire.diversity_Igclones$IGK_ratio<-PAAD.VAL.repertoire.diversity_Igclones$clones_IGK/PAAD.VAL.repertoire.diversity_Igclones$IGK_filter
PAAD.VAL.repertoire.diversity_Igclones$IGL_ratio<-PAAD.VAL.repertoire.diversity_Igclones$clones_IGL/PAAD.VAL.repertoire.diversity_Igclones$IGL_filter

Ig_clones<-melt(PAAD.VAL.repertoire.diversity_Igclones[,c("sample","IGH_ratio","IGK_ratio","IGL_ratio","outcome")])
tiff("Results/boxplot_Ig_clones_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_clones, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#Entropy recon
Ig_entropy<-melt(PAAD.VAL.repertoire.diversity_Igclones[,c("sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL","outcome")])
tiff("Results/boxplot_Ig_entropy_recon_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#T expression
T_expr<-melt(PAAD.VAL.repertoire.diversity_tclones[,c("sample","TRA_expression_filter","TRB_expression_filter","TRD_expression_filter","TRG_expression_filter","outcome")])
T_expr$value<-log10(T_expr$value)
tiff("Results/boxplot_T_expression_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(T_expr, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#T expression VST
T_expr<-melt(PAAD.VAL.repertoire.diversity_tclones[,c("sample","TRA_expression_vst","TRB_expression_vst","TRD_expression_vst","TRG_expression_vst","outcome")])
T_expr$value<-log10(T_expr$value)
tiff("Results/boxplot_T_expression_TCGA_VAL_VST.tiff",res=300,h=2500,w=3500)
ggboxplot(T_expr, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#Entropy
T_entropy<-melt(PAAD.VAL.repertoire.diversity_tclones[,c("sample","entropy_TRA","entropy_TRB","entropy_TRD","entropy_TRG","outcome")])
tiff("Results/boxplot_T_entropy_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(T_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()


#Entropy
T_entropy<-melt(PAAD.VAL.repertoire.diversity_tclones[,c("sample","entropy_TRA","entropy_TRB","entropy_TRD","entropy_TRG","outcome")])
tiff("Results/boxplot_T_entropy_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(T_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
ddev.off()

#Clones
T_clones<-melt(PAAD.VAL.repertoire.diversity_tclones[,c("sample","clones_TRA","clones_TRB","clones_TRD","clones_TRG","outcome")])
tiff("Results/boxplot_T_clones_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(T_clones, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#Entropy recon
T_entropy<-melt(PAAD.VAL.repertoire.diversity_tclones[,c("sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG","outcome")])
tiff("Results/boxplot_T_entropy_recon_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(T_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#Clones
T_clones<-melt(PAAD.VAL.repertoire.diversity_tclones[,c("sample","clones_recon_TRA","clones_recon_TRB","clones_recon_TRD","clones_recon_TRG","outcome")])
tiff("Results/boxplot_T_clones_recon_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(T_clones, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#Alpha_Beta_ratio_expression
Alpha_Beta_ratio_expression<-melt(PAAD.VAL.repertoire.diversity_tclones[,c("sample","Alpha_Beta_ratio_expression","outcome")])
tiff("Results/boxplot_Alpha_Beta_ratio_expression_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(Alpha_Beta_ratio_expression, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#KappaLambda_ratio_expression
KappaLambda_ratio_expression<-melt(PAAD.VAL.repertoire.diversity_Igclones[,c("sample","KappaLambda_ratio_expression","outcome")])
tiff("Results/boxplot_KappaLambda_ratio_expression_TCGA_VAL.tiff",res=300,h=2500,w=3500)
ggboxplot(KappaLambda_ratio_expression, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1],cols[1],cols[2],cols[3],cols[3]), labels = c("normal-pancreas (TCGA)",
                                                                                     "normal-pancreas (Val)", "pseudonormal-pancreas (TCGA)", "tumor-pancreas (TCGA)", "tumor-pancreas (Val)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (TCGA)","normal-pancreas (Val)"),c("tumor-pancreas (TCGA)","tumor-pancreas (Val)"),
                      c("normal-pancreas (TCGA)","tumor-pancreas (TCGA)"),c("normal-pancreas (Val)","tumor-pancreas (Val)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

#########################################################
### 5. TCGA vs. GTEX pancreas with clones call together ####
#########################################################
load("Data/PAAD_GTEx/PAAD_GTEx_FullData.Rdata")
PAAD.GTEx.repertoire.diversity$outcome<-factor(PAAD.GTEx.repertoire.diversity$outcome)

###Log transform the expression
PAAD.GTEx.repertoire.diversity$IGH_expression_filter_log<-log10(PAAD.GTEx.repertoire.diversity$IGH_expression_filter+0.00000001)
PAAD.GTEx.repertoire.diversity$IGK_expression_filter_log<-log10(PAAD.GTEx.repertoire.diversity$IGK_expression_filter+0.00000001)
PAAD.GTEx.repertoire.diversity$IGL_expression_filter_log<-log10(PAAD.GTEx.repertoire.diversity$IGL_expression_filter+0.00000001)
PAAD.GTEx.repertoire.diversity$TRA_expression_filter_log<-log10(PAAD.GTEx.repertoire.diversity$TRA_expression_filter+0.00000001)
PAAD.GTEx.repertoire.diversity$TRB_expression_filter_log<-log10(PAAD.GTEx.repertoire.diversity$TRB_expression_filter+0.00000001)
PAAD.GTEx.repertoire.diversity$TRD_expression_filter_log<-log10(PAAD.GTEx.repertoire.diversity$TRD_expression_filter+0.00000001)
PAAD.GTEx.repertoire.diversity$TRG_expression_filter_log<-log10(PAAD.GTEx.repertoire.diversity$TRG_expression_filter+0.00000001)

##Filter by number of reads and clones <100
PAAD.GTEx.repertoire.diversity_treads<-PAAD.GTEx.repertoire.diversity[which(PAAD.GTEx.repertoire.diversity$T_Reads_filter>100),]
PAAD.GTEx.repertoire.diversity_Igreads<-PAAD.GTEx.repertoire.diversity[which(PAAD.GTEx.repertoire.diversity$Ig_Reads_filter>100),]
PAAD.GTEx.repertoire.diversity$T_clones<-PAAD.GTEx.repertoire.diversity$clones_TRA+PAAD.GTEx.repertoire.diversity$clones_TRB+
  PAAD.GTEx.repertoire.diversity$clones_TRD+PAAD.GTEx.repertoire.diversity$clones_TRG
PAAD.GTEx.repertoire.diversity$Ig_clones<-PAAD.GTEx.repertoire.diversity$clones_IGH+PAAD.GTEx.repertoire.diversity$clones_IGK+
  PAAD.GTEx.repertoire.diversity$clones_IGL
PAAD.GTEx.repertoire.diversity_tclones<-PAAD.GTEx.repertoire.diversity[which(PAAD.GTEx.repertoire.diversity$T_clones>100),]
PAAD.GTEx.repertoire.diversity_Igclones<-PAAD.GTEx.repertoire.diversity[which(PAAD.GTEx.repertoire.diversity$Ig_clones>100),]

PAAD.GTEx.repertoire.diversity$Ig_clone_down<-PAAD.GTEx.repertoire.diversity$clones_IGH_down+PAAD.GTEx.repertoire.diversity$clones_IGK_down+
  PAAD.GTEx.repertoire.diversity$clones_IGL_down
PAAD.GTEx.repertoire.diversity_Igclones_down<-PAAD.GTEx.repertoire.diversity[which(PAAD.GTEx.repertoire.diversity$Ig_clone_down>100),]

#Normalization measures
cols=brewer.pal(3,name = "Set1")

tiff("Results/boxplot_norm_expression_TCGA_GTEX.tiff",res=300,h=2000,w=3500)
boxplot(t(PAAD.GTEx.repertoire.diversity[, c("IGH_expression_filter_log",   "IGK_expression_filter_log",   "IGL_expression_filter_log", 
                                            "TRA_expression_filter_log",  "TRB_expression_filter_log", "TRD_expression_filter_log", "TRG_expression_filter_log")]),
        col=c(cols[2],cols[1],cols[1],cols[1])[PAAD.GTEx.repertoire.diversity$outcome],main="norm_expression")
legend("topright",c("TCGA","GTEx"),fill=cols)
dev.off()
tiff("Results/boxplot_norm_expression_vst_TCGA_GTEx.tiff",res=300,h=2000,w=3500)
boxplot(t(PAAD.GTEx.repertoire.diversity[, c("IGH_expression_vst",   "IGK_expression_vst",   "IGL_expression_vst", 
                                            "TRA_expression_vst",  "TRB_expression_vst", "TRD_expression_vst", "TRG_expression_vst")]),
        col=c(cols[2],cols[1],cols[1],cols[1])[PAAD.GTEx.repertoire.diversity$outcome],main="norm_expression_vst")
legend("topright",c("TCGA","GTEx"),fill=cols)
dev.off()
tiff("Results/boxplot_norm_entropy_TCGA_GTEx.tiff",res=300,h=2000,w=3500)
boxplot(t(PAAD.GTEx.repertoire.diversity[, c("entropy_IGH", "entropy_IGK", "entropy_IGL", "entropy_TRA",
                                            "entropy_TRB", "entropy_TRD", "entropy_TRG")]),
        col=c(cols[2],cols[1],cols[1],cols[1])[PAAD.GTEx.repertoire.diversity$outcome],main="norm_entropy")
legend("topright",c("TCGA","GTEx"),fill=cols)
dev.off()

tiff("Results/boxplot_norm_entropy_recon_TCGA_GTEx.tiff",res=300,h=2000,w=3500)
boxplot(t(PAAD.GTEx.repertoire.diversity[,c("entropy_recon_IGH", "entropy_recon_IGK", "entropy_recon_IGL", 
                                           "entropy_recon_TRA","entropy_recon_TRB", "entropy_recon_TRD", "entropy_recon_TRG")]),
        col=c(cols[2],cols[1],cols[1],cols[1])[PAAD.GTEx.repertoire.diversity$outcome],main="norm_entropy_recon")
legend("topright",c("TCGA","GTEx"),fill=cols)
dev.off()


########Analysis 

brewer.pal(3,name = "Accent")
cols=c( "#7FC97F","#386CB0","#FDC086", "#BEAED4")

#Ig expression
Ig_expr<-melt(PAAD.GTEx.repertoire.diversity_Igclones[,c("sample","IGH_expression_filter","IGK_expression_filter",
                                                         "IGL_expression_filter","outcome")])
Ig_expr<-Ig_expr[which(Ig_expr$value!=0),]
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_TCGA_GTEX.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_expr, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                             "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  #stat_compare_means(label = "p.format")
  stat_compare_means(
  comparisons =list(c("normal-pancreas (GTEx)","normal-pancreas (TCGA)"),c("normal-pancreas (GTEx)","pseudonormal-pancreas (TCGA)"),
                    c("normal-pancreas (GTEx)","tumor-pancreas (TCGA)"),
                    c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))

dev.off()

summary(glm(IGH_expression ~ outcome + sex + age + race, data = PAAD.GTEx.repertoire.diversity_Igreads))
   
#Ig expression vst
Ig_expr<-melt(PAAD.GTEx.repertoire.diversity_Igclones[,c("sample","IGH_expression_vst","IGK_expression_vst","IGL_expression_vst","outcome")])
Ig_expr<-Ig_expr[which(Ig_expr$value!=0),]
#Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_TCGA_GTEX_VST.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_expr, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                  "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  stat_compare_means(
    comparisons =list(c("normal-pancreas (GTEx)","normal-pancreas (TCGA)"),c("normal-pancreas (GTEx)","pseudonormal-pancreas (TCGA)"),
                      c("normal-pancreas (GTEx)","tumor-pancreas (TCGA)")))
dev.off()

summary(glm(IGH_expression ~ outcome + sex + age + race, data = PAAD.GTEx.repertoire.diversity_Igreads))


#Entropy 
Ig_entropy<-melt(PAAD.GTEx.repertoire.diversity_Igclones[,c("sample","entropy_IGH","entropy_IGK","entropy_IGL","outcome")])
Ig_entropy<-Ig_entropy[which(Ig_entropy$value!=0),]
tiff("Results/boxplot_Ig_entropy_TCGA_GTEx.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                  "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  #stat_compare_means(label = "p.format")
  stat_compare_means(
    comparisons =list(c("normal-pancreas (GTEx)","normal-pancreas (TCGA)"),c("normal-pancreas (GTEx)","pseudonormal-pancreas (TCGA)"),
                      c("normal-pancreas (GTEx)","tumor-pancreas (TCGA)"),
                      c("pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")))
dev.off()

summary(glm(entropy_recon_IGH ~ outcome + sex + age + race, data = PAAD.GTEx.repertoire.diversity_Igclones))


#Entropy recon
Ig_entropy<-melt(PAAD.GTEx.repertoire.diversity_Igclones[,c("sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL","outcome")])
Ig_entropy<-Ig_entropy[which(Ig_entropy$value!=0),]
tiff("Results/boxplot_Ig_entropy_recon_TCGA_GTEx_ALL.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                  "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  #stat_compare_means(label = "p.format")
  stat_compare_means(
    comparisons =list(c("normal-pancreas (GTEx)","normal-pancreas (TCGA)"),c("normal-pancreas (GTEx)","pseudonormal-pancreas (TCGA)"),
                      c("normal-pancreas (GTEx)","tumor-pancreas (TCGA)")))
dev.off()

summary(glm(entropy_recon_IGH ~ outcome + sex + age + race, data = PAAD.GTEx.repertoire.diversity_Igclones))

#Entropy down
Ig_entropy<-melt(PAAD.GTEx.repertoire.diversity_Igclones_down[,c("sample","entropy_IGH_down","entropy_IGK_down","entropy_IGL_down","outcome")])
Ig_entropy<-Ig_entropy[which(Ig_entropy$value!=0),]
tiff("Results/boxplot_Ig_entropy_down_TCGA_GTEx_ALL.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                  "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  #stat_compare_means(label = "p.format")
  stat_compare_means(
    comparisons =list(c("normal-pancreas (GTEx)","normal-pancreas (TCGA)"),c("normal-pancreas (GTEx)","pseudonormal-pancreas (TCGA)"),
                      c("normal-pancreas (GTEx)","tumor-pancreas (TCGA)")))
dev.off()

summary(glm(entropy_recon_IGH_down ~ outcome + sex + age + race, data = PAAD.GTEx.repertoire.diversity_Igclones_down))

#T expression
T_expr<-melt(PAAD.GTEx.repertoire.diversity_tclones[,c("sample","TRA_expression_filter","TRB_expression_filter",
                                                      "TRD_expression_filter","TRG_expression_filter","outcome")])
T_expr<-T_expr[which(T_expr$value!=0),]
T_expr$value<-log10(T_expr$value)
tiff("Results/boxplot_T_expression_TCGA_GTEx.tiff",res=300,h=2500,w=3500)
ggboxplot(T_expr, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                  "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  #stat_compare_means(label = "p.format")
  stat_compare_means(
    comparisons =list(c("normal-pancreas (GTEx)","normal-pancreas (TCGA)"),c("normal-pancreas (GTEx)","pseudonormal-pancreas (TCGA)"),
                      c("normal-pancreas (GTEx)","tumor-pancreas (TCGA)")))
dev.off()

summary(glm(TRB_expression ~ outcome + sex + age + race, data = PAAD.GTEx.repertoire.diversity_treads))


#Entropy 
T_entropy<-melt(PAAD.GTEx.repertoire.diversity_tclones[,c("sample","entropy_TRA","entropy_TRB","entropy_TRD","entropy_TRG","outcome")])
T_entropy<-T_entropy[which(T_entropy$value!=0),]
tiff("Results/boxplot_T_entropy_TCGA_GTEx.tiff",res=300,h=2500,w=3500)
ggboxplot(T_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                  "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  #stat_compare_means(label = "p.format")
  stat_compare_means(
    comparisons =list(c("normal-pancreas (GTEx)","normal-pancreas (TCGA)"),c("normal-pancreas (GTEx)","pseudonormal-pancreas (TCGA)"),
                      c("normal-pancreas (GTEx)","tumor-pancreas (TCGA)")))
dev.off()


#KappaLambda_ratio_expression
KappaLambda_ratio_expression<-melt(PAAD.GTEx.repertoire.diversity_Igclones[,c("sample","KappaLambda_ratio_expression_filter","outcome")])
tiff("Results/boxplot_KappaLambda_ratio_expression_TCGA_GTEx.tiff",res=300,h=2500,w=3500)
ggboxplot(KappaLambda_ratio_expression, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                  "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  #stat_compare_means(label = "p.format")
  stat_compare_means(
    comparisons =list(c("normal-pancreas (GTEx)","normal-pancreas (TCGA)"),c("normal-pancreas (GTEx)","pseudonormal-pancreas (TCGA)"),
                      c("normal-pancreas (GTEx)","tumor-pancreas (TCGA)")))
dev.off()

###################################
### 6. Compare with GTEX Blood  ###
###################################
load("Data/PAAD_GTEX_blood_diversity.Rdata")

PAAD.GTEx.Blood.repertoire.diversity$outcome<-factor(PAAD.GTEx.Blood.repertoire.diversity$outcome)
PAAD.GTEx.Blood.repertoire.diversity$sample<-rownames(PAAD.GTEx.Blood.repertoire.diversity)

##Normalization measures
cols=brewer.pal(3,name = "Set1")

tiff("Results/boxplot_norm_expression_TCGA_GTEX_blood.tiff",res=300,h=2000,w=3500)
boxplot(t(PAAD.GTEx.Blood.repertoire.diversity[, c("IGH_expression",   "IGK_expression",   "IGL_expression", 
                                             "TRA_expression",  "TRB_expression", "TRD_expression", "TRG_expression")]),
        col=c(cols[2],cols[2],cols[1],cols[1],cols[1])[PAAD.GTEx.Blood.repertoire.diversity$outcome],main="norm_expression")
legend("topright",c("TCGA","GTEx"),fill=cols)
dev.off()
tiff("Results/boxplot_norm_expression_vst_TCGA_GTEx_blood.tiff",res=300,h=2000,w=3500)
boxplot(t(PAAD.GTEx.Blood.repertoire.diversity[, c("IGH_expression_vst",   "IGK_expression_vst",   "IGL_expression_vst", 
                                             "TRA_expression_vst",  "TRB_expression_vst", "TRD_expression_vst", "TRG_expression_vst")]),
        col=c(cols[2],cols[2],cols[1],cols[1],cols[1])[PAAD.GTEx.Blood.repertoire.diversity$outcome],main="norm_expression_vst")
legend("topright",c("TCGA","GTEx"),fill=cols)
dev.off()
tiff("Results/boxplot_norm_entropy_TCGA_GTEx_blood.tiff",res=300,h=2000,w=3500)
boxplot(t(PAAD.GTEx.Blood.repertoire.diversity[, c("entropy_IGH", "entropy_IGK", "entropy_IGL", "entropy_TRA",
                                             "entropy_TRB", "entropy_TRD", "entropy_TRG")]),
        col=c(cols[2],cols[2],cols[1],cols[1],cols[1])[PAAD.GTEx.Blood.repertoire.diversity$outcome],main="norm_entropy")
legend("topright",c("TCGA","GTEx"),fill=cols)
dev.off()

tiff("Results/boxplot_norm_entropy_recon_TCGA_GTEx_blood.tiff",res=300,h=2000,w=3500)
boxplot(t(PAAD.GTEx.Blood.repertoire.diversity[,c("entropy_recon_IGH", "entropy_recon_IGK", "entropy_recon_IGL", 
                                            "entropy_recon_TRA","entropy_recon_TRB", "entropy_recon_TRD", "entropy_recon_TRG")]),
        col=c(cols[2],cols[2],cols[1],cols[1],cols[1])[PAAD.GTEx.Blood.repertoire.diversity$outcome],main="norm_entropy_recon")
legend("topright",c("TCGA","GTEx"),fill=cols)
dev.off()


#reads
PAAD.GTEx.Blood.repertoire.diversity_treads<-PAAD.GTEx.Blood.repertoire.diversity[which(PAAD.GTEx.Blood.repertoire.diversity$T_Reads>100),]
PAAD.GTEx.Blood.repertoire.diversity_Igreads<-PAAD.GTEx.Blood.repertoire.diversity[which(PAAD.GTEx.Blood.repertoire.diversity$Ig_Reads>100),]

##Clones 
PAAD.GTEx.Blood.repertoire.diversity$T_clones<-PAAD.GTEx.Blood.repertoire.diversity$clones_TRA+PAAD.GTEx.Blood.repertoire.diversity$clones_TRB+
  PAAD.GTEx.Blood.repertoire.diversity$clones_TRD+PAAD.GTEx.Blood.repertoire.diversity$clones_TRG
PAAD.GTEx.Blood.repertoire.diversity$Ig_clones<-PAAD.GTEx.Blood.repertoire.diversity$clones_IGH+PAAD.GTEx.Blood.repertoire.diversity$clones_IGK+
  PAAD.GTEx.Blood.repertoire.diversity$clones_IGL
PAAD.GTEx.Blood.repertoire.diversity_tclones<-PAAD.GTEx.Blood.repertoire.diversity[which(PAAD.GTEx.Blood.repertoire.diversity$T_clones>100),]
PAAD.GTEx.Blood.repertoire.diversity_Igclones<-PAAD.GTEx.Blood.repertoire.diversity[which(PAAD.GTEx.Blood.repertoire.diversity$Ig_clones>100),]


brewer.pal(4,name = "Accent")
brewer.pal(4,name = "Pastel1")
cols=c( "#FBB4AE", "#7FC97F","#386CB0","#FDC086", "#BEAED4")

#Ig expression
Ig_expr<-melt(PAAD.GTEx.Blood.repertoire.diversity_Igclones[,c("sample","IGH_expression_filter","IGK_expression_filter","IGL_expression_filter","outcome")])
Ig_expr<-Ig_expr[which(Ig_expr$value!=0),]
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_TCGA_GTEX_blood.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_expr, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c( "blood (GTEx)","normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                  "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  #stat_compare_means(label = "p.format")
  stat_compare_means(
  comparisons =list(c( "blood (GTEx)","normal-pancreas (GTEx)"),c("blood (GTEx)","normal-pancreas (TCGA)"),
                    c("blood (GTEx)","pseudonormal-pancreas (TCGA)"),c("blood (GTEx)","tumor-pancreas (TCGA)")))

dev.off()

#T expression
T_expr<-melt(PAAD.GTEx.Blood.repertoire.diversity_treads[,c("sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression","outcome")])
T_expr<-T_expr[which(T_expr$value!=0),]
T_expr$value<-log10(T_expr$value)
tiff("Results/boxplot_T_expression_TCGA_GTEx_blood.tiff",res=300,h=2500,w=3500)
ggboxplot(T_expr, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c( "blood (GTEx)","normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                   "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  #stat_compare_means(label = "p.format")
  stat_compare_means(
    comparisons =list(c( "blood (GTEx)","normal-pancreas (GTEx)"),c("blood (GTEx)","normal-pancreas (TCGA)"),
                      c("blood (GTEx)","pseudonormal-pancreas (TCGA)"),c("blood (GTEx)","tumor-pancreas (TCGA)")))
dev.off()

#Entropy
Ig_entropy<-melt(PAAD.GTEx.Blood.repertoire.diversity_Igclones[,c("sample","entropy_IGH","entropy_IGK","entropy_IGL","outcome")])
Ig_entropy<-Ig_entropy[which(Ig_entropy$value!=0),]
tiff("Results/boxplot_Ig_entropy_TCGA_GTEx_blood.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c( "blood (GTEx)","normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                   "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  #stat_compare_means(label = "p.format")
  stat_compare_means(
    comparisons =list(c( "blood (GTEx)","normal-pancreas (GTEx)"),c("blood (GTEx)","normal-pancreas (TCGA)"),
                      c("blood (GTEx)","pseudonormal-pancreas (TCGA)"),c("blood (GTEx)","tumor-pancreas (TCGA)")))
dev.off()

#Clones
Ig_clones<-melt(PAAD.GTEx.Blood.repertoire.diversity_Igclones[,c("sample","clones_IGH","clones_IGK","clones_IGL","outcome")])
Ig_clones<-Ig_clones[which(Ig_clones$value!=0),]
tiff("Results/boxplot_Ig_clones_TCGA_GTEx_blood.tiff",res=300,h=2500,w=3500)
ggboxplot(Ig_clones, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c( "blood (GTEx)","normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                   "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  #stat_compare_means(label = "p.format")
  stat_compare_means(
    comparisons =list(c( "blood (GTEx)","normal-pancreas (GTEx)"),c("blood (GTEx)","normal-pancreas (TCGA)"),
                      c("blood (GTEx)","pseudonormal-pancreas (TCGA)"),c("blood (GTEx)","tumor-pancreas (TCGA)")))
dev.off()


#Entropy recon
T_entropy<-melt(PAAD.GTEx.Blood.repertoire.diversity_Igreads[,c("sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG","outcome")])
T_entropy<-T_entropy[which(T_entropy$value!=0),]
tiff("Results/boxplot_T_entropy_recon_TCGA_GTEx_blood.tiff",res=300,h=2500,w=3500)
ggboxplot(T_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value,color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c( "blood (GTEx)","normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                                   "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)")) +
  #stat_compare_means(label = "p.format")
  stat_compare_means(
    comparisons =list(c( "blood (GTEx)","normal-pancreas (GTEx)"),c("blood (GTEx)","normal-pancreas (TCGA)"),
                      c("blood (GTEx)","pseudonormal-pancreas (TCGA)"),c("blood (GTEx)","tumor-pancreas (TCGA)")))
dev.off()
