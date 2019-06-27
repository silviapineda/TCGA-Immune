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

load("Data/PAAD/PAAD_FullData.Rdata")
##################
####Descriptive analysis to see if there are differences by tumor and adjacent_normal
##################
cols=brewer.pal(3,name = "Accent")

###T markers
PAAD.repertoire.diversity_treads<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$T_Reads>100),]
T_markers<-c("TRA_expression","TRB_expression","TRD_expression","TRG_expression","Alpha_Beta_ratio_expression",
             "clones_recon_TRA","clones_recon_TRB","clones_recon_TRD","clones_recon_TRG","entropy_recon_TRA",
             "entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG")

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

p.T_markers=NULL
for(i in 1:length(T_markers)){
  p.T_markers[i]<-coef(summary(glm(PAAD.repertoire.diversity_treads[,T_markers[i]]~PAAD.repertoire.diversity_treads$Tumor_type_2categ)))[2,4]
  PAAD.repertoire.diversity_treads$Marker<-PAAD.repertoire.diversity_treads[,T_markers[i]]
  tiff(paste0("Results/boxplot_",T_markers[i],".tiff"),res=300,h=2500,w=3000)
  print(ggplot(PAAD.repertoire.diversity_treads) + 
          geom_boxplot(aes(y=Marker, color=Tumor_type_2categ, x=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) +
          geom_point(aes(y=Marker, color=Tumor_type_2categ,x=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +     
          scale_y_continuous(name=T_markers[i]) + stat_compare_means(aes(x=Tumor_type_2categ, y=Marker, color=Tumor_type_2categ),label.x.npc = "center") +
          geom_text(aes(y=Marker,x=Tumor_type_2categ,label=ifelse(is_outlier(PAAD.repertoire.diversity_treads$Marker)==T,
                                                                  as.character(PAAD.repertoire.diversity_treads$TCGA_sample),""),hjust=-0.1),size = 3) +
          scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal_pseudonormal_pancreas", "Tumor_pancres")) )
  
  dev.off()
  
}
T_markers[which(p.T_markers<0.05)] ## ""TRD_expression"    "TRG_expression"    "clones_recon_TRG"  "entropy_recon_TRG"

##Ig markers
PAAD.repertoire.diversity_Igreads<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$IG_Reads>100),]
Ig_markers<-c("IGH_expression","IGK_expression","IGL_expression","KappaLambda_ratio_expression",
              "clones_recon_IGH","clones_recon_IGK","clones_recon_IGL","entropy_recon_IGH",
              "entropy_recon_IGK","entropy_recon_IGL")

p.Ig_markers=NULL
for(i in 1:length(Ig_markers)){
  p.Ig_markers[i]<-coef(summary(glm(PAAD.repertoire.diversity_Igreads[,Ig_markers[i]]~PAAD.repertoire.diversity_Igreads$Tumor_type_2categ)))[2,4]
  PAAD.repertoire.diversity_Igreads$Marker<-PAAD.repertoire.diversity_Igreads[,Ig_markers[i]]
  is_outlier(PAAD.repertoire.diversity_Igreads$Marker)
  tiff(paste0("Results/boxplot_",Ig_markers[i],".tiff"),res=300,h=2500,w=3000)
  print(ggplot(PAAD.repertoire.diversity_Igreads) + 
          geom_boxplot(aes(y=Marker, color=Tumor_type_2categ, x=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) +
          geom_point(aes(y=Marker, color=Tumor_type_2categ,x=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +     
          scale_y_continuous(name=Ig_markers[i]) + stat_compare_means(aes(x=Tumor_type_2categ, y=Marker, color=Tumor_type_2categ),label.x.npc = "center") +
          geom_text(aes(y=Marker,x=Tumor_type_2categ,label=ifelse(is_outlier(PAAD.repertoire.diversity_Igreads$Marker)==T,
                                                                  as.character(PAAD.repertoire.diversity_Igreads$TCGA_sample),""),hjust=-0.1),size = 3) +
          scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal_pseudonormal_pancreas", "Tumor_pancres")) )
  
  dev.off()
}
Ig_markers[which(p.Ig_markers<0.05)] ## "clones_recon_IGH"  "clones_recon_IGL"  "entropy_recon_IGH"

####Summary plots
Ig_expr<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","IGH_expression","IGK_expression","IGL_expression","Tumor_type_2categ")])
Ig_expr<-Ig_expr[which(Ig_expr$value!=0),]
tiff("Results/boxplot_Ig_expression.tiff",res=300,h=2500,w=3000)
ggplot(Ig_expr) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal_pseudonormal_pancreas", "Tumor_pancres")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()

TR_expr<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression","Tumor_type_2categ")])
TR_expr<-TR_expr[which(TR_expr$value!=0),]
tiff("Results/boxplot_TR_expression.tiff",res=300,h=2500,w=3500)
ggplot(TR_expr) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal_pseudonormal_pancreas", "Tumor_pancres")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()

Ig_clones<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","clones_recon_IGH","clones_recon_IGK","clones_recon_IGL","Tumor_type_2categ")])
Ig_clones<-Ig_clones[which(Ig_clones$value!=0),]
tiff("Results/boxplot_Ig_clones.tiff",res=300,h=2500,w=3000)
ggplot(Ig_clones) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal_pseudonormal_pancreas", "Tumor_pancres")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()

TR_clones<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","clones_recon_TRA","clones_recon_TRB","clones_recon_TRD","clones_recon_TRG","Tumor_type_2categ")])
TR_clones<-TR_clones[which(TR_clones$value!=0),]
tiff("Results/boxplot_TR_clones.tiff",res=300,h=2500,w=3500)
ggplot(TR_clones) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal_pseudonormal_pancreas", "Tumor_pancres")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()

Ig_entropy<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL","Tumor_type_2categ")])
Ig_entropy<-Ig_entropy[which(Ig_entropy$value!=0),]
tiff("Results/boxplot_Ig_entropy.tiff",res=300,h=2500,w=3000)
ggplot(Ig_entropy) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal_pseudonormal_pancreas", "Tumor_pancres")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()

TR_entropy<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG","Tumor_type_2categ")])
TR_entropy<-TR_entropy[which(TR_entropy$value!=0),]
tiff("Results/boxplot_TR_entropy.tiff",res=300,h=2500,w=3500)
ggplot(TR_entropy) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal_pseudonormal_pancreas", "Tumor_pancres")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ))
dev.off()

kappa_lambda<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","KappaLambda_ratio_expression","Tumor_type_2categ")])
tiff("Results/boxplot_kappa_lambda.tiff",res=300,h=2500,w=3000)
ggplot(kappa_lambda) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal_pseudonormal_pancreas", "Tumor_pancres")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()

alpha_beta_ratio<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","Alpha_Beta_ratio_expression","Tumor_type_2categ")])
alpha_beta_ratio<-alpha_beta_ratio[which(alpha_beta_ratio$value!=1),]
tiff("Results/boxplot_alpha_beta.tiff",res=300,h=2500,w=3500)
ggplot(alpha_beta_ratio) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal_pseudonormal_pancreas", "Tumor_pancres")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()

Ig_cdr3length<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","cdr3_length_IGH","cdr3_length_IGK","cdr3_length_IGL","Tumor_type_2categ")])
tiff("Results/boxplot_Ig_cdr3_length.tiff",res=300,h=2500,w=3000)
ggplot(Ig_cdr3length) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal_pseudonormal_pancreas", "Tumor_pancres")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()

TCR_cdr3length<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","cdr3_length_TRA","cdr3_length_TRB","cdr3_length_TRD","cdr3_length_TRG","Tumor_type_2categ")])
TCR_cdr3length<-TCR_cdr3length[which(TCR_cdr3length$value!=0),]
tiff("Results/boxplot_TCR_cdr3_length.tiff",res=300,h=2500,w=3000)
ggplot(TCR_cdr3length) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal_pseudonormal_pancreas", "Tumor_pancres")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()
