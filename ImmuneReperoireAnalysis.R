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
brewer.pal(3,name = "Accent")
cols=c( "#7FC97F", "#FDC086", "#BEAED4" )

#################
### 1. TCGA #####
#################

###T markers
#Barplot
tiff("Results/barplot_Treads.tiff",res=300,h=2500,w=4000)
barplot(PAAD.repertoire.diversity$T_Reads,col=cols[PAAD.repertoire.diversity$Tumor_type_3categ],main="Number of T-Reads",xlab = "Samples", ylab = "Reads",las=2)
abline(h=100)
legend("topright", legend=levels(PAAD.repertoire.diversity$Tumor_type_3categ),col=cols,pch=15, cex=0.8)
dev.off()

tiff("Results/Corr_plot_reads.tiff",res=300,h=2500,w=4000)
plot(PAAD.repertoire.diversity$Total_Reads,PAAD.repertoire.diversity$totalSeqReads,col=cols[PAAD.repertoire.diversity$Tumor_type_3categ],pch=19,
     xlab = "Total B- and T- reads aligned", ylab = "Total sequencing reads",
     main = paste0("rho = ",round(cor(PAAD.repertoire.diversity$Total_Reads,PAAD.repertoire.diversity$totalSeqReads),2)))
legend("bottomright", legend=levels(PAAD.repertoire.diversity$Tumor_type_3categ),col=cols,pch=15, cex=0.8)
dev.off()

##Differences by tumor categ 
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
#Barplot
tiff("Results/barplot_Igreads.tiff",res=300,h=2500,w=4000)
barplot(PAAD.repertoire.diversity$IG_Reads,col=cols[PAAD.repertoire.diversity$Tumor_type_3categ],main="Number of Ig-Reads",xlab = "Samples", ylab = "Reads",las=2)
abline(h=100)
legend("topright", legend=levels(PAAD.repertoire.diversity$Tumor_type_3categ),col=cols,pch=15, cex=0.8)
dev.off()

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
Ig_expr<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","IGH_expression","IGK_expression","IGL_expression","Tumor_type_3categ")])
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_PAAD.tiff",res=300,h=2500,w=3000)
ggplot(Ig_expr) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_3categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = cols, labels = c("normal_pancreas","pseudonormal_pancreas", "tumor_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_3categ)) 
dev.off()

TR_expr<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression","Tumor_type_3categ")])
TR_expr<-TR_expr[which(TR_expr$value!=0),]
TR_expr$value<-log10(TR_expr$value)
tiff("Results/boxplot_TR_expression_PAAD.tiff",res=300,h=2500,w=3500)
ggplot(TR_expr) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_3categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = cols, labels = c("normal_pancreas","pseudonormal_pancreas", "tumor_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_3categ)) 
dev.off()

Ig_entropy<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL","Tumor_type_3categ")])
tiff("Results/boxplot_Ig_entropy_PAAD.tiff",res=300,h=2500,w=3000)
ggplot(Ig_entropy) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_3categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = cols, labels = c("normal_pancreas","pseudonormal_pancreas", "tumor_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_3categ)) 
dev.off()

TR_entropy<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG","Tumor_type_3categ")])
TR_entropy<-TR_entropy[which(TR_entropy$value!=0),]
tiff("Results/boxplot_TR_entropy_PAAD.tiff",res=300,h=2500,w=3500)
ggplot(TR_entropy) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_3categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = cols, labels = c("normal_pancreas","pseudonormal_pancreas", "tumor_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_3categ))
dev.off()

kappa_lambda<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","KappaLambda_ratio_expression","Tumor_type_3categ")])
tiff("Results/boxplot_kappa_lambda_PAAD.tiff",res=300,h=2500,w=3000)
ggplot(kappa_lambda) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_3categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = cols, labels =  c("normal_pancreas","pseudonormal_pancreas", "tumor_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_3categ)) 
dev.off()

alpha_beta_ratio<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","Alpha_Beta_ratio_expression","Tumor_type_3categ")])
tiff("Results/boxplot_alpha_beta_PAAD.tiff",res=300,h=2500,w=3500)
ggplot(alpha_beta_ratio) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_3categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = cols, labels = c("normal_pancreas","pseudonormal_pancreas", "tumor_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_3categ)) 
dev.off()

Ig_cdr3length<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","cdr3_length_IGH","cdr3_length_IGK","cdr3_length_IGL","Tumor_type_3categ")])
tiff("Results/boxplot_Ig_cdr3_length.tiff",res=300,h=2500,w=3000)
ggplot(Ig_cdr3length) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_3categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = cols, labels = c("normal_pancreas","pseudonormal_pancreas", "tumor_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_3categ)) 
dev.off()

TCR_cdr3length<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","cdr3_length_TRA","cdr3_length_TRB","cdr3_length_TRD","cdr3_length_TRG","Tumor_type_3categ")])
TCR_cdr3length<-TCR_cdr3length[which(TCR_cdr3length$value!=0),]
tiff("Results/boxplot_TCR_cdr3_length.tiff",res=300,h=2500,w=3000)
ggplot(TCR_cdr3length) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_3categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_3categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = cols, labels = c("normal_pancreas","pseudonormal_pancreas", "tumor_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_3categ)) 
dev.off()

##########################
### 2. Compare with GTEX ###
#########################
load("Data/GTEx/Pancreas/GTEx_FullData.Rdata")

####Summary plots
#IgExpression
Ig_expr_PAAD<-melt(PAAD.repertoire.diversity[,c("TCGA_sample","IGH_expression","IGK_expression","IGL_expression","Tumor_type_2categ")])
Pancreas.repertoire.diversity$TCGA_sample<-rownames(Pancreas.repertoire.diversity)
Pancreas.repertoire.diversity$Tumor_type_2categ<-rep("Normal_pancreas")
Ig_expr_Pancreas<-melt(Pancreas.repertoire.diversity[,c("TCGA_sample","IGH_expression","IGK_expression","IGL_expression","Tumor_type_2categ")])

Ig_expr<-rbind(Ig_expr_PAAD,Ig_expr_Pancreas)
Ig_expr<-Ig_expr[which(Ig_expr$value!=0),]
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression.tiff",res=300,h=2500,w=3500)
ggplot(Ig_expr) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2],cols[3]), labels = c("adjacent_normal_pseudonormal_tumor_pancreas", "Tumor_pancreas","Normal_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()

#Entropy
Ig_entropy_PAAD<-melt(PAAD.repertoire.diversity[,c("TCGA_sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL","Tumor_type_2categ")])
Ig_entropy_Pancreas<-melt(Pancreas.repertoire.diversity[,c("TCGA_sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL","Tumor_type_2categ")])

Ig_expr<-rbind(Ig_entropy_PAAD,Ig_entropy_Pancreas)
Ig_expr<-Ig_expr[which(Ig_expr$value!=0),]
tiff("Results/boxplot_Ig_entropy.tiff",res=300,h=2500,w=3500)
ggplot(Ig_expr) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2],cols[3]), labels = c("adjacent_normal_pseudonormal_tumor_pancreas", "Tumor_pancreas","Normal_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()

#Texpression
T_expr_PAAD<-melt(PAAD.repertoire.diversity[,c("TCGA_sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression","Tumor_type_2categ")])
T_expr_Pancreas<-melt(Pancreas.repertoire.diversity[,c("TCGA_sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression","Tumor_type_2categ")])

T_expr<-rbind(T_expr_PAAD,T_expr_Pancreas)
T_expr<-T_expr[which(T_expr$value!=0),]
T_expr$value<-log10(T_expr$value)
tiff("Results/boxplot_T_expression.tiff",res=300,h=2500,w=3500)
ggplot(T_expr) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2],cols[3]), labels = c("adjacent_normal_pseudonormal_tumor_pancreas", "Tumor_pancreas","Normal_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()


#Entropy (T)
T_entropy_PAAD<-melt(PAAD.repertoire.diversity[,c("TCGA_sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG","Tumor_type_2categ")])
T_entropy_Pancreas<-melt(Pancreas.repertoire.diversity[,c("TCGA_sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG","Tumor_type_2categ")])

T_expr<-rbind(T_entropy_PAAD,T_entropy_Pancreas)
T_expr<-T_expr[which(T_expr$value!=0),]
tiff("Results/boxplot_T_entropy.tiff",res=300,h=2500,w=3500)
ggplot(T_expr) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2],cols[3]), labels = c("adjacent_normal_pseudonormal_tumor_pancreas", "Tumor_pancreas","Normal_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()

#KappaLambda
kappa_lambda_PAAD<-melt(PAAD.repertoire.diversity[,c("TCGA_sample","KappaLambda_ratio_expression","Tumor_type_2categ")])
kappa_lambda_Pancreas<-melt(Pancreas.repertoire.diversity[,c("TCGA_sample","KappaLambda_ratio_expression","Tumor_type_2categ")])

kappa_lambda<-rbind(kappa_lambda_PAAD,kappa_lambda_Pancreas)
kappa_lambda<-kappa_lambda[which(kappa_lambda$value!=0),]
tiff("Results/boxplot_kappa_lambda.tiff",res=300,h=2500,w=3500)
ggplot(kappa_lambda) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2],cols[3]), labels = c("adjacent_normal_pseudonormal_tumor_pancreas", "Tumor_pancreas","Normal_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()


#alphabeta
alpha_beta_ratio_PAAD<-melt(PAAD.repertoire.diversity[,c("TCGA_sample","Alpha_Beta_ratio_expression","Tumor_type_2categ")])
alpha_beta_ratio_Pancreas<-melt(Pancreas.repertoire.diversity[,c("TCGA_sample","Alpha_Beta_ratio_expression","Tumor_type_2categ")])

alpha_beta_ratio<-rbind(alpha_beta_ratio_PAAD,alpha_beta_ratio_Pancreas)
alpha_beta_ratio<-alpha_beta_ratio[which(alpha_beta_ratio$value!=1),]
tiff("Results/boxplot_alpha_beta.tiff",res=300,h=2500,w=3500)
ggplot(alpha_beta_ratio) + geom_boxplot(aes(x=variable, y=value, color=Tumor_type_2categ), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=Tumor_type_2categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2],cols[3]), labels = c("adjacent_normal_pseudonormal_tumor_pancreas", "Tumor_pancreas","Normal_pancreas")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=Tumor_type_2categ)) 
dev.off()

############################
### Validation Pancreas ####
############################
load("Data/Pancreas_Validation/Pancreas_Validation_RepertoireResults_diversity.Rdata")
cols=brewer.pal(5,name = "Set1")

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

Pancreas.Validation.repertoire.diversity_treads<-Pancreas.Validation.repertoire.diversity[which(Pancreas.Validation.repertoire.diversity$T_Reads>100),]

##Ig markers
#Barplot
Pancreas.Validation.repertoire.diversity$IG_Reads<-Pancreas.Validation.repertoire.diversity$IGH+Pancreas.Validation.repertoire.diversity$IGK+
  Pancreas.Validation.repertoire.diversity$IGL

tiff("Results/barplot_Igreads_PancreasValidation.tiff",res=300,h=2500,w=4000)
barplot(Pancreas.Validation.repertoire.diversity$IG_Reads,col=cols[Pancreas.Validation.repertoire.diversity$tissue],main="Number of Ig-Reads",xlab = "Samples", ylab = "Reads",las=2)
abline(h=100)
legend("topright", legend=levels(Pancreas.Validation.repertoire.diversity$tissue),col=cols,pch=15, cex=0.8)
dev.off()

Pancreas.Validation.repertoire.diversity_Igreads<-Pancreas.Validation.repertoire.diversity[which(Pancreas.Validation.repertoire.diversity$IG_Reads>100),]

####Summary plots
#IgExpression
Ig_expr1<-melt(Pancreas.Validation.repertoire.diversity_Igreads[,c("sample","IGH_expression","IGK_expression","IGL_expression","tissue")])
Ig_expr2<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","IGH_expression","IGK_expression","IGL_expression","Tumor_type_3categ")])
colnames(Ig_expr2)[1:2]<-c("sample","tissue")

Ig_expr<-rbind(Ig_expr1,Ig_expr2)
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/boxplot_Ig_expression_All.tiff",res=300,h=2500,w=3500)
ggplot(Ig_expr) + geom_boxplot(aes(x=variable, y=value, color=tissue), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=tissue), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[1],cols[2],cols[2],cols[2]), labels = c("normal pancreas (val)", "pancreas tumor (val)","normal pancreas (TCGA)",
                                               "pseudonormal (TCGA)","pancreas tumor (TCGA)")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=tissue)) 
dev.off()

#Entropy
Ig_entropy1<-melt(Pancreas.Validation.repertoire.diversity_Igreads[,c("sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL","tissue")])
Ig_entropy2<-melt(PAAD.repertoire.diversity_Igreads[,c("TCGA_sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL","Tumor_type_3categ")])
colnames(Ig_entropy2)[1:2]<-c("sample","tissue")

Ig_entropy<-rbind(Ig_entropy1,Ig_entropy2)
tiff("Results/boxplot_Ig_entropy_All.tiff",res=300,h=2500,w=3500)
ggplot(Ig_entropy) + geom_boxplot(aes(x=variable, y=value, color=tissue), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=tissue), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[1],cols[2],cols[2],cols[2]), labels = c("normal pancreas (val)", "pancreas tumor (val)","normal pancreas (TCGA)",
                                                                                      "pseudonormal (TCGA)","pancreas tumor (TCGA)")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=tissue)) 
dev.off()

#Texpression
T_expr1<-melt(Pancreas.Validation.repertoire.diversity_treads[,c("sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression","tissue")])
T_expr2<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression","Tumor_type_3categ")])
colnames(T_expr2)[1:2]<-c("sample","tissue")

T_expr<-rbind(T_expr1,T_expr2)
T_expr<-T_expr[which(T_expr$value!=0),]
T_expr$value<-log10(T_expr$value)
tiff("Results/boxplot_T_expression_All.tiff",res=300,h=2500,w=3500)
ggplot(T_expr) + geom_boxplot(aes(x=variable, y=value, color=tissue), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=tissue), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[1],cols[2],cols[2],cols[2]), labels = c("normal pancreas (val)", "pancreas tumor (val)","normal pancreas (TCGA)",
                                                                                      "pseudonormal (TCGA)","pancreas tumor (TCGA)")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=tissue)) 
dev.off()


#Entropy (T)
T_entropy1<-melt(Pancreas.Validation.repertoire.diversity_treads[,c("sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG","tissue")])
T_entropy2<-melt(PAAD.repertoire.diversity_treads[,c("TCGA_sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG","Tumor_type_3categ")])
colnames(T_entropy2)[1:2]<-c("sample","tissue")

T_entropy<-rbind(T_entropy1,T_entropy2)
T_entropy<-T_entropy[which(T_entropy$value!=0),]
tiff("Results/boxplot_T_entropy_All.tiff",res=300,h=2500,w=3500)
ggplot(T_entropy) + geom_boxplot(aes(x=variable, y=value, color=tissue), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=tissue), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[1],cols[2],cols[2],cols[2]), labels = c("normal pancreas (val)", "pancreas tumor (val)","normal pancreas (TCGA)",
                                                                                      "pseudonormal (TCGA)","pancreas tumor (TCGA)")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=tissue)) 
dev.off()

#KappaLambda
kappa_lambda<-melt(Pancreas.Validation.repertoire.diversity_Igreads[,c("sample","KappaLambda_ratio_expression","tissue")])

tiff("Results/boxplot_kappa_lambda_PancreasValidation.tiff",res=300,h=2500,w=3500)
ggplot(kappa_lambda) + geom_boxplot(aes(x=variable, y=value, color=tissue), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=tissue), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal pancreas", "pancreas tumor")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=tissue)) 
dev.off()


#alphabeta
alpha_beta_ratio<-melt(Pancreas.Validation.repertoire.diversity_treads[,c("sample","Alpha_Beta_ratio_expression","tissue")])

tiff("Results/boxplot_alpha_beta_PancreasValidation.tiff",res=300,h=2500,w=3500)
ggplot(alpha_beta_ratio) + geom_boxplot(aes(x=variable, y=value, color=tissue), alpha = 0, position = position_dodge(width = .8)) + 
  geom_point(aes(x=variable, y=value, color=tissue), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols[1], cols[2]), labels = c("normal pancreas", "pancreas tumor")) +
  theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+
  stat_compare_means(aes(x=variable, y=value, color=tissue)) 
dev.off()






############################
### GTEX Blood           ###
############################
load("Data/GTEx/Blood/MIXCR/GTEX_FullData.Rdata")
totalReads<-read.table("Data/GTEx/Blood/MIXCR/total_reads_GTEX.txt",sep=";")
id<-match(rownames(GTEX.repertoire.diversity),totalReads$V1)
GTEX.repertoire.diversity$TotalSeq<-totalReads[id,"V2"]


##Me he quedado aqui
sra<-read.csv("Data/GTEx/Blood/SraRunTableBloodRNAseq.csv")
rownames(GTEX.repertoire.diversity)<-unlist(strsplit(as.character(rownames(GTEX.repertoire.diversity)), "\\}"))
id<-match(sra$Run,rownames(GTEX.repertoire.diversity))
GTEX.repertoire.diversity<-GTEX.repertoire.diversity[na.omit(id),] ##439
id<-match(rownames(GTEX.repertoire.diversity),sra$Run)
sra[na.omit(id),""]
sra$Run[id]
cols=brewer.pal(5,name = "Set1")

