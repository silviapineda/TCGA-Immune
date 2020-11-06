rm(list = ls(all = TRUE))
x <-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy outcomes
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Network analysis
###         
###
### Author: Silvia Pineda
### Date: January, 2019
############################################################################################
library(ggplot2)
library("igraph")
library("ineq")
library("dplyr")
library("RColorBrewer")
library("ggpubr")

setwd("~/TCGA-Immune/")
load("Data/PAAD/PAAD_FullData.Rdata")
load("Data/GTEx/Pancreas/GTEx_FullData.Rdata")
load("Data/Pancreas_Validation/Pancreas_Validation_FullData.Rdata")
load("Data/Validation_Normal_pancreas/Pancreas_Normal_Validation_FullData.Rdata")
#load("Data/GTEx/Blood/GTEx_FullData_OnlyDiversity.Rdata")

########################
##4.Plot the network
########################
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")

#### TCGA - Only PDAC samples
PAAD.repertoire.tumor<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$Tumor_type_4categ=="PDAC"),] #144
PAAD.repertoire.tumor$TCGA_sample<-substr(PAAD.repertoire.tumor$TCGA_sample,1,15)
#### GTEX - Pancreas
GTEX.repertoire.normal<-Pancreas.repertoire.diversity
##Tumor Validation
Validation.repertoire.tumor<-Pancreas.Validation.repertoire.diversity[which(Pancreas.Validation.repertoire.diversity$tissue=="pancreas tumor"),]
##Normal Validation
Validation.repertoire.normal<-Pancreas.Normal.Validation.repertoire.diversity

Cluster_vertex_gini_distribution<-rbind(PAAD.repertoire.tumor[,c("cluster_gini_IGH","vertex_gini_IGH","cluster_gini_IGK","vertex_gini_IGK","cluster_gini_IGL","vertex_gini_IGL",
                                                                 "cluster_gini_TRA","vertex_gini_TRA","cluster_gini_TRB","vertex_gini_TRB")],
                                        GTEX.repertoire.normal[,c("cluster_gini_IGH","vertex_gini_IGH","cluster_gini_IGK","vertex_gini_IGK","cluster_gini_IGL","vertex_gini_IGL",
                                                                  "cluster_gini_TRA","vertex_gini_TRA","cluster_gini_TRB","vertex_gini_TRB")],
                                        Validation.repertoire.tumor[,c("cluster_gini_IGH","vertex_gini_IGH","cluster_gini_IGK","vertex_gini_IGK","cluster_gini_IGL","vertex_gini_IGL",
                                                                       "cluster_gini_TRA","vertex_gini_TRA","cluster_gini_TRB","vertex_gini_TRB")],
                                        Validation.repertoire.normal[,c("cluster_gini_IGH","vertex_gini_IGH","cluster_gini_IGK","vertex_gini_IGK","cluster_gini_IGL","vertex_gini_IGL",
                                                                        "cluster_gini_TRA","vertex_gini_TRA","cluster_gini_TRB","vertex_gini_TRB")])
Cluster_vertex_gini_distribution$outcome<-c(rep("TCGA-PDAC",nrow(PAAD.repertoire.tumor)),rep("GTEX-Normal",nrow(GTEX.repertoire.normal)),
                                           rep("Validation-PDAC",nrow(Validation.repertoire.tumor)),rep("Validation-Normal",nrow(Validation.repertoire.normal)))
Cluster_vertex_gini_distribution$outcome<-factor(Cluster_vertex_gini_distribution$outcome)

#Cluster_vertex_gini_distribution<-Cluster_vertex_gini_distribution[which(Cluster_vertex_gini_distribution$cluster_gini_TRB!=0),]
chainType="TRB"
tiff(paste0("Results/network_vertex_cluster_gini_",chainType,".tiff"),h=1700,w=1700,res=300)
par(fig=c(0,0.8,0,0.8))
plot(Cluster_vertex_gini_distribution[,paste0("cluster_gini_",chainType)], 
     Cluster_vertex_gini_distribution[,paste0("vertex_gini_",chainType)],
     col = cols[factor(Cluster_vertex_gini_distribution$outcome)],
     pch=20,ylab = c("Gini (Vextex)"),xlab = c("Gini (Cluster)"),ylim=c(0,0.9),xlim=c(0,0.8))
legend("bottomright",legend=c("GTEX-Normal",
                              "TCGA-PDAC",
                              "Validation-Normal",
                              "Validation-PDAC"), 
       col=cols,pch=20,cex=0.8)

par(fig=c(0,0.8,0.5,1), new=TRUE)
#summary(glm(Cluster_vertex_gini_distribution[,paste0("cluster_gini_",chainType)]~Cluster_vertex_gini_distribution$outcome))
boxplot(Cluster_vertex_gini_distribution[,paste0("cluster_gini_",chainType)]~Cluster_vertex_gini_distribution$outcome,
        col=cols, horizontal=TRUE, axes=FALSE,ann=FALSE,ylim=c(0,0.8))



par(fig=c(0.6,1,0,0.8),new=TRUE)
#summary(glm(Cluster_vertex_gini_distribution[,paste0("vertex_gini_",chainType)]~Cluster_vertex_gini_distribution$outcome))
boxplot(Cluster_vertex_gini_distribution[,paste0("vertex_gini_",chainType)]~Cluster_vertex_gini_distribution$outcome,
        col=cols, axes=FALSE,ann=FALSE,ylim=c(0,0.9))

dev.off()


tiff(paste0("Results/network_cluster_gini_",chainType,"_ALL.tiff"),h=2000,w=2000,res=300)
ggboxplot(Cluster_vertex_gini_distribution, x = "outcome" , y =  "cluster_gini_TRB",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y= cluster_gini_TRB, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Normal",
                                                  "TCGA-PDAC",
                                                  "Validation-Normal",
                                                  "Validation-PDAC")) +
  
  stat_compare_means(
    comparisons =list(c("GTEX-Normal","Validation-Normal"),c("GTEX-Normal","TCGA-PDAC"),
                      c("TCGA-PDAC","Validation-PDAC"),c("Validation-Normal","Validation-PDAC"),
                      c("TCGA-PDAC","Validation-Normal")))

dev.off()

tiff(paste0("Results/network_vertex_gini_",chainType,"_ALL.tiff"),h=2000,w=2000,res=300)
ggboxplot(Cluster_vertex_gini_distribution, x = "outcome" , y =  "vertex_gini_TRB",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y= vertex_gini_TRB, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Normal",
                                                  "TCGA-PDAC",
                                                  "Validation-Normal",
                                                  "Validation-PDAC")) +
  
  stat_compare_means(
    comparisons =list(c("GTEX-Normal","Validation-Normal"),c("GTEX-Normal","TCGA-PDAC"),
                      c("TCGA-PDAC","Validation-PDAC"),c("Validation-Normal","Validation-PDAC"),
                      c("TCGA-PDAC","Validation-Normal")))
dev.off()

###TCR
tiff(paste0("Results/network_vertex_gini_TRA.tiff"),h=1700,w=1700,res=300)
ggboxplot(Cluster_vertex_gini_distribution, x = "outcome" , y =  "cluster_gini_TRA",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y= cluster_gini_TRA, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Normal",
                                                            "TCGA-PDAC",
                                                            "Validation-Normal",
                                                            "Validation-PDAC")) +
  
  stat_compare_means(
    comparisons =list(c("GTEX-Normal","TCGA-PDAC"),
                      c("TCGA-PDAC","Validation-PDAC")))
dev.off()


######################################
##4.Plot the network with GTEX blood
######################################
#### TCGA - 
PAAD.repertoire<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type_4categ!="PAC-Other"),] #131
PAAD.repertoire$Tumor_type_4categ<-factor(PAAD.repertoire$Tumor_type_4categ)
PAAD.repertoire$TCGA_sample<-substr(PAAD.repertoire$TCGA_sample,1,15)
#### GTEX - Pancreas
GTEX.repertoire.normal<-Pancreas.repertoire.diversity
##GTEx - Blood
GTEX.blood.repertoire.diversity<-GTEX.blood.repertoire.diversity[which(is.na(GTEX.blood.repertoire.diversity$totalReads)==F),]

Cluster_vertex_gini_distribution<-rbind(PAAD.repertoire[,c("cluster_gini_IGK","vertex_gini_IGK")],
                                        GTEX.repertoire.normal[,c("cluster_gini_IGK","vertex_gini_IGK")],
                                        GTEX.blood.repertoire.diversity[,c("cluster_gini_IGK","vertex_gini_IGK")])

Cluster_vertex_gini_distribution$outcome<-c(as.character(PAAD.repertoire$Tumor_type_4categ),rep("GTEX-Normal",nrow(GTEX.repertoire.normal)),
                                rep("GTEX-Blood",nrow(GTEX.blood.repertoire.diversity)))
Cluster_vertex_gini_distribution$outcome<-replace(Cluster_vertex_gini_distribution$outcome,Cluster_vertex_gini_distribution$outcome=="normal_pancreas","TCGA-normal-adj-pancreas")
Cluster_vertex_gini_distribution$outcome<-replace(Cluster_vertex_gini_distribution$outcome,Cluster_vertex_gini_distribution$outcome=="PDAC","TCGA-PDAC")
Cluster_vertex_gini_distribution$outcome<-replace(Cluster_vertex_gini_distribution$outcome,Cluster_vertex_gini_distribution$outcome=="pseudonormal_pancreas","TCGA-pseudonormal_pancreas")

Cluster_vertex_gini_distribution$outcome<-factor(Cluster_vertex_gini_distribution$outcome)


chainType="IGK"
cols= c("#FBB4AE","#7FC97F","#FDC086","#BEAED4","#B3CDE3")
tiff(paste0("Results/network_vertex_cluster_gini_",chainType,"_ALL_blood.tiff"),h=2000,w=2000,res=300)
par(fig=c(0,0.8,0,0.8))
plot(Cluster_vertex_gini_distribution[,paste0("cluster_gini_",chainType)], 
     Cluster_vertex_gini_distribution[,paste0("vertex_gini_",chainType)],
     col = cols[factor(Cluster_vertex_gini_distribution$outcome)],
     pch=20,ylab = c("Gini (Vextex)"),xlab = c("Gini (Cluster)"))
legend("topleft",legend = c("GTEX-Blood",
                                       "GTEX-Normal",
                                       "TCGA-normal-adj-pancreas",
                                       "TCGA-PDAC",
                                       "TCGA-pseudonormal_pancreas"), 
       col=cols,pch=20,cex=0.8)

par(fig=c(0,0.8,0.55,1), new=TRUE)
summary(glm(Cluster_vertex_gini_distribution[,paste0("cluster_gini_",chainType)]~Cluster_vertex_gini_distribution$outcome))
boxplot(Cluster_vertex_gini_distribution[,paste0("cluster_gini_",chainType)]~Cluster_vertex_gini_distribution$outcome,
        col=cols, horizontal=TRUE, axes=FALSE)



par(fig=c(0.65,1,0,0.8),new=TRUE)
summary(glm(Cluster_vertex_gini_distribution[,paste0("vertex_gini_",chainType)]~Cluster_vertex_gini_distribution$outcome))
boxplot(Cluster_vertex_gini_distribution[,paste0("vertex_gini_",chainType)]~Cluster_vertex_gini_distribution$outcome,
        col=cols, axes=FALSE)

dev.off()

tiff(paste0("Results/network_cluster_gini_",chainType,"_ALL_blood.tiff"),h=2000,w=2000,res=300)
ggboxplot(Cluster_vertex_gini_distribution, x = "outcome" , y =  "cluster_gini_IGK",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y= cluster_gini_IGK, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Blood",
                                                  "GTEX-Normal",
                                                  "TCGA-normal-adj-pancreas",
                                                  "TCGA-PDAC",
                                                  "TCGA-pseudonormal_pancreas")) +
  
  stat_compare_means(
    comparisons =list(c("GTEX-Blood","GTEX-Normal"),c("GTEX-Blood","TCGA-normal-adj-pancreas"),
                      c("GTEX-Blood","TCGA-PDAC"),c("GTEX-Blood","TCGA-pseudonormal_pancreas"),
                      c("TCGA-PDAC","TCGA-normal-adj-pancreas"),c("TCGA-PDAC","TCGA-pseudonormal_pancreas"),
                      c("TCGA-normal-adj-pancreas","TCGA-pseudonormal_pancreas")))

dev.off()

tiff(paste0("Results/network_vertex_gini_",chainType,"_ALL_blood.tiff"),h=2000,w=2000,res=300)
ggboxplot(Cluster_vertex_gini_distribution, x = "outcome" , y =  "vertex_gini_IGK",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y= vertex_gini_IGK, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Blood",
                                                  "GTEX-Normal",
                                                  "TCGA-normal-adj-pancreas",
                                                  "TCGA-PDAC",
                                                  "TCGA-pseudonormal_pancreas")) +
  
  stat_compare_means(
    comparisons =list(c("GTEX-Blood","GTEX-Normal"),c("GTEX-Blood","TCGA-normal-adj-pancreas"),
                      c("GTEX-Blood","TCGA-PDAC"),c("GTEX-Blood","TCGA-pseudonormal_pancreas"),
                      c("TCGA-PDAC","TCGA-normal-adj-pancreas"),c("TCGA-PDAC","TCGA-pseudonormal_pancreas"),
                      c("TCGA-normal-adj-pancreas","TCGA-pseudonormal_pancreas")))

dev.off()
