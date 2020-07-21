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
library(zCompositions)
library(easyCODA)
library(ellipse)
library(factoextra)
library(RColorBrewer)
library(sva)

##Colors: brewer.pal(4, "Accent")
##"#BEAED4" (IGH) "#7FC97F" (IGK) "#FDC086" (IGL)
setwd("~/TCGA-Immune/")

load("Data/PAAD_GTEx/PAAD_GTEx_FullData.Rdata")

##Filter for those that are pancreas or normal
PAAD.GTEx.repertoire.diversity.tumor.normal<-PAAD.GTEx.repertoire.diversity[which(PAAD.GTEx.repertoire.diversity$outcome=="normal-pancreas (GTEx)" |
                                                                                  PAAD.GTEx.repertoire.diversity$outcome=="tumor-pancreas (TCGA)"),]
data_merge_qc<-data_merge[which(is.na(match(data_merge$sample,rownames(PAAD.GTEx.repertoire.diversity.tumor.normal)))==F),]

############
## 1. Build the matrix with the clones by samples
###########

data_qc_chain_IGK<-data_merge_qc[which(data_merge_qc$chainType=="IGK"),]
clone_type_IGK<-t(as.data.frame(unclass(table(data_qc_chain_IGK$V_J_lenghCDR3_CloneId,factor(data_qc_chain_IGK$sample))))) 
data_qc_chain_IGL<-data_merge_qc[which(data_merge_qc$chainType=="IGL"),]
clone_type_IGL<-t(as.data.frame(unclass(table(data_qc_chain_IGL$V_J_lenghCDR3_CloneId,factor(data_qc_chain_IGL$sample))))) 
data_qc_chain_IGH<-data_merge_qc[which(data_merge_qc$chainType=="IGH"),]
clone_type_IGH<-t(as.data.frame(unclass(table(data_qc_chain_IGH$V_J_lenghCDR3_CloneId,factor(data_qc_chain_IGH$sample))))) 
data_qc_chain_IG<-data_merge_qc[which(data_merge_qc$chainType=="IGH" |
                                        data_merge_qc$chainType=="IGK" |
                                        data_merge_qc$chainType=="IGL"),]
clone_type_IG<-t(as.data.frame(unclass(table(data_qc_chain_IG$V_J_lenghCDR3_CloneId,factor(data_qc_chain_IG$sample))))) 

data_qc_chain_TCR<-data_merge_qc[which(data_merge_qc$chainType=="TRA" | 
                                         data_merge_qc$chainType=="TRB" |
                                         data_merge_qc$chainType=="TRD" | 
                                         data_merge_qc$chainType=="TRG"),]
clone_type_TCR<-t(as.data.frame(unclass(table(data_qc_chain_TCR$V_J_lenghCDR3_CloneId,factor(data_qc_chain_TCR$sample))))) 

### Calculate relative abundances and logtransform
clone_type_IG_relative <- clone_type_IG / colSums(clone_type_IG)
min<-min(clone_type_IG_relative[clone_type_IG_relative!=0])/2
clone_type_IG_relative <- clone_type_IG_relative+min
clone_type_IG_relative_log<-log2(clone_type_IG_relative)

#####################################
### Principal Components Analysis ###
#####################################
pca <- prcomp(t(clone_type_IG_relative_log), scale = TRUE)

SPP <- factor(PAAD.GTEx.repertoire.diversity.tumor.normal$outcome)
levels.SPP <- factor(c("normal-pancreas (GTEx)", "tumor-pancreas (TCGA)"))

brewer.pal(4,name = "Accent")
cols=c( "#7FC97F","#BEAED4")


pc <- c(1,2)
tiff(paste0("Results/PCA_clones_relative_abundance_IG.tiff"),width = 2000, height = 2000, res = 300)
plot(pca$x[,pc[1]], pca$x[,pc[2]], col=cols[SPP],pch=20,xlab="PCA1",ylab="PCA2")
legend("bottomright", legend=levels(levels.SPP), col=cols,pch=20,cex=0.8)
dev.off()


### Correct by batch effect using Combat
combat_edata1 = ComBat(dat=clone_type_IG_relative_log, batch=PAAD.GTEx.repertoire.diversity.tumor.normal$outcome, mod=NULL, par.prior=TRUE, prior.plots=FALSE)



###Filter by clones share at least in 2 samples
clone_type_filter_IGK<-clone_type_IGK[,which(colSums(clone_type_IGK!=0)>=2)] #
clone_type_filter_IGL<-clone_type_IGL[,which(colSums(clone_type_IGL!=0)>=2)] #
clone_type_filter_IGH<-clone_type_IGH[,which(colSums(clone_type_IGH!=0)>=2)] #
clone_type_filter_IG<-clone_type_IG[,which(colSums(clone_type_IG!=0)>=2)] #
clone_type_filter_TCR<-clone_type_TCR[,which(colSums(clone_type_TCR!=0)>=2)] #



###################
### Zcomposition ###
###################
###Applied Geometric Bayesian Multiplicative to substitute the zeros
clone_type_IG.GBM <- cmultRepl(t(clone_type_filter_IG), output="p-counts")
clone_type_IG.GBM.pro<-as.matrix(t(clone_type_IG.GBM))/colSums(clone_type_IG.GBM)

###Apllied the LRA (logratios Analysis)
clone_type_IG.GBM.lra<-LRA(clone_type_IG.GBM.pro)
#Total inertia
sum(clone_type_IG.GBM.lra$sv^2) ##Total Variance = 1.65
#Percentages of inertia
100*clone_type_IG.GBM.lra$sv^2/sum(clone_type_IG.GBM.lra$sv^2) ##amount of variance explain by components #10.5 (first) 1.9 (second)

#Plot results from LRA
COLORS =c(brewer.pal(3,"Accent")[1],brewer.pal(3,"Accent")[2])
COLORS2=c(brewer.pal(3,"Dark2")[1],brewer.pal(3,"Dark2")[3])
tiff("Results/LRA_plot_weighted.tiff",width = 2000, height = 2000, res = 300)
plot(clone_type_IG.GBM.lra$rowpcoord,col=COLORS[PAAD.repertoire.diversity$Tumor_type_2categ],
     xlab="Dimension 1 (10.5%)", ylab="Dimension 2 (1.9%)",pch=20,cex.axis=1.2,cex.lab=1.2)
CIplot_biv(clone_type_IG.GBM.lra$rowpcoor[,1], clone_type_IG.GBM.lra$rowpcoor[,2], group=PAAD.repertoire.diversity$Tumor_type_2categ, 
           groupcols=COLORS2,add=T,shade=T)
dev.off()

############################################
####Perform the analysis in only pancreas###
###########################################
PAAD.repertoire.tumor<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type_2categ=="Tumor_pancreas"),]
clone_type_filter_IG_tumor<-clone_type_filter_IG[match(rownames(PAAD.repertoire.tumor),rownames(clone_type_filter_IG)),]
clone_type_filter_IG_tumor<-clone_type_filter_IG_tumor[,which(colSums(clone_type_filter_IG_tumor!=0)>=2)] #


###Applied Geometric Bayesian Multiplicative to substitute the zeros
clone_type_IG_tumor.GBM <- cmultRepl(t(clone_type_filter_IG_tumor), output="p-counts")
clone_type_IG_tumor.GBM.pro<-as.matrix(t(clone_type_IG_tumor.GBM))/colSums(clone_type_IG_tumor.GBM)

###Apllied the LRA (logratios Analysis)
clone_type_IG_tumor.GBM.lra<-LRA(clone_type_IG_tumor.GBM.pro)
#Total inertia
sum(clone_type_IG_tumor.GBM.lra$sv^2) ##Total Variance = 1.69
#Percentages of inertia
100*clone_type_IG_tumor.GBM.lra$sv^2/sum(clone_type_IG_tumor.GBM.lra$sv^2) ##amount of variance explain by components #10.5 (first) 1.9 (second)

##Extract Clinical outcomes only for tumor
clinical.patient.tumor<-clinical.patient[match(substr(PAAD.repertoire.tumor$TCGA_sample,1,12),clinical.patient$bcr_patient_barcode),]
#Plot results from LRA
COLORS =c(brewer.pal(3,"Accent")[1],brewer.pal(3,"Accent")[2])
COLORS2=c(brewer.pal(3,"Dark2")[1],brewer.pal(3,"Dark2")[3])
tiff("Results/LRA_plot_gender.tiff",width = 2000, height = 2000, res = 300)
plot(clone_type_IG_tumor.GBM.lra$rowpcoord,col=COLORS[clinical.patient.tumor$gender],
     xlab="Dimension 1 (10.5%)", ylab="Dimension 2 (1.9%)",pch=20,cex.axis=1.2,cex.lab=1.2)
CIplot_biv(clone_type_IG_tumor.GBM.lra$rowpcoor[,1], clone_type_IG.GBM.lra$rowpcoor[,2], group=clinical.patient.tumor$gender, 
           groupcols=COLORS2,add=T,shade=T)
dev.off()


###############################
## Correspondance analysis  ##

# sample profiles (relative abundances), also transposing matrix
# no zero substitution
clone_type_filter_IG.pro <- clone_type_filter_IG / colSums(clone_type_filter_IG)

# correspondence analysis using ca() function in ca package
# also CA() function in easyCODA can be used 
# (at the moment, both are available in easyCODA)
clone_type_IG.ca <- CA(clone_type_filter_IG.pro)
#Total inertia
sum(clone_type_IG.ca$sv^2) ##Total Variance = 12.o
#Percentages of inertia
100*clone_type_IG.ca$sv^2/sum(clone_type_IG.ca$sv^2) ##amount of variance explain by components #6.5 (first) 5.1 (second)

#Plot results from LRA
COLORS =c(brewer.pal(3,"Accent")[1],brewer.pal(3,"Accent")[2])
COLORS2=c(brewer.pal(3,"Dark2")[1],brewer.pal(3,"Dark2")[3])
tiff("Results/CA_plot.tiff",width = 2000, height = 2000, res = 300)
plot(clone_type_IG.ca$rowpcoord,col=COLORS[PAAD.repertoire.diversity$Tumor_type_2categ],
     xlab="Dimension 1 (6.5%)", ylab="Dimension 2 (5.1%)",pch=20,cex.axis=1.2,cex.lab=1.2)
CIplot_biv(clone_type_IG.ca$rowpcoor[,1], clone_type_IG.GBM.lra$rowpcoor[,2], group=PAAD.repertoire.diversity$Tumor_type_2categ, 
           groupcols=COLORS2,add=T,shade=T)
dev.off()


######
##Transformation to logratios 
clone_type_filter_IG.lr <- LR(clone_type_IG.GBM.pro)


#####Clustering
clone_type_IG.clus <- hclust(t(clone_type_filter_IG.pro), method="ward.D2")
plot(marina.clus, cex=0.5)

marina.sub.LR <- LR(genera0.GBM.pro[,marina.ind], weight=FALSE)$LR
dim(marina.sub.LR)

                              
