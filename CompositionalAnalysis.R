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
### Date: July, 2019
############################################################################################
library(zCompositions)
library(easyCODA)
library(ellipse)

##Colors: brewer.pal(4, "Accent")
##"#BEAED4" (IGH) "#7FC97F" (IGK) "#FDC086" (IGL)
setwd("~/TCGA-Immune/")

load("Data/PAAD/PAAD_FullData.Rdata")

data_merge_Ig_aa<-read.csv("Data/PAAD/ClonesInfered_Ig_aa.csv")
data_merge_TCR_aa<-read.csv("Data/PAAD/ClonesInfered_TCR_aa.csv")
data_merge_aa<-rbind(data_merge_Ig_aa,data_merge_TCR_aa)
data_merge<-merge(data_merge,data_merge_aa[,c("seqID","CloneId")],by=c("seqID"))
data_merge$V_J_lenghCDR3aa_CloneId<-paste(data_merge$V_J_lenghCDR3aa,data_merge$CloneId.y,sep="_")

##Filter for those that are pancreas or normal-pseudonormal
data_merge_qc<-data_merge[which(is.na(match(data_merge$sample,rownames(PAAD.repertoire.diversity)))==F),]

##Filter for those that are only pancreas
##Filter for those that are pancreas or normal-pseudonormal
PAAD.repertoire.diversity.tumor<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type_2categ=="Tumor_pancreas"),]
data_merge_pancreas<-data_merge[which(is.na(match(data_merge$sample,rownames(PAAD.repertoire.diversity.tumor)))==F),]

############
## 1. Build the matrix with the clones by samples
###########
## Pancreas ans pseudonormal
data_qc_chain_IGK<-data_merge_qc[which(data_merge_qc$chainType=="IGK"),]
clone_type_IGK<-t(as.data.frame(unclass(table(data_qc_chain_IGK$V_J_lenghCDR3aa_CloneId,factor(data_qc_chain_IGK$sample))))) 
data_qc_chain_IGL<-data_merge_qc[which(data_merge_qc$chainType=="IGL"),]
clone_type_IGL<-t(as.data.frame(unclass(table(data_qc_chain_IGL$V_J_lenghCDR3aa_CloneId,factor(data_qc_chain_IGL$sample))))) 
data_qc_chain_IGH<-data_merge_qc[which(data_merge_qc$chainType=="IGH"),]
clone_type_IGH<-t(as.data.frame(unclass(table(data_qc_chain_IGH$V_J_lenghCDR3aa_CloneId,factor(data_qc_chain_IGH$sample))))) 
data_qc_chain_IG<-data_merge_qc[which(data_merge_qc$chainType=="IGH" |
                                        data_merge_qc$chainType=="IGK" |
                                        data_merge_qc$chainType=="IGL"),]
clone_type_IG<-t(as.data.frame(unclass(table(data_qc_chain_IG$V_J_lenghCDR3aa_CloneId,factor(data_qc_chain_IG$sample))))) 

data_qc_chain_TCR<-data_merge_qc[which(data_merge_qc$chainType=="TRA" | 
                                         data_merge_qc$chainType=="TRB" |
                                         data_merge_qc$chainType=="TRD" | 
                                         data_merge_qc$chainType=="TRG"),]
clone_type_TCR<-t(as.data.frame(unclass(table(data_qc_chain_TCR$V_J_lenghCDR3aa_CloneId,factor(data_qc_chain_TCR$sample))))) 


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

                              