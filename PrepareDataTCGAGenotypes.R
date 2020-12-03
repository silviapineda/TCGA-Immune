
rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. Genotypes~PAAD 
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Preparing Immune and Clinical for PAAD 
###         
###
### Author: Silvia Pineda
### Date: April, 2019
############################################################################################

setwd("~/TCGA-Immune/")
load("Data/PAAD/PAAD_clonotypes_TCGA.Rdata")
load("Data/PAAD/PAAD_FullData.Rdata")
SNPsPAAD<-read.table("Data/PAAD/Genotypes/callSNPsPAAD.txt",header = T)
colnames(SNPsPAAD)<-SNPsPAAD[1,]
rownames(SNPsPAAD)<-SNPsPAAD[,1]
SNPsPAAD<-SNPsPAAD[-1,-1]

id<-match(substr(PAAD.repertoire.tumor.filter$TCGA_sample,1,12),substr(colnames(SNPsPAAD),1,12))
Genotypes.PAAD.filter<-SNPsPAAD[,na.omit(id)]

id<-match(substr(colnames(SNPsPAAD),1,12),substr(PAAD.repertoire.tumor.filter$TCGA_sample,1,12))
PAAD.repertoire.tumor.filter.snp<-PAAD.repertoire.tumor.filter[na.omit(id),]

id<-match(rownames(PAAD.repertoire.tumor.filter.snp),rownames(clone_type_filter_IGK))
clone_type_filter_IGK_snp<-clone_type_filter_IGK[na.omit(id),]
rownames(clone_type_filter_IGK_snp)<-PAAD.repertoire.tumor.filter.snp$TCGA_sample

id<-match(rownames(PAAD.repertoire.tumor.filter.snp),rownames(clone_type_filter_IGL))
clone_type_filter_IGL_snp<-clone_type_filter_IGL[na.omit(id),]
rownames(clone_type_filter_IGL_snp)<-PAAD.repertoire.tumor.filter.snp$TCGA_sample

id<-match(rownames(PAAD.repertoire.tumor.filter.snp),rownames(clone_type_filter_IGH))
clone_type_filter_IGH_snp<-clone_type_filter_IGH[na.omit(id),]
rownames(clone_type_filter_IGH_snp)<-PAAD.repertoire.tumor.filter.snp$TCGA_sample


save(PAAD.repertoire.tumor.filter.snp,Genotypes.PAAD.filter,clone_type_filter_IGH_snp,clone_type_filter_IGK_snp,clone_type_filter_IGL_snp,
     file="Data/PAAD/PAAD_Repertoire_Genotypes.Rdata")

