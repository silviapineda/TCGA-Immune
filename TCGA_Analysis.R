rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. Prepare data 
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

load("Data/PAAD/PAAD_FullData.Rdata")
PAAD_repertoire_diversity$Tumor_type_2categ<-ifelse(PAAD_repertoire_diversity$Tumor_type=="Tumor_pancreas","Tumor_pancres",
                                             ifelse(PAAD_repertoire_diversity$Tumor_type=="Solid_tissue_normal","Adjacent_normal_pancreas",
                                             ifelse(PAAD_repertoire_diversity$Tumor_type=="Adjacent_normal_pancreas","Adjacent_normal_pancreas","Other")))

####Descriptive analysis to see if there are differences by tumor and adjacent_normal
p.value<-NULL
for(i in 1:(ncol(PAAD_repertoire_diversity)-3)){
  p.value[i]<-coef(summary(glm(PAAD_repertoire_diversity[,i]~PAAD_repertoire_diversity$Tumor_type_2categ)))[3,4]
}
boxplot(PAAD_repertoire_diversity$entropy_recon_IGK~PAAD_repertoire_diversity$Tumor_type_2categ)

##Only tumor pancreas
PAAD_repertoire_tumor<-PAAD_repertoire_diversity[which(PAAD_repertoire_diversity$Tumor_type=="Tumor_pancreas"),]
PAAD_repertoire_tumor$bcr_patient_barcode<-substr(PAAD_repertoire_tumor$TCGA_sample,1,12)

##Merge with xCell
xCell.data.PAAD<-data.frame(t(xCell.data.PAAD))
xCell.data.PAAD$bcr_patient_barcode<-substr(rownames(xCell.data.PAAD),1,12)
xcell.repertoire.data<-merge(PAAD_repertoire_tumor,xCell.data.PAAD,by="bcr_patient_barcode")

xcell.data<-xcell.repertoire.data[,54:length(colnames(xcell.repertoire.data))]
p.value=NULL
for (i in 1:length(colnames(xcell.data))){
  p.value[i]<-coef(summary(glm(xcell.repertoire.data$clones_recon_TRA~xcell.data[,i])))[2,4]
}
xcell.data[,which(p.adjust(p.value)<0.05)]

plot(xcell.repertoire.data$IG_expression~xcell.data$B.cells)





###Merge with Clinical data
repertoire.patient.data<-merge(PAAD_repertoire_tumor,clinical.patient,by="bcr_patient_barcode")
plot(repertoire.patient.data$IG_expression~repertoire.patient.data$lymph_node_examined_count)
summary(glm(repertoire.patient.data$IG_expression~repertoire.patient.data$primary_lymph_node_presentation_assessment))


