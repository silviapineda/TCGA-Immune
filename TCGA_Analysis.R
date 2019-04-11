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

load("Data/PAAD/PAAD_repertoire_xCell_clinical.Rdata")

boxplot(PAAD_repertoire$Alpha_Beta_ratio_expression~PAAD_repertoire$Tumor_type,las=2)

##Only tumor pancreas
PAAD_repertoire_tumor<-PAAD_repertoire[which(PAAD_repertoire$Tumor_type=="Tumor_pancreas"),]
PAAD_repertoire_tumor$bcr_patient_barcode<-substr(rownames(PAAD_repertoire_tumor),1,12)

##Merge with xCell
xCell.data.PAAD<-data.frame(t(xCell.data.PAAD))
xCell.data.PAAD$bcr_patient_barcode<-substr(rownames(xCell.data.PAAD),1,12)
xcell.repertoire.data<-merge(PAAD_repertoire_tumor,xCell.data.PAAD,by="bcr_patient_barcode")

xcell.data<-xcell.repertoire.data[,24:length(colnames(xcell.repertoire.data))]
p.value=NULL
for (i in 1:length(colnames(xcell.data))){
  p.value[i]<-coef(summary(glm(xcell.repertoire.data$IG_expression~xcell.data[,i])))[2,4]
}
xcell.data[,which(p.adjust(p.value)<0.05)]

plot(xcell.repertoire.data$IG_expression~xcell.data$B.cells)





###Merge with Clinical data
repertoire.patient.data<-merge(PAAD_repertoire_tumor,clinical.patient,by="bcr_patient_barcode")
plot(repertoire.patient.data$IG_expression~repertoire.patient.data$lymph_node_examined_count)
summary(glm(repertoire.patient.data$IG_expression~repertoire.patient.data$primary_lymph_node_presentation_assessment))


