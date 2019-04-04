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

####Read repertoire data
PAAD_repertoire<-readRDS("Data/TCGA_Immune_Rep/PAAD_RepertoireResults.rds")
file_ids = rownames(PAAD_repertoire)
annotation = TCGAtranslateID(file_ids)
annotation$patient_barcode = substr(annotation$submitter_id,1,12)
write.csv(annotation,"Data/PAAD/annotation_PAAD.csv") ###To mark the non pancreas based on http://clincancerres.aacrjournals.org/content/clincanres/early/2018/05/08/1078-0432.CCR-18-0290.full.pdf
annotation<-read.csv("Data/PAAD/annotation_PAAD.csv")
id<-match(rownames(PAAD_repertoire),annotation$file_id)
rownames(PAAD_repertoire)<-annotation$submitter_id[id]
PAAD_repertoire$Tumor_type<-annotation$tumor_type[id]

####CLinical data
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-PAAD", 
                  data.category = "Clinical", 
                  file.type = "xml")
GDCdownload(query)
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
clinical.folow_up <- GDCprepare_clinic(query, clinical.info = "follow_up")
clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")
clinical.patient <- GDCprepare_clinic(query, clinical.info = "patient")
clinical.stage_event <- GDCprepare_clinic(query, clinical.info = "stage_event") 
clinical.new_tumor_event <- GDCprepare_clinic(query, clinical.info = "new_tumor_event")

###Biospecimen data
query <- GDCquery(project = "TCGA-PAAD", 
                  data.category = "Biospecimen", 
                  file.type = "xml")
GDCdownload(query)
biospecimen.sample <- GDCprepare_clinic(query, clinical.info = "sample")
biospecimen.bio_patient <- GDCprepare_clinic(query, clinical.info = "bio_patient")
biospecimen.slide <- GDCprepare_clinic(query, clinical.info = "slide")


#################################
###Cleaning the clinical data ###
################################
library(dplyr)
##Clinical.drug
clinical.drug<-clinical.drug[,-c(3,6:10,16:18)]
clinical.drug<-clinical.drug[!duplicated(clinical.drug), ]
match("Immunization",clinical.drug$clinical_trail_drug_classification) #212
clinical.drug$drug_name<-as.character(clinical.drug$drug_name)
clinical.drug$drug_name[212]<-"Immunization"
clinical.drug<-clinical.drug[,-9]
write.csv(clinical.drug,"Data/PAAD/Clinical/clinical_drug_PAAD.csv")

##Clinical.radiation
clinical.radiation<-clinical.radiation[!duplicated(clinical.radiation),]
clinical.radiation<-clinical.radiation[,-c(7,12:13,15)]
write.csv(clinical.radiation,"Data/PAAD/Clinical/clinical_radiation_PAAD.csv")

##Clinical.patient
clinical.patient<-clinical.patient[!duplicated(clinical.patient),]
clinical.patient$tobacco_smoking_history_master<-ifelse(clinical.patient$tobacco_smoking_history==1,"Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime)",
                                            ifelse(clinical.patient$tobacco_smoking_history==2,"Current smoker (includes daily smokers and non-daily smokers or occasional smokers)",
                                            ifelse(clinical.patient$tobacco_smoking_history==3,"Current reformed smoker for > 15 years (greater than 15 years)",
                                            ifelse(clinical.patient$tobacco_smoking_history==4,"Current reformed smoker for ≤15 years (less than or equal to 15 years)",
                                            ifelse(clinical.patient$tobacco_smoking_history==5,"Current reformed smoker, duration not specified",NA)))))
clinical.patient<-clinical.patient[,-c(2,72,75:80)]
write.csv(clinical.patient,"Data/PAAD/Clinical/clinical_patient_PAAD.csv")

##CLinical.stage has the same information than the one in clinical.patient

##clinical.new_tumor_event
clinical.new_tumor_event<-clinical.new_tumor_event[!duplicated(clinical.new_tumor_event),]
write.csv(clinical.new_tumor_event,"Data/PAAD/Clinical/clinical_new_tumor_event_PAAD.csv")


##CLinical.follow.up
clinical.folow_up<-read.csv("Data/PAAD/Clinical/Clinical_follow_up_fromSuppTableS1_PAAD.csv")
write.csv(clinical.folow_up,"Data/PAAD/Clinical/clinical_folow_up_PAAD.csv")


#biospecimen.sample
biospecimen.slide<-biospecimen.slide[!duplicated(biospecimen.slide),]
biospecimen.slide<-biospecimen.slide[,-c(2,4,10,13,15)]
write.csv(biospecimen.slide,"Data/PAAD/Clinical/biospecimen_slide_PAAD.csv")

###Read xCELL data
xCell.data<-read.csv("Data/xCell_TCGA_RSEM.csv")
rownames(xCell.data)<-xCell.data$X
xCell.data<-xCell.data[,-1]
xx<-str_replace(colnames(xCell.data),"\\.","-")
xx<-str_replace(xx,"\\.","-")
xx<-str_replace(xx,"\\.","-")
colnames(xCell.data)<-xx
write.csv(xCell.data,"Data/xCell_TCGA.csv")
id.xcell<-match(substring(rownames(PAAD_repertoire),1,15),colnames(xCell.data))
xCell.data.PAAD<-xCell.data[,na.omit(id.xcell)]

save(PAAD_repertoire,xCell.data.PAAD,clinical.drug,clinical.patient,clinical.radiation,clinical.new_tumor_event,clinical.folow_up,biospecimen.slide,annotation,
     file="Data/PAAD/PAAD_repertoire_xCell_clinical.Rdata")
