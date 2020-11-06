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
library(pgirmess)
library(stringr)
library(SciViews)

setwd("~/TCGA-Immune/")

########## Read the alignment files
files <- list.files("Data/PAAD/MIXCR_PAAD/alignments/")
data<-c()
for(i in files) {
  cat(i, "\n")
  t <- read.delim(paste0("Data/PAAD/MIXCR_PAAD/alignments/",i))
  t$sample<-substr(i, 1, nchar(i)-15)
  data <- rbind(data, t)
}
save(data,file="Data/PAAD/PAAD_MIXCR_results.Rdata")

#####
load("Data/PAAD/PAAD_MIXCR_results.Rdata")
###Chain
data$chainType<-ifelse(data$bestVGene!="", substr(data$bestVGene,1,3),
                       ifelse(data$bestVGene=="",substr(data$bestJGene,1,3),NA))

###Add a column with the reads per chain for Ig and TRA
##Reads per chain
read_count <- table(data$sample)
read_count_chain <- table(data$sample, data$chainType)
reads <- data.frame(cbind(read_count,read_count_chain))

####### The data needs to be normalized by the unmapped reads 
totalReads<-read.table("Data/PAAD/MIXCR_PAAD/total_reads.txt",sep=";") ##We need to extract this number from the MIXCR report with the python script
id<-match(rownames(reads),totalReads$V1)
reads$totalReads<-totalReads[id,2]

##Total reads
reads$Ig_Reads<-reads$IGH+reads$IGK+reads$IGL
reads$T_Reads<- reads$TRA+reads$TRB+reads$TRD+reads$TRG

####Normalize the nuber of reads
reads$IG_expression<-(reads$IGH+reads$IGK+reads$IGL)/reads$totalReads
reads$IGH_expression<-reads$IGH/reads$totalReads
reads$IGK_expression<-reads$IGK/reads$totalReads
reads$IGL_expression<-reads$IGL/reads$totalReads

reads$T_expression<-(reads$TRA+reads$TRB+reads$TRD+reads$TRG)/reads$totalReads
reads$TRA_expression<-reads$TRA/reads$totalReads
reads$TRB_expression<-reads$TRB/reads$totalReads
reads$TRD_expression<-reads$TRD/reads$totalReads
reads$TRG_expression<-reads$TRG/reads$totalReads
###Ratio
reads$Alpha_Beta_ratio_expression<-(reads$TRA_expression+reads$TRB_expression)/reads$T_expression
reads$KappaLambda_ratio_expression <- (reads$IGK_expression / reads$IGL_expression)

####Restriction to CDR3 to extract the clones
data_full_cdr3<-data[which(data$nSeqCDR3!=""),] #2,087,012
data_full_cdr3$CDR3_length<-nchar(as.character(data_full_cdr3$nSeqCDR3)) 
data_full_cdr3<-data_full_cdr3[which(data_full_cdr3$CDR3_length!=3),] #2,087,010

##Obtain the ID for the clone call
data_full_cdr3$seqID<-seq(1,nrow(data_full_cdr3))
data_full_cdr3$SEQUENCE_ID<-paste(data_full_cdr3$sample,data_full_cdr3$seqID,sep="_")

##Variable to build clones
data_full_cdr3$V_J_lenghCDR3 = paste(data_full_cdr3$bestVGene,data_full_cdr3$bestJGene,data_full_cdr3$CDR3_length,sep="_")

###Chain
data_full_cdr3_Ig<-data_full_cdr3[which(data_full_cdr3$chainType=="IGH" |
                                          data_full_cdr3$chainType=="IGK" |
                                          data_full_cdr3$chainType=="IGL"),]
data_full_cdr3_TCR<-data_full_cdr3[which(data_full_cdr3$chainType=="TRA" | 
                                           data_full_cdr3$chainType=="TRB" |
                                           data_full_cdr3$chainType=="TRD" | 
                                           data_full_cdr3$chainType=="TRG"),]

###save the data to call the clones by all samples using the nucleotides.py
data_clonesInference_Ig<-data_full_cdr3_Ig[,c("SEQUENCE_ID","sample","nSeqCDR3","CDR3_length","bestVGene","bestJGene","V_J_lenghCDR3")]
data_clonesInference_TCR<-data_full_cdr3_TCR[,c("SEQUENCE_ID","sample","nSeqCDR3","CDR3_length","bestVGene","bestJGene","V_J_lenghCDR3")]
write.table(data_clonesInference_Ig,file="Data/PAAD/data_for_cloneInfered_Ig_PAAD.txt",row.names = F,sep="\t")
write.table(data_clonesInference_TCR,file="Data/PAAD/data_for_cloneInfered_TCR_PAAD.txt",row.names = F,sep="\t")

### After passing the nucleotides.py
##Read the clones and merge with the data
nucleotides_Ig<-read.csv("Data/PAAD/ClonesInfered_PAAD_Ig.csv")
nucleotides_TCR<-read.csv("Data/PAAD/ClonesInfered_PAAD_TCR.csv")
nucleotides<-rbind(nucleotides_Ig,nucleotides_TCR)
data_merge<-merge(data_full_cdr3,nucleotides[,c("SEQUENCE_ID","CloneId")],by=c("SEQUENCE_ID"))

####Reads per chain only for clonal information
##Reads per chain
read_count <- table(data_merge$sample)
read_count_chain <- table(data_merge$sample, data_merge$chainType)
reads_filter <- data.frame(cbind(read_count,read_count_chain))
####### The data needs to be normalized by the unmapped reads 
totalReads<-read.table("Data/PAAD/MIXCR_PAAD/total_reads.txt",sep=";") ##We need to extract this number from the MIXCR report with the python script
id<-match(rownames(reads_filter),totalReads$V1)
reads_filter$totalReads<-totalReads[id,2]

##Total reads
reads_filter$Ig_Reads<-reads_filter$IGH+reads_filter$IGK+reads_filter$IGL
reads_filter$T_reads_filter<- reads_filter$TRA+reads_filter$TRB+reads_filter$TRD+reads_filter$TRG

####Normalize the nuber of reads_filter
reads_filter$IG_expression<-(reads_filter$IGH+reads_filter$IGK+reads_filter$IGL)/reads_filter$totalReads
reads_filter$IGH_expression<-reads_filter$IGH/reads_filter$totalReads
reads_filter$IGK_expression<-reads_filter$IGK/reads_filter$totalReads
reads_filter$IGL_expression<-reads_filter$IGL/reads_filter$totalReads

reads_filter$T_expression<-(reads_filter$TRA+reads_filter$TRB+reads_filter$TRD+reads_filter$TRG)/reads_filter$totalReads
reads_filter$TRA_expression<-reads_filter$TRA/reads_filter$totalReads
reads_filter$TRB_expression<-reads_filter$TRB/reads_filter$totalReads
reads_filter$TRD_expression<-reads_filter$TRD/reads_filter$totalReads
reads_filter$TRG_expression<-reads_filter$TRG/reads_filter$totalReads
###Ratio
reads_filter$Alpha_Beta_ratio_expression<-(reads_filter$TRA_expression+reads_filter$TRB_expression)/reads_filter$T_expression
reads_filter$KappaLambda_ratio_expression <- (reads_filter$IGK_expression / reads_filter$IGL_expression)

colnames(reads_filter)<-c("read_count_filter","IGH_filter","IGK_filter", "IGL_filter","TRA_filter","TRB_filter","TRD_filter","TRG_filter", 
                          "totalReads","Ig_Reads_filter","T_Reads_filter","IG_expression_filter", "IGH_expression_filter","IGK_expression_filter",
                          "IGL_expression_filter", "T_expression_filter" ,               
                          "TRA_expression_filter", "TRB_expression_filter","TRD_expression_filter","TRG_expression_filter",
                          "Alpha_Beta_ratio_expression_filter",  "KappaLambda_ratio_expression_filter")



##Clones per chain
data_merge$V_J_lenghCDR3_CloneId = paste(data_merge$V_J_lenghCDR3,data_merge$CloneId,sep="_")

clones_count<- unique(data_merge[,c("sample","V_J_lenghCDR3_CloneId","chainType")])
clones<-data.frame(cbind(table(clones_count$sample,clones_count$chainType)))
colnames(clones)<-c("clones_IGH","clones_IGK","clones_IGL","clones_TRA","clones_TRB","clones_TRD","clones_TRG")

##Diversity measures
sample<-rownames(clones)
entropy_IGH<-NULL
entropy_IGK<-NULL
entropy_IGL<-NULL
entropy_TRA<-NULL
entropy_TRB<-NULL
entropy_TRD<-NULL
entropy_TRG<-NULL
for (i in 1:length(sample)){
  print(i)
  data_sample_unique<-data_merge[which(data_merge$sample==sample[i]),]
  clones_sample<-data_sample_unique[,"V_J_lenghCDR3_CloneId"]
  clones_sample_IGH<-data_sample_unique[which(data_sample_unique$chainType=="IGH"),"V_J_lenghCDR3_CloneId"]
  clones_sample_IGK<-data_sample_unique[which(data_sample_unique$chainType=="IGK"),"V_J_lenghCDR3_CloneId"]
  clones_sample_IGL<-data_sample_unique[which(data_sample_unique$chainType=="IGL"),"V_J_lenghCDR3_CloneId"]
  clones_sample_TRA<-data_sample_unique[which(data_sample_unique$chainType=="TRA"),"V_J_lenghCDR3_CloneId"]
  clones_sample_TRB<-data_sample_unique[which(data_sample_unique$chainType=="TRB"),"V_J_lenghCDR3_CloneId"]
  clones_sample_TRD<-data_sample_unique[which(data_sample_unique$chainType=="TRD"),"V_J_lenghCDR3_CloneId"]
  clones_sample_TRG<-data_sample_unique[which(data_sample_unique$chainType=="TRG"),"V_J_lenghCDR3_CloneId"]
  
  #To write file to run with Recon
  #write.delim(data.frame(table(table(clones_sample_IGH))),file=paste("clones_sample_IGH_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  #write.delim(data.frame(table(table(clones_sample_IGK))),file=paste("clones_sample_IGK_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  #write.delim(data.frame(table(table(clones_sample_IGL))),file=paste("clones_sample_IGL_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.delim(data.frame(table(table(clones_sample_TRA))),file=paste("clones_sample_TRA_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.delim(data.frame(table(table(clones_sample_TRB))),file=paste("clones_sample_TRB_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.delim(data.frame(table(table(clones_sample_TRD))),file=paste("clones_sample_TRD_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.delim(data.frame(table(table(clones_sample_TRG))),file=paste("clones_sample_TRG_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  
  fi_IGH<-as.numeric(table(clones_sample_IGH))/length(clones_sample_IGH)
  fi_IGK<-as.numeric(table(clones_sample_IGK))/length(clones_sample_IGK)
  fi_IGL<-as.numeric(table(clones_sample_IGL))/length(clones_sample_IGL)
  fi_TRA<-as.numeric(table(clones_sample_TRA))/length(clones_sample_TRA)
  fi_TRB<-as.numeric(table(clones_sample_TRB))/length(clones_sample_TRB)
  fi_TRD<-as.numeric(table(clones_sample_TRD))/length(clones_sample_TRD)
  fi_TRG<-as.numeric(table(clones_sample_TRG))/length(clones_sample_TRG)
  
  hi_IGH<-fi_IGH*log2(fi_IGH)
  hi_IGK<-fi_IGK*log2(fi_IGK)
  hi_IGL<-fi_IGL*log2(fi_IGL)
  hi_TRA<-fi_TRA*log2(fi_TRA)
  hi_TRB<-fi_TRB*log2(fi_TRB)
  hi_TRD<-fi_TRD*log2(fi_TRD)
  hi_TRG<-fi_TRG*log2(fi_TRG)
  
  entropy_IGH[i]=-sum(hi_IGH)
  entropy_IGK[i]=-sum(hi_IGK)
  entropy_IGL[i]=-sum(hi_IGL)
  entropy_TRA[i]=-sum(hi_TRA)
  entropy_TRB[i]=-sum(hi_TRB)
  entropy_TRD[i]=-sum(hi_TRD)
  entropy_TRG[i]=-sum(hi_TRG)
  
}

diversity<-cbind(clones,entropy_IGH,entropy_IGK,entropy_IGL,entropy_TRA,entropy_TRB,entropy_TRD,entropy_TRG)

####After runing recon
recon<-read.table("Data/PAAD/RECON/test_D_number_table.txt",header=T)
chain<-substr(recon$sample_name,61,63)
sample<-substr(recon$sample_name,65,100)
##0.0D is species richness (Number of clones)
##Entropy is ln(1.0.D)
diversity<-as.data.frame(diversity)
chain_list<-unique(chain)
for(i in chain_list){
  recon_chain<-recon[which(chain==i),]
  sample_chain<-substr(recon_chain$sample_name,65,100)
  id<-match(rownames(diversity),sample_chain)
  clone_chain<-ifelse(is.na(id)==F,recon_chain[id,"est_0.0D"],0)
  assign(paste0("clones_recon_",i),clone_chain)
  entroy_chain<-ifelse(is.na(id)==F,ln(recon_chain[id,"est_1.0D"]),0)
  assign(paste0("entropy_recon_",i),entroy_chain)
}
diversity$clones_recon_IGH<-clones_recon_IGH
diversity$clones_recon_IGK<-clones_recon_IGK
diversity$clones_recon_IGL<-clones_recon_IGL
diversity$clones_recon_TRA<-clones_recon_TRA
diversity$clones_recon_TRB<-clones_recon_TRB
diversity$clones_recon_TRD<-clones_recon_TRD
diversity$clones_recon_TRG<-clones_recon_TRG

diversity$entropy_recon_IGH<-entropy_recon_IGH
diversity$entropy_recon_IGK<-entropy_recon_IGK
diversity$entropy_recon_IGL<-entropy_recon_IGL
diversity$entropy_recon_TRA<-entropy_recon_TRA
diversity$entropy_recon_TRB<-entropy_recon_TRB
diversity$entropy_recon_TRD<-entropy_recon_TRD
diversity$entropy_recon_TRG<-entropy_recon_TRG
#Simpson Index (1/2D) and BPI (1/∞D)

##length of CDR3
cdr3_length<-aggregate(data_merge$CDR3_length,by=list(data_merge$sample,data_merge$chainType), FUN=mean)

cdr3_length_IGH<-cdr3_length[which(cdr3_length$Group.2=="IGH"),]
cdr3_length_IGK<-cdr3_length[which(cdr3_length$Group.2=="IGK"),]
cdr3_length_IGL<-cdr3_length[which(cdr3_length$Group.2=="IGL"),]
cdr3_length_TRA<-cdr3_length[which(cdr3_length$Group.2=="TRA"),]
cdr3_length_TRB<-cdr3_length[which(cdr3_length$Group.2=="TRB"),]
cdr3_length_TRD<-cdr3_length[which(cdr3_length$Group.2=="TRD"),]
cdr3_length_TRG<-cdr3_length[which(cdr3_length$Group.2=="TRG"),]

id.IGH<-match(rownames(diversity),cdr3_length_IGH$Group.1)
diversity$cdr3_length_IGH<-cdr3_length$x[id.IGH]
diversity$cdr3_length_IGH<-replace(diversity$cdr3_length_IGH,is.na(diversity$cdr3_length_IGH)==T,0)
id.IGK<-match(rownames(diversity),cdr3_length_IGK$Group.1)
diversity$cdr3_length_IGK<-cdr3_length$x[id.IGK]
diversity$cdr3_length_IGK<-replace(diversity$cdr3_length_IGK,is.na(diversity$cdr3_length_IGK)==T,0)
id.IGL<-match(rownames(diversity),cdr3_length_IGL$Group.1)
diversity$cdr3_length_IGL<-cdr3_length$x[id.IGL]
diversity$cdr3_length_IGL<-replace(diversity$cdr3_length_IGL,is.na(diversity$cdr3_length_IGL)==T,0)
id.TRA<-match(rownames(diversity),cdr3_length_TRA$Group.1)
diversity$cdr3_length_TRA<-cdr3_length$x[id.TRA]
diversity$cdr3_length_TRA<-replace(diversity$cdr3_length_TRA,is.na(diversity$cdr3_length_TRA)==T,0)
id.TRB<-match(rownames(diversity),cdr3_length_TRB$Group.1)
diversity$cdr3_length_TRB<-cdr3_length$x[id.TRB]
diversity$cdr3_length_TRB<-replace(diversity$cdr3_length_TRB,is.na(diversity$cdr3_length_TRB)==T,0)
id.TRD<-match(rownames(diversity),cdr3_length_TRD$Group.1)
diversity$cdr3_length_TRD<-cdr3_length$x[id.TRD]
diversity$cdr3_length_TRD<-replace(diversity$cdr3_length_TRD,is.na(diversity$cdr3_length_TRD)==T,0)
id.TRG<-match(rownames(diversity),cdr3_length_TRG$Group.1)
diversity$cdr3_length_TRG<-cdr3_length$x[id.TRG]
diversity$cdr3_length_TRG<-replace(diversity$cdr3_length_TRG,is.na(diversity$cdr3_length_TRG)==T,0)


PAAD_repertoire_diversity<-cbind(reads,reads_filter,diversity)

save(data_merge,PAAD_repertoire_diversity,file="Data/PAAD/PAAD_RepertoireResults_diversity.Rdata")

###Annotation with the IDs
file_ids = rownames(PAAD_repertoire_diversity)
annotation = TCGAtranslateID(file_ids)
annotation$patient_barcode = substr(annotation$submitter_id,1,12)
write.csv(annotation,"Data/PAAD/annotation_PAAD.csv") ###To mark the non pancreas based on http://clincancerres.aacrjournals.org/content/clincanres/early/2018/05/08/1078-0432.CCR-18-0290.full.pdf
annotation<-read.csv("Data/PAAD/annotation_PAAD.csv")
id<-match(rownames(PAAD_repertoire_diversity),annotation$file_id)
PAAD_repertoire_diversity$TCGA_sample<-annotation$submitter_id[id]
PAAD_repertoire_diversity$Tumor_type<-annotation$tumor_type[id]

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
clinical.drug<-read.csv("Data/PAAD/Clinical/clinical_drug_PAAD.csv")
##Clinical.radiation
clinical.radiation<-clinical.radiation[!duplicated(clinical.radiation),]
clinical.radiation<-clinical.radiation[,-c(7,12:13,15)]
write.csv(clinical.radiation,"Data/PAAD/Clinical/clinical_radiation_PAAD.csv")
clinical.radiation<-read.csv("Data/PAAD/Clinical/clinical_radiation_PAAD.csv")
##Clinical.patient
clinical.patient<-clinical.patient[!duplicated(clinical.patient),]
clinical.patient$tobacco_smoking_history_master<-ifelse(clinical.patient$tobacco_smoking_history==1,"Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime)",
                                            ifelse(clinical.patient$tobacco_smoking_history==2,"Current smoker (includes daily smokers and non-daily smokers or occasional smokers)",
                                            ifelse(clinical.patient$tobacco_smoking_history==3,"Current reformed smoker for > 15 years (greater than 15 years)",
                                            ifelse(clinical.patient$tobacco_smoking_history==4,"Current reformed smoker for ≤15 years (less than or equal to 15 years)",
                                            ifelse(clinical.patient$tobacco_smoking_history==5,"Current reformed smoker, duration not specified",NA)))))
clinical.patient<-clinical.patient[,-c(2,72,75:80)]
write.csv(clinical.patient,"Data/PAAD/Clinical/clinical_patient_PAAD.csv")
clinical.patient<-read.csv("Data/PAAD/Clinical/clinical_patient_PAAD.csv")
##CLinical.stage has the same information than the one in clinical.patient

##clinical.new_tumor_event
clinical.new_tumor_event<-clinical.new_tumor_event[!duplicated(clinical.new_tumor_event),]
write.csv(clinical.new_tumor_event,"Data/PAAD/Clinical/clinical_new_tumor_event_PAAD.csv")
clinical.new_tumor_event<-read.csv("Data/PAAD/Clinical/clinical_new_tumor_event_PAAD.csv")

##CLinical.follow.up
clinical.folow_up<-read.csv("Data/PAAD/Clinical/Clinical_follow_up_fromSuppTableS1_PAAD.csv")
write.csv(clinical.folow_up,"Data/PAAD/Clinical/clinical_folow_up_PAAD.csv")
clinical.folow_up<-read.csv("Data/PAAD/Clinical/clinical_folow_up_PAAD.csv")

#biospecimen.sample
biospecimen.slide<-biospecimen.slide[!duplicated(biospecimen.slide),]
biospecimen.slide<-biospecimen.slide[,-c(2,4,10,13,15)]
write.csv(biospecimen.slide,"Data/PAAD/Clinical/biospecimen_slide_PAAD.csv")
biospecimen.slide<-read.csv("Data/PAAD/Clinical/biospecimen_slide_PAAD.csv")

##Subtypes of PAAD
paad.subtype <- TCGAquery_subtype(tumor = "paad")
write.csv(paad.subtype,"Data/PAAD/Clinical/paad.subtype.csv")
paad.subtype<-read.csv("Data/PAAD/Clinical/paad.subtype.csv")

###Read xCELL data
xCell.data<-read.table("Data/xCELL/xCell_PAAD_rsem_xCell_0643041519.txt",header = T,sep="\t")
xcell.pvalue<-read.table("Data/xCELL/xCell_PAAD_rsem_xCell_0643041519.pvals.txt",header=T,sep="\t")
rownames(xCell.data)<-xCell.data$X
rownames(xcell.pvalue)<-xcell.pvalue$X
xCell.data<-xCell.data[,-1]
xcell.pvalue<-xcell.pvalue[,-1]
xx<-str_replace(colnames(xCell.data),"\\.","-")
xx<-str_replace(xx,"\\.","-")
xx<-str_replace(xx,"\\.","-")
colnames(xCell.data)<-xx
colnames(xcell.pvalue)<-xx
id.xcell<-match(substr(PAAD_repertoire_diversity$TCGA_sample,1,15),colnames(xCell.data))
xCell.data.PAAD<-xCell.data[,na.omit(id.xcell)] ##181 
xCell.pvalue.PAAD<-xcell.pvalue[,na.omit(id.xcell)] ##181

PAAD.repertoire.diversity<-PAAD_repertoire_diversity
##Adding the tumor_type_2categ variable 
PAAD.repertoire.diversity$Tumor_type_2categ<-ifelse(PAAD.repertoire.diversity$Tumor_type=="Tumor_pancreas","Tumor_pancreas",
                                                    ifelse(PAAD.repertoire.diversity$Tumor_type=="Solid_tissue_normal","normal_pseudonormal_pancreas",
                                                           ifelse(PAAD.repertoire.diversity$Tumor_type=="Adjacent_normal_pancreas","normal_pseudonormal_pancreas",
                                                                  ifelse(PAAD.repertoire.diversity$Tumor_type=="Pseudonormal (<1% neoplastic cellularity)","normal_pseudonormal_pancreas",NA))))
PAAD.repertoire.diversity$Tumor_type_2categ<-as.factor(PAAD.repertoire.diversity$Tumor_type_2categ)
PAAD.repertoire.diversity<-PAAD.repertoire.diversity[which(is.na(PAAD.repertoire.diversity$Tumor_type_2categ)==F),]

##Adding the tumor_type_3categ variable 
PAAD.repertoire.diversity$Tumor_type_3categ<-ifelse(PAAD.repertoire.diversity$Tumor_type=="Tumor_pancreas","Tumor_pancreas",
                                                    ifelse(PAAD.repertoire.diversity$Tumor_type=="Solid_tissue_normal","normal_pancreas",
                                                           ifelse(PAAD.repertoire.diversity$Tumor_type=="Adjacent_normal_pancreas","normal_pancreas",
                                                                  ifelse(PAAD.repertoire.diversity$Tumor_type=="Pseudonormal (<1% neoplastic cellularity)","pseudonormal_pancreas",NA))))
PAAD.repertoire.diversity$Tumor_type_3categ<-as.factor(PAAD.repertoire.diversity$Tumor_type_3categ)
PAAD.repertoire.diversity<-PAAD.repertoire.diversity[which(is.na(PAAD.repertoire.diversity$Tumor_type_3categ)==F),]


##Adding the tumor_type_4categ variable
id<-match(substr(PAAD.repertoire.diversity$TCGA_sample,1,12),clinical.patient$bcr_patient_barcode)
PAAD.repertoire.diversity$Tumor_type_4categ<-ifelse(PAAD.repertoire.diversity$Tumor_type=="Tumor_pancreas" & 
                                                    clinical.patient$histological_type[id]=="Pancreas-Adenocarcinoma Ductal Type","PDAC",
                                                    ifelse(PAAD.repertoire.diversity$Tumor_type=="Tumor_pancreas" & 
                                                             clinical.patient$histological_type[id]=="Pancreas-Adenocarcinoma-Other Subtype","PDAC",
                                                    ifelse(PAAD.repertoire.diversity$Tumor_type=="Solid_tissue_normal","normal_pancreas",
                                                           ifelse(PAAD.repertoire.diversity$Tumor_type=="Adjacent_normal_pancreas","normal_pancreas",
                                                                  ifelse(PAAD.repertoire.diversity$Tumor_type=="Pseudonormal (<1% neoplastic cellularity)","pseudonormal_pancreas",NA)))))
PAAD.repertoire.diversity$Tumor_type_4categ<-as.factor(PAAD.repertoire.diversity$Tumor_type_4categ)

save(data_merge,PAAD.repertoire.diversity,xCell.data.PAAD,xCell.pvalue.PAAD,paad.subtype,clinical.drug,clinical.patient,clinical.radiation,clinical.new_tumor_event,clinical.folow_up,biospecimen.slide,annotation,
     file="Data/PAAD/PAAD_FullData.Rdata")



########################################################################
### Perform the subsampled data to obtain the diversity measures   ####
######################################################################
load("Data/PAAD/PAAD_FullData.Rdata")

##Diversity measures
Obtain_diversity<-function(data,PAAD.repertoire.diversity,chainType,sample_order){
  sample<-rownames(PAAD.repertoire.diversity)
  entropy<-NULL
  data<-data[which(data$chainType==chainType),]
  for (i in 1:length(sample)){
    print(i)
    data_sample_unique<-data[which(data$sample==sample[i]),]
    data_subsample<-data_sample_unique[sample(dim(data_sample_unique)[1], round(dim(data_sample_unique)[1]*sample_order)),]
    clones_sample<-data_subsample[,"V_J_lenghCDR3_CloneId"]
  
    fi<-as.numeric(table(clones_sample))/length(clones_sample)
    hi<-fi*log2(fi)
    entropy[i]=-sum(hi)
  }
  return(entropy)
}

chainType="IGH"
for (i in 1:5){
  assign(paste0("entropy_",i,"_",chainType),Obtain_diversity(data_merge,PAAD.repertoire.diversity,chainType,0.2))
}
PAAD.repertoire.diversity[,paste0("entropy_20_",chainType)]<-rowMeans(cbind(get(paste0("entropy_1_",chainType)),
                                                                                   get(paste0("entropy_2_",chainType)),
                                                                                   get(paste0("entropy_3_",chainType)),
                                                                                   get(paste0("entropy_4_",chainType)),
                                                                                   get(paste0("entropy_5_",chainType))))

for (i in 1:5){
  assign(paste0("entropy_",i,"_",chainType),Obtain_diversity(data_merge,PAAD.repertoire.diversity,chainType,0.4))
}
PAAD.repertoire.diversity[,paste0("entropy_40_",chainType)]<-rowMeans(cbind(get(paste0("entropy_1_",chainType)),
                                                                              get(paste0("entropy_2_",chainType)),
                                                                              get(paste0("entropy_3_",chainType)),
                                                                              get(paste0("entropy_4_",chainType)),
                                                                              get(paste0("entropy_5_",chainType))))

for (i in 1:5){
  assign(paste0("entropy_",i,"_",chainType),Obtain_diversity(data_merge,PAAD.repertoire.diversity,chainType,0.6))
}
PAAD.repertoire.diversity[,paste0("entropy_60_",chainType)]<-rowMeans(cbind(get(paste0("entropy_1_",chainType)),
                                                                              get(paste0("entropy_2_",chainType)),
                                                                              get(paste0("entropy_3_",chainType)),
                                                                              get(paste0("entropy_4_",chainType)),
                                                                              get(paste0("entropy_5_",chainType))))


for (i in 1:5){
  assign(paste0("entropy_",i,"_",chainType),Obtain_diversity(data_merge,PAAD.repertoire.diversity,chainType,0.8))
}
PAAD.repertoire.diversity[,paste0("entropy_80_",chainType)]<-rowMeans(cbind(get(paste0("entropy_1_",chainType)),
                                                                              get(paste0("entropy_2_",chainType)),
                                                                              get(paste0("entropy_3_",chainType)),
                                                                              get(paste0("entropy_4_",chainType)),
                                                                              get(paste0("entropy_5_",chainType))))

save(data_merge,PAAD.repertoire.diversity,xCell.data.PAAD,xCell.pvalue.PAAD,paad.subtype,clinical.drug,clinical.patient,clinical.radiation,clinical.new_tumor_event,clinical.folow_up,biospecimen.slide,annotation,
     file="Data/PAAD/PAAD_FullData_2.Rdata")


