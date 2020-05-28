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
files <- list.files("Data/GTEx/Blood/MIXCR/alignments/")
data<-c()
for(i in files) {
  cat(i, "\n")
  t <- read.delim(paste0("Data/GTEx/Blood/MIXCR/alignments/",i))
  t$sample<-substr(i, 2, nchar(i)-16)
  data <- rbind(data, t)
}

sra<-read.csv("Data/GTEx/Blood/SraRunTableBloodRNAseq.csv")
id_rna_blood<-match(data$sample,sra$Run)
data<-data[which(is.na(id_rna_blood)==F),] ##done in the server
save(data,file="Data/GTEx/Blood/GTEX_data_blood.Rdata")

#####
load("Data/GTEx/Blood/GTEX_data_blood.Rdata")


###Chain
data$chainType<-ifelse(data$bestVGene!="", substr(data$bestVGene,1,3),
                       ifelse(data$bestVGene=="",substr(data$bestJGene,1,3),NA))

###Add a column with the reads per chain for Ig and TRA
##Reads per chain
read_count <- table(data$sample)
read_count_chain <- table(data$sample, data$chainType)
reads <- data.frame(cbind(read_count,read_count_chain))

####### The data needs to be normalized by the unmapped reads 
totalReads<-read.table("Data/GTEx/Blood/MIXCR/report/total_reads.txt",sep=";") ##We need to extract this number from the MIXCR report with the python script
totalReads$V1<-substr(totalReads$V1,2,11)
totalReads$V1<-unlist(strsplit(totalReads$V1, "}"))
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

####Restrcition to CDR3 to obatin the clones
data_full_cdr3<-data[which(data$nSeqCDR3!=""),] #81277
data_full_cdr3$CDR3_length<-nchar(as.character(data_full_cdr3$nSeqCDR3)) 
data_full_cdr3<-data_full_cdr3[which(data_full_cdr3$CDR3_length!=3),] #81277

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
write.table(data_clonesInference_Ig,file="Data/GTEx/Blood/data_for_cloneInfered_Ig_GTEX_Blood.txt",row.names = F,sep="\t")
write.table(data_clonesInference_TCR,file="Data/GTEx/Blo/data_for_cloneInfered_TCR_GTEX_Blood.txt",row.names = F,sep="\t")

### After passing the nucleotides.py
##Read the clones and merge with the data
nucleotides_Ig<-read.csv("Data/GTEx/Blood//ClonesInfered_Ig_GTEx_Blood.csv")
nucleotides_TCR<-read.csv("Data/GTEx/Blood/ClonesInfered_TCR_GTEx_Blood.csv")
nucleotides<-rbind(nucleotides_Ig,nucleotides_TCR)
data_merge<-merge(data_full_cdr3,nucleotides[,c("SEQUENCE_ID","CloneId")],by=c("SEQUENCE_ID"))

####Reads per chain only for clonal information
##Reads per chain
read_count <- table(data_merge$sample)
read_count_chain <- table(data_merge$sample, data_merge$chainType)
reads_filter <- data.frame(cbind(read_count,read_count_chain))

####### The data needs to be normalized by the unmapped reads 
totalReads<-read.table("Data/GTEx/Blood/MIXCR/report/total_reads.txt",sep=";") ##We need to extract this number from the MIXCR report with the python script
totalReads$V1<-substr(totalReads$V1,2,11)
totalReads$V1<-unlist(strsplit(totalReads$V1, "}"))
id<-match(rownames(reads),totalReads$V1)
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
  # write.table(data.frame(table(table(clones_sample_IGH))),file=paste("clones_sample_IGH_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.table(data.frame(table(table(clones_sample_IGK))),file=paste("clones_sample_IGK_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.table(data.frame(table(table(clones_sample_IGL))),file=paste("clones_sample_IGL_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.table(data.frame(table(table(clones_sample_TRA))),file=paste("clones_sample_TRA_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.table(data.frame(table(table(clones_sample_TRB))),file=paste("clones_sample_TRB_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.table(data.frame(table(table(clones_sample_TRD))),file=paste("clones_sample_TRD_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.table(data.frame(table(table(clones_sample_TRG))),file=paste("clones_sample_TRG_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # # 
  
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

diversity<-cbind(reads,reads_filter,clones,entropy_IGH,entropy_IGK,entropy_IGL,entropy_TRA,entropy_TRB,entropy_TRD,entropy_TRG)

save(diversity,data_merge,"Data/GTEx/Blood/GTEx_Blood_RepertoireResults_diversity.Rdata")

####After runing recon
recon_IG<-read.table("Data/GTEx/Blood/RECON/test_D_number_table_IG.txt",header=T)
recon_TR<-read.table("Data/GTEx/Blood/RECON/test_D_number_table_TR.txt",header=T)

recon<-rbind(recon_IG,recon_TR)
chain<-substr(recon$sample_name,67,69)
sample<-substr(recon$sample_name,71,80)
sample<-unlist(strsplit(sample, "\\."))

##0.0D is species richness (Number of clones)
##Entropy is ln(1.0.D)
diversity<-as.data.frame(diversity)
chain_list<-unique(chain)
for(i in chain_list){
  recon_chain<-recon[which(chain==i),]
  id<-match(rownames(diversity),sample)
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
cdr3_length_IGH_2<-ifelse(is.na(id.IGH)==T,0,cdr3_length_IGH$x)
id.IGK<-match(rownames(diversity),cdr3_length_IGK$Group.1)
cdr3_length_IGK_2<-ifelse(is.na(id.IGK)==T,0,cdr3_length_IGK$x)
id.IGL<-match(rownames(diversity),cdr3_length_IGL$Group.1)
cdr3_length_IGL_2<-ifelse(is.na(id.IGL)==T,0,cdr3_length_IGL$x)
id.TRA<-match(rownames(diversity),cdr3_length_TRA$Group.1)
cdr3_length_TRA_2<-ifelse(is.na(id.TRA)==T,0,cdr3_length_TRA$x)
id.TRB<-match(rownames(diversity),cdr3_length_TRB$Group.1)
cdr3_length_TRB_2<-ifelse(is.na(id.TRB)==T,0,cdr3_length_TRB$x)
id.TRD<-match(rownames(diversity),cdr3_length_TRD$Group.1)
cdr3_length_TRD_2<-ifelse(is.na(id.TRD)==T,0,cdr3_length_TRD$x)
id.TRG<-match(rownames(diversity),cdr3_length_TRG$Group.1)
cdr3_length_TRG_2<-ifelse(is.na(id.TRG)==T,0,cdr3_length_TRG$x)

diversity<-cbind(diversity,cdr3_length_IGH_2,cdr3_length_IGK_2,cdr3_length_IGL_2,cdr3_length_TRA_2,cdr3_length_TRB_2,cdr3_length_TRD_2,cdr3_length_TRG_2)
colnames(diversity)[51:57]<-c("cdr3_length_IGH","cdr3_length_IGK","cdr3_length_IGL","cdr3_length_TRA","cdr3_length_TRB",
                              "cdr3_length_TRD","cdr3_length_TRG")
data_merge_blood<-data_merge
GTEX.blood.repertoire.diversity<-diversity
save(data_merge_blood,GTEX.blood.repertoire.diversity,file="Data/GTEx/Pancreas/GTEx_Blood_RepertoireResults_diversity.Rdata")

####Annotation 
sra<-read.csv("Data/GTEx/Blood/SraRunTableBloodRNAseq.csv")

###Samples that has passed QC for gene expression
samples<-read.csv("Data/GTEx/samples_QC_GTEX.csv")
id<-match(samples$x, sra$Sample_Name)
sra_qc<-sra[na.omit(id),]

##Match with GTEX blood diversity
id<-match(rownames(GTEX.blood.repertoire.diversity),sra_qc$Run)
GTEX.blood.repertoire.diversity<-GTEX.blood.repertoire.diversity[which(is.na(id)==F),] ##391
GTEX.blood.repertoire.diversity$SUBJID<-sra_qc$submitted_subject_id[na.omit(id)]

annotation_gtex<-read.csv("Data/GTEx/GTEX_annotation_phenotypes.csv")
id<-match(GTEX.blood.repertoire.diversity$SUBJID,annotation_gtex$SUBJID)
annotation_gtex_blood<-annotation_gtex[id,]

id<-match(data_merge_blood$sample,rownames(GTEX.blood.repertoire.diversity))
data_merge_blood2<-data_merge_blood[which(is.na(id)==F),]

save(data_merge_blood,GTEX.blood.repertoire.diversity,annotation_gtex_blood,file="Data/GTEx/Blood/GTEx_blood_FullData.Rdata")
