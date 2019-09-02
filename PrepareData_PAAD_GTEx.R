rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. 
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Preparing Pancreas TCGA and Validation together 
###         
###
### Author: Silvia Pineda
### Date: August, 2019
############################################################################################
library(pgirmess)
library(stringr)
library(SciViews)

setwd("~/TCGA-Immune/")

load("Data/GTEx/Pancreas/GTEX_data_pancreas.Rdata")
data_gtex<-data
data_gtex$dataset<-"GTEx_pancreas"
data_gtex<-data_gtex[,c(2,1,3:ncol(data_gtex))]
colnames(data_gtex)[2]<-"readSequence"

load("Data/PAAD/PAAD_MIXCR_results.Rdata")
data$dataset<-"TCGA_PAAD"

data<-rbind(data_gtex,data)

###Chain
data$chainType<-ifelse(data$bestVGene!="", substr(data$bestVGene,1,3),
                       ifelse(data$bestVGene=="",substr(data$bestJGene,1,3),NA))

##Reads per chain
read_count <- table(data$sample)
read_count_chain <- table(data$sample, data$chainType)
reads <- data.frame(cbind(read_count,read_count_chain))

####### The data needs to be normalized by the unmapped reads 
totalReads_PAAD<-read.table("Data/PAAD/MIXCR_PAAD/total_reads.txt",sep=";") ##We need to extract this number from the MIXCR report with the python script
totalReads_gtex<-read.table("Data/GTEx/Pancreas/MIXCR/report/total_reads.txt",sep=";") ##We need to extract this number from the MIXCR report with the python script
totalReads_gtex$V1<-substr(totalReads_gtex$V1,4,13)
totalReads_gtex$V1<-unlist(strsplit(totalReads_gtex$V1, "_"))
totalReads<-rbind(totalReads_PAAD,totalReads_gtex)
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
data_full_cdr3<-data[which(data$nSeqCDR3!=""),] #1082293
data_full_cdr3$CDR3_length<-nchar(as.character(data_full_cdr3$nSeqCDR3)) 
data_full_cdr3<-data_full_cdr3[which(data_full_cdr3$CDR3_length>3),] 

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
data_clonesInference_Ig<-data_full_cdr3_Ig[,c("SEQUENCE_ID","sample","nSeqCDR3","aaSeqCDR3","CDR3_length","bestVGene","bestJGene","V_J_lenghCDR3")]
data_clonesInference_TCR<-data_full_cdr3_TCR[,c("SEQUENCE_ID","sample","nSeqCDR3","aaSeqCDR3","CDR3_length","bestVGene","bestJGene","V_J_lenghCDR3")]
write.table(data_clonesInference_Ig,file="Data/data_for_cloneInfered_Ig_PAAD_GTEx.txt",row.names = F,sep="\t")
write.table(data_clonesInference_TCR,file="Data/data_for_cloneInfered_TCR_PAAD_GTEx.txt",row.names = F,sep="\t")

### After passing the nucleotides.py
##Read the clones and merge with the data
nucleotides_Ig<-read.csv("Data/PAAD_GTEx/ClonesInfered_Ig_PAAD_GTEx.csv")
nucleotides_TCR<-read.csv("Data/PAAD_GTEx//ClonesInfered_TCR_PAAD_GTEx.csv")
nucleotides<-rbind(nucleotides_Ig,nucleotides_TCR)
data_merge<-merge(data_full_cdr3,nucleotides[,c("SEQUENCE_ID","CloneId")],by=c("SEQUENCE_ID"))

####Reads per chain only for clonal information
##Reads per chain
read_count <- table(data_merge$sample)
read_count_chain <- table(data_merge$sample, data_merge$chainType)
reads_filter <- data.frame(cbind(read_count,read_count_chain))
####### The data needs to be normalized by the unmapped reads 
totalReads_PAAD<-read.table("Data/PAAD/MIXCR_PAAD/total_reads.txt",sep=";") ##We need to extract this number from the MIXCR report with the python script
totalReads_gtex<-read.table("Data/GTEx/Pancreas/MIXCR/report/total_reads.txt",sep=";") ##We need to extract this number from the MIXCR report with the python script
totalReads_gtex$V1<-substr(totalReads_gtex$V1,4,13)
totalReads_gtex$V1<-unlist(strsplit(totalReads_gtex$V1, "_"))
totalReads<-rbind(totalReads_PAAD,totalReads_gtex)
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
  # write.delim(data.frame(table(table(clones_sample_IGH))),file=paste("Data/PAAD_GTEx//RECON/clones_sample_IGH_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.delim(data.frame(table(table(clones_sample_IGK))),file=paste("Data/PAAD_GTEx/RECON/clones_sample_IGK_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.delim(data.frame(table(table(clones_sample_IGL))),file=paste("Data/PAAD_GTEx/RECON/clones_sample_IGL_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.delim(data.frame(table(table(clones_sample_TRA))),file=paste("Data/PAAD_GTEx/RECON/clones_sample_TRA_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.delim(data.frame(table(table(clones_sample_TRB))),file=paste("Data/PAAD_GTEx/RECON/clones_sample_TRB_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.delim(data.frame(table(table(clones_sample_TRD))),file=paste("Data/PAAD_GTEx/RECON/clones_sample_TRD_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # write.delim(data.frame(table(table(clones_sample_TRG))),file=paste("Data/PAAD_GTEx/RECON/clones_sample_TRG_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  # 
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
recon<-read.table("Data/PAAD_GTEx/RECON/test_D_number_table.txt",header=T)
chain<-substr(recon$sample_name,66,68)
sample<-substr(recon$sample_name,70,105)
sample<-unlist(strsplit(sample, "\\."))
sample<-sample[which(nchar(sample)>7)]
##0.0D is species richness (Number of clones)
##Entropy is ln(1.0.D)
diversity<-as.data.frame(diversity)
chain_list<-unique(chain)
for(i in chain_list){
  recon_chain<-recon[which(chain==i),]
  sample_chain<-substr(recon_chain$sample_name,70,105)
  sample_chain<-unlist(strsplit(sample_chain, "\\."))
  sample_chain<-sample_chain[which(nchar(sample_chain)>7)]
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
colnames(diversity)[29:35]<-c("cdr3_length_IGH","cdr3_length_IGK","cdr3_length_IGL","cdr3_length_TRA","cdr3_length_TRB",
                              "cdr3_length_TRD","cdr3_length_TRG")


PAAD.GTEx.repertoire.diversity<-cbind(reads,reads_filter,diversity)

###Para que no se aplaste el data_merge llamando a los otros FULLDATA
data_merge_PAAD_GTEX<-data_merge

load("Data/GTEx/Pancreas/GTEx_FullData.Rdata")
load("Data/PAAD/PAAD_FullData.Rdata")

#PAAD.repertoire.diversity<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type_3categ=="Tumor_pancreas"),]
id.paad<-match(rownames(PAAD.repertoire.diversity),rownames(PAAD.GTEx.repertoire.diversity))
id.gtex<-match(rownames(Pancreas.repertoire.diversity),rownames(PAAD.GTEx.repertoire.diversity))

##Outcome
PAAD.GTEx.repertoire.diversity$outcome<-NA
PAAD.GTEx.repertoire.diversity$outcome[na.omit(id.paad)]<-as.character(PAAD.repertoire.diversity$Tumor_type_3categ)
PAAD.GTEx.repertoire.diversity$outcome[na.omit(id.gtex)]<-c("Normal-pancreas (GTEx)")

PAAD.GTEx.repertoire.diversity$outcome<-ifelse(PAAD.GTEx.repertoire.diversity$outcome=="normal_pancreas","normal-pancreas (TCGA)",
                                        ifelse(PAAD.GTEx.repertoire.diversity$outcome=="pseudonormal_pancreas","pseudonormal-pancreas (TCGA)",
                                        ifelse(PAAD.GTEx.repertoire.diversity$outcome=="Tumor_pancreas","tumor-pancreas (TCGA)",
                                        ifelse(PAAD.GTEx.repertoire.diversity$outcome=="Normal-pancreas (GTEx)","normal-pancreas (GTEx)",NA))))

##Annotation sample
PAAD.GTEx.repertoire.diversity$sample2<-NA
PAAD.GTEx.repertoire.diversity$sample2[na.omit(id.paad)]<-as.character(substr(PAAD.repertoire.diversity$TCGA_sample,1,12))
PAAD.GTEx.repertoire.diversity$sample2[na.omit(id.gtex)]<-as.character(Pancreas.repertoire.diversity$SUBJID)

##Sex
id.paad<-match(clinical.patient$bcr_patient_barcode,PAAD.GTEx.repertoire.diversity$sample2)
PAAD.GTEx.repertoire.diversity$sex<-NA
PAAD.GTEx.repertoire.diversity$sex[na.omit(id.paad)]<-as.character(clinical.patient$gender[which(na.omit(id.paad)!=F)])
id.gtex<-match(annotation_gtex_pancreas$SUBJID,PAAD.GTEx.repertoire.diversity$sample2)
PAAD.GTEx.repertoire.diversity$sex[na.omit(id.gtex)]<-as.character(annotation_gtex_pancreas$SEX[which(na.omit(id.gtex)!=F)])

PAAD.GTEx.repertoire.diversity$sex<-ifelse(PAAD.GTEx.repertoire.diversity$sex=="Female" | PAAD.GTEx.repertoire.diversity$sex=="FEMALE","Female",
                                    ifelse(PAAD.GTEx.repertoire.diversity$sex=="Male"| PAAD.GTEx.repertoire.diversity$sex=="MALE","Male",NA))

##Age
PAAD.GTEx.repertoire.diversity$age<-NA
PAAD.GTEx.repertoire.diversity$age[na.omit(id.paad)]<-clinical.patient$age_at_initial_pathologic_diagnosis[which(na.omit(id.paad)!=F)]
PAAD.GTEx.repertoire.diversity$age[na.omit(id.gtex)]<-annotation_gtex_pancreas$AGE[which(na.omit(id.gtex)!=F)]

##Race
PAAD.GTEx.repertoire.diversity$race<-NA
PAAD.GTEx.repertoire.diversity$race[na.omit(id.paad)]<-as.character(clinical.patient$race_list[which(na.omit(id.paad)!=F)])
PAAD.GTEx.repertoire.diversity$race[na.omit(id.gtex)]<-as.character(annotation_gtex_pancreas$RACE[which(na.omit(id.gtex)!=F)])

PAAD.GTEx.repertoire.diversity$race<-ifelse(PAAD.GTEx.repertoire.diversity$race=="Asian" | PAAD.GTEx.repertoire.diversity$race=="ASIAN","Asian",
                                           ifelse(PAAD.GTEx.repertoire.diversity$race=="BLACK OR AFRICAN AMERICAN"| PAAD.GTEx.repertoire.diversity$race=="Black/AA","Black/AA",
                                                  ifelse(PAAD.GTEx.repertoire.diversity$race=="White"| PAAD.GTEx.repertoire.diversity$race=="WHITE","White",NA)))

PAAD.GTEx.repertoire.diversity<-PAAD.GTEx.repertoire.diversity[which(is.na(PAAD.GTEx.repertoire.diversity$outcome)==F),]


PAAD.GTEx.repertoire.diversity$sample<-rownames(PAAD.GTEx.repertoire.diversity)

data_merge<-data_merge_PAAD_GTEX
save(data_merge,PAAD.GTEx.repertoire.diversity,file="Data/PAAD_GTEx/PAAD_GTEx_FullData.Rdata")


######################
### Downsampling ####
#####################
load("Data/PAAD_GTEx/PAAD_GTEx_FullData.Rdata")

##Filter by clones > 200
PAAD.GTEx.repertoire.diversity$Ig_clones<-PAAD.GTEx.repertoire.diversity$clones_IGH+PAAD.GTEx.repertoire.diversity$clones_IGK+
  PAAD.GTEx.repertoire.diversity$clones_IGL
PAAD.GTEx.repertoire.diversity_Igclones<-PAAD.GTEx.repertoire.diversity[which(PAAD.GTEx.repertoire.diversity$Ig_clones>200),]

###Dowwn for IG
sample<-PAAD.GTEx.repertoire.diversity_Igclones$sample
id<-match(data_clonesInference_Ig$sample,sample)
data_clonesInference_Ig_qc<-data_clonesInference_Ig[which(is.na(id)==F),]
min<-min(table(data_clonesInference_Ig_qc$sample))  ##311 reads
set.seed(27)
for (j in 1:5){
  data_clonesInference_Ig_down<-NULL
  for (i in 1:length(sample)){
    print(i)
    data_clonesInference_Ig_sample<-data_clonesInference_Ig_qc[which(data_clonesInference_Ig_qc$sample==sample[i]),]
    data_clonesInference_Ig_down<-rbind(data_clonesInference_Ig_down,data_clonesInference_Ig_sample[sample(nrow(data_clonesInference_Ig_sample),min),])
  }
  write.table(data_clonesInference_Ig_down,file=paste0("Data/PAAD_GTEx/data_for_cloneInfered_Ig_PAAD_GTEx_down",j,".txt"),row.names = F,sep="\t")
}

### After passing the nucleotides.py
##Read the clones and merge with the data
for (j in 1:5){
  nucleotides_Ig_down<-read.csv(paste0("Data/PAAD_GTEx//ClonesInfered_Ig_PAAD_GTEx_down",j,".csv"))
  data_merge_down1<-merge(data_full_cdr3,nucleotides_Ig_down[,c("SEQUENCE_ID","CloneId")],by=c("SEQUENCE_ID"))
  
  ##Clones per chain
  data_merge_down1$V_J_lenghCDR3_CloneId = paste(data_merge_down1$V_J_lenghCDR3,data_merge_down1$CloneId,sep="_")
  
  clones_count<- unique(data_merge_down1[,c("sample","V_J_lenghCDR3_CloneId","chainType")])
  clones<-data.frame(cbind(table(clones_count$sample,clones_count$chainType)))
  colnames(clones)<-c("clones_IGH","clones_IGK","clones_IGL")
  
  ##Diversity measures
  sample<-rownames(clones)
  entropy_IGH<-NULL
  entropy_IGK<-NULL
  entropy_IGL<-NULL
  for (i in 1:length(sample)){
    print(i)
    data_sample_unique<-data_merge_down1[which(data_merge_down1$sample==sample[i]),]
    clones_sample<-data_sample_unique[,"V_J_lenghCDR3_CloneId"]
    clones_sample_IGH<-data_sample_unique[which(data_sample_unique$chainType=="IGH"),"V_J_lenghCDR3_CloneId"]
    clones_sample_IGK<-data_sample_unique[which(data_sample_unique$chainType=="IGK"),"V_J_lenghCDR3_CloneId"]
    clones_sample_IGL<-data_sample_unique[which(data_sample_unique$chainType=="IGL"),"V_J_lenghCDR3_CloneId"]
    
    #To write file to run with Recon
    #write.delim(data.frame(table(table(clones_sample_IGH))),file=paste("Data/PAAD_GTEx//RECON/Down/clones_sample_IGH_down_",j,sample[i],".txt",sep=""),sep="\t",col.names=F)
    #write.delim(data.frame(table(table(clones_sample_IGK))),file=paste("Data/PAAD_GTEx/RECON/Down/clones_sample_IGK_down_",j,sample[i],".txt",sep=""),sep="\t",col.names=F)
    #write.delim(data.frame(table(table(clones_sample_IGL))),file=paste("Data/PAAD_GTEx/RECON/Down/clones_sample_IGL_down_",j,sample[i],".txt",sep=""),sep="\t",col.names=F)
    
    fi_IGH<-as.numeric(table(clones_sample_IGH))/length(clones_sample_IGH)
    fi_IGK<-as.numeric(table(clones_sample_IGK))/length(clones_sample_IGK)
    fi_IGL<-as.numeric(table(clones_sample_IGL))/length(clones_sample_IGL)
    
    hi_IGH<-fi_IGH*log2(fi_IGH)
    hi_IGK<-fi_IGK*log2(fi_IGK)
    hi_IGL<-fi_IGL*log2(fi_IGL)
   
    entropy_IGH[i]=-sum(hi_IGH)
    entropy_IGK[i]=-sum(hi_IGK)
    entropy_IGL[i]=-sum(hi_IGL)
    
  }
  assign(paste0("diversity_down",j),cbind(clones,entropy_IGH,entropy_IGK,entropy_IGL))
  
  ####After runing recon
  recon<-read.table(paste0("Data/PAAD_GTEx/RECON/Down/test_D_number_table_down_",j,".txt"),header=T)
  chain<-substr(recon$sample_name,71,73)
  sample<-substr(recon$sample_name,81,116)
  sample<-unlist(strsplit(sample, "\\."))
  sample<-sample[which(nchar(sample)>7)]
  ##0.0D is species richness (Number of clones)
  ##Entropy is ln(1.0.D)
  chain_list<-unique(chain)
  for(i in chain_list){
    recon_chain<-recon[which(chain==i),]
    sample_chain<-substr(recon_chain$sample_name,81,116)
    sample_chain<-unlist(strsplit(sample_chain, "\\."))
    sample_chain<-sample_chain[which(nchar(sample_chain)>7)]
    id<-match(rownames(get(paste0("diversity_down",j))),sample_chain)
    clone_chain<-ifelse(is.na(id)==F,recon_chain[id,"est_0.0D"],0)
    assign(paste0("clones_recon_",i),clone_chain)
    entroy_chain<-ifelse(is.na(id)==F,ln(recon_chain[id,"est_1.0D"]),0)
    assign(paste0("entropy_recon_",i),entroy_chain)
  }
  assign(paste0("diversity_down",j),cbind(get(paste0("diversity_down",j)),clones_recon_IGH,clones_recon_IGK,clones_recon_IGL,entropy_recon_IGH,entropy_recon_IGK,entropy_recon_IGL))
}  


id<-match(rownames(diversity_down1),rownames(PAAD.GTEx.repertoire.diversity))

PAAD.GTEx.repertoire.diversity$clones_IGH_down<-NA
PAAD.GTEx.repertoire.diversity$clones_IGK_down<-NA
PAAD.GTEx.repertoire.diversity$clones_IGL_down<-NA
PAAD.GTEx.repertoire.diversity$entropy_IGH_down<-NA
PAAD.GTEx.repertoire.diversity$entropy_IGK_down<-NA
PAAD.GTEx.repertoire.diversity$entropy_IGL_down<-NA

PAAD.GTEx.repertoire.diversity$clones_IGH_down[na.omit(id)]<-rowMeans(cbind(diversity_down1$clones_IGH,diversity_down2$clones_IGH,diversity_down3$clones_IGH,
                                                               diversity_down4$clones_IGH,diversity_down5$clones_IGH))
PAAD.GTEx.repertoire.diversity$clones_IGK_down[na.omit(id)]<-rowMeans(cbind(diversity_down1$clones_IGK,diversity_down2$clones_IGK,diversity_down3$clones_IGK,
                                                                    diversity_down4$clones_IGK,diversity_down5$clones_IGK))
PAAD.GTEx.repertoire.diversity$clones_IGL_down[na.omit(id)]<-rowMeans(cbind(diversity_down1$clones_IGL,diversity_down2$clones_IGL,diversity_down3$clones_IGL,
                                                                            diversity_down4$clones_IGL,diversity_down5$clones_IGL))
PAAD.GTEx.repertoire.diversity$entropy_IGH_down[na.omit(id)]<-rowMeans(cbind(diversity_down1$entropy_IGH,diversity_down2$entropy_IGH,diversity_down3$entropy_IGH,
                                                                             diversity_down4$entropy_IGH,diversity_down5$entropy_IGH))
PAAD.GTEx.repertoire.diversity$entropy_IGK_down[na.omit(id)]<-rowMeans(cbind(diversity_down1$entropy_IGK,diversity_down2$entropy_IGK,diversity_down3$entropy_IGK,
                                                                             diversity_down4$entropy_IGK,diversity_down5$entropy_IGK))
PAAD.GTEx.repertoire.diversity$entropy_IGL_down[na.omit(id)]<-rowMeans(cbind(diversity_down1$entropy_IGL,diversity_down2$entropy_IGL,diversity_down3$entropy_IGL,
                                                                             diversity_down4$entropy_IGL,diversity_down5$entropy_IGL))


PAAD.GTEx.repertoire.diversity$clones_recon_IGH_down<-NA
PAAD.GTEx.repertoire.diversity$clones_recon_IGK_down<-NA
PAAD.GTEx.repertoire.diversity$clones_recon_IGL_down<-NA
PAAD.GTEx.repertoire.diversity$entropy_recon_IGH_down<-NA
PAAD.GTEx.repertoire.diversity$entropy_recon_IGK_down<-NA
PAAD.GTEx.repertoire.diversity$entropy_recon_IGL_down<-NA

PAAD.GTEx.repertoire.diversity$clones_recon_IGH_down[na.omit(id)]<-rowMeans(cbind(diversity_down1$clones_recon_IGH,diversity_down2$clones_recon_IGH,diversity_down3$clones_recon_IGH,
                                                                                  diversity_down4$clones_recon_IGH,diversity_down5$clones_recon_IGH))
PAAD.GTEx.repertoire.diversity$clones_recon_IGK_down[na.omit(id)]<-rowMeans(cbind(diversity_down1$clones_recon_IGK,diversity_down2$clones_IGK,diversity_down3$clones_recon_IGK,
                                                                            diversity_down4$clones_recon_IGK,diversity_down5$clones_recon_IGK))
PAAD.GTEx.repertoire.diversity$clones_recon_IGL_down[na.omit(id)]<-rowMeans(cbind(diversity_down1$clones_recon_IGL,diversity_down2$clones_recon_IGL,diversity_down3$clones_recon_IGL,
                                                                            diversity_down4$clones_recon_IGL,diversity_down5$clones_recon_IGL))
PAAD.GTEx.repertoire.diversity$entropy_recon_IGH_down[na.omit(id)]<-rowMeans(cbind(diversity_down1$entropy_recon_IGH,diversity_down2$entropy_recon_IGH,diversity_down3$entropy_recon_IGH,
                                                                             diversity_down4$entropy_recon_IGH,diversity_down5$entropy_recon_IGH))
PAAD.GTEx.repertoire.diversity$entropy_recon_IGK_down[na.omit(id)]<-rowMeans(cbind(diversity_down1$entropy_recon_IGK,diversity_down2$entropy_recon_IGK,diversity_down3$entropy_recon_IGK,
                                                                             diversity_down4$entropy_recon_IGK,diversity_down5$entropy_recon_IGK))
PAAD.GTEx.repertoire.diversity$entropy_recon_IGL_down[na.omit(id)]<-rowMeans(cbind(diversity_down1$entropy_recon_IGL,diversity_down2$entropy_recon_IGL,diversity_down3$entropy_recon_IGL,
                                                                             diversity_down4$entropy_recon_IGL,diversity_down5$entropy_recon_IGL))
save(data_merge,PAAD.GTEx.repertoire.diversity,file="Data/PAAD_GTEx/PAAD_GTEx_FullData.Rdata")
