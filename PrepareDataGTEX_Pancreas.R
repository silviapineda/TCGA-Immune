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

#####Filter by CDR3
data<-read.csv("Data/GTEx/Pancreas/GTEX_data_pancreas.csv")
data_full_cdr3<-data[which(data$nSeqCDR3!=""),] #79600
data_full_cdr3$CDR3_length<-nchar(as.character(data_full_cdr3$nSeqCDR3)) 
data_full_cdr3<-data_full_cdr3[which(data_full_cdr3$CDR3_length!=3),] #79600

##Obtain the ID for the clone call
data_full_cdr3$seqID<-seq(1,nrow(data_full_cdr3))
data_full_cdr3$SEQUENCE_ID<-paste(data_full_cdr3$sample,data_full_cdr3$seqID,sep="_")

##Variable to build clones
data_full_cdr3$V_J_lenghCDR3 = paste(data_full_cdr3$bestVGene,data_full_cdr3$bestJGene,data_full_cdr3$CDR3_length,sep="_")

###Chain
data_full_cdr3$chainType<-substr(data_full_cdr3$bestVGene,1,3)
data_full_cdr3_Ig<-data_full_cdr3[which(data_full_cdr3$chainType=="IGH" |
                                          data_full_cdr3$chainType=="IGK" |
                                          data_full_cdr3$chainType=="IGL"),]
data_full_cdr3_TCR<-data_full_cdr3[which(data_full_cdr3$chainType=="TRA" | 
                                           data_full_cdr3$chainType=="TRB" |
                                           data_full_cdr3$chainType=="TRD" | 
                                           data_full_cdr3$chainType=="TRG"),]


###save the data to call the clones by all samples using the nucleotides.py
data_clonesInference_Ig<-data_full_cdr3_Ig[,c("SEQUENCE_ID","sample_SRR","nSeqCDR3","CDR3_length","bestVGene","bestJGene","V_J_lenghCDR3")]
data_clonesInference_TCR<-data_full_cdr3_TCR[,c("SEQUENCE_ID","sample_SRR","nSeqCDR3","CDR3_length","bestVGene","bestJGene","V_J_lenghCDR3")]
write.table(data_clonesInference_Ig,file="Data/GTEx/Pancreas/data_for_cloneInfered_Ig_GTEX_PAAD.txt",row.names = F,sep="\t")
write.table(data_clonesInference_TCR,file="Data/GTEx/Pancreas/data_for_cloneInfered_TCR_GTEX_PAAD.txt",row.names = F,sep="\t")

### After passing the nucleotides.py
##Read the clones and merge with the data
nucleotides_Ig<-read.csv("Data/GTEx/Pancreas/ClonesInfered_Ig_GTEX_PAAD.csv")
nucleotides_TCR<-read.csv("Data/GTEx/Pancreas/ClonesInfered_TCR_GTEX_PAAD.csv")
nucleotides<-rbind(nucleotides_Ig,nucleotides_TCR)
data_merge<-merge(data_full_cdr3,nucleotides[,c("SEQUENCE_ID","CloneId")],by=c("SEQUENCE_ID"))

##Clones per chain
data_merge$V_J_lenghCDR3_CloneId = paste(data_merge$V_J_lenghCDR3,data_merge$CloneId,sep="_")

clones_count<- unique(data_merge[,c("sample_SRR","V_J_lenghCDR3_CloneId","chainType")])
clones<-data.frame(cbind(table(clones_count$sample_SRR,clones_count$chainType)))
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
  data_sample_unique<-data_merge[which(data_merge$sample_SRR==sample[i]),]
  clones_sample<-data_sample_unique[,"V_J_lenghCDR3_CloneId"]
  clones_sample_IGH<-data_sample_unique[which(data_sample_unique$chainType=="IGH"),"V_J_lenghCDR3_CloneId"]
  clones_sample_IGK<-data_sample_unique[which(data_sample_unique$chainType=="IGK"),"V_J_lenghCDR3_CloneId"]
  clones_sample_IGL<-data_sample_unique[which(data_sample_unique$chainType=="IGL"),"V_J_lenghCDR3_CloneId"]
  clones_sample_TRA<-data_sample_unique[which(data_sample_unique$chainType=="TRA"),"V_J_lenghCDR3_CloneId"]
  clones_sample_TRB<-data_sample_unique[which(data_sample_unique$chainType=="TRB"),"V_J_lenghCDR3_CloneId"]
  clones_sample_TRD<-data_sample_unique[which(data_sample_unique$chainType=="TRD"),"V_J_lenghCDR3_CloneId"]
  clones_sample_TRG<-data_sample_unique[which(data_sample_unique$chainType=="TRG"),"V_J_lenghCDR3_CloneId"]
  
  #To write file to run with Recon
  write.delim(data.frame(table(table(clones_sample_IGH))),file=paste("clones_sample_IGH_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  write.delim(data.frame(table(table(clones_sample_IGK))),file=paste("clones_sample_IGK_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  write.delim(data.frame(table(table(clones_sample_IGL))),file=paste("clones_sample_IGL_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  write.delim(data.frame(table(table(clones_sample_TRA))),file=paste("clones_sample_TRA_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  write.delim(data.frame(table(table(clones_sample_TRB))),file=paste("clones_sample_TRB_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  write.delim(data.frame(table(table(clones_sample_TRD))),file=paste("clones_sample_TRD_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  write.delim(data.frame(table(table(clones_sample_TRG))),file=paste("clones_sample_TRG_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  
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
recon<-read.table("Data/GTEx/Pancreas/Recon/test_D_number_table.txt",header=T)
chain<-substr(recon$sample_name,70,72)
sample<-substr(recon$sample_name,74,83)
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
colnames(diversity)[29:35]<-c("cdr3_length_IGH","cdr3_length_IGK","cdr3_length_IGL","cdr3_length_TRA","cdr3_length_TRB",
                              "cdr3_length_TRD","cdr3_length_TRG")
data_merge_pancreas<-data_merge
diversity_pancreas<-diversity
save(data_merge_pancreas,diversity_pancreas,file="Data/GTEx/Pancreas/GTEx_Pancreas_RepertoireResults_diversity.Rdata")


####MIXED WITH THE DATA EXTRACTED FROM INES
load("Data/GTEx/Pancreas/GTEx_Pancreas_RepertoireResults_diversity.Rdata")

####### The data needs to be normalized by the unmapped reads 
read_count <- table(data_merge_pancreas$sample_SRR)
read_count_chain <- table(data_merge_pancreas$sample_SRR, data_merge_pancreas$chainType)
reads <- data.frame(cbind(read_count,read_count_chain))

totalReads<-read.table("Ines//total_reads.csv",sep=";") ##We need to extract this number from the MIXCR report with the python script
totalReads$V1<-unlist(strsplit(as.character(totalReads$V1), "\\_"))
id<-match(totalReads$V1,rownames(reads))
reads$totalSeqReads<-totalReads[id,2]

####Normalize the nuber of reads
reads$IG_expression<-(reads$IGH+reads$IGK+reads$IGL)/reads$totalSeqReads
reads$IGH_expression<-reads$IGH/reads$totalSeqReads
reads$IGK_expression<-reads$IGK/reads$totalSeqReads
reads$IGL_expression<-reads$IGL/reads$totalSeqReads

reads$T_expression<-(reads$TRA+reads$TRB+reads$TRD+reads$TRG)/reads$totalSeqReads
reads$TRA_expression<-reads$TRA/reads$totalSeqReads
reads$TRB_expression<-reads$TRB/reads$totalSeqReads
reads$TRD_expression<-reads$TRD/reads$totalSeqReads
reads$TRG_expression<-reads$TRG/reads$totalSeqReads
###Ratio
reads$Alpha_Beta_ratio_expression<-(reads$TRA_expression+reads$TRB_expression)/reads$T_expression
reads$KappaLambda_ratio_expression <- (reads$IGK_expression / reads$IGL_expression)

Pancreas.repertoire.diversity<-cbind(reads[,-9],diversity_pancreas)
save(data_merge_pancreas,Pancreas.repertoire.diversity,file="Data/GTEx/Pancreas/GTEx_FullData.Rdata")