rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. Prepare GTEX data 
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: 
###         
###
### Author: Silvia Pineda
### Date: May, 2019
############################################################################################
library(pgirmess)
library(stringr)
library(SciViews)

setwd("~/TCGA-Immune/")
load("Data/GTEX/MIXCR/dataset_MIXCR_GTEX.Rdata") 

#####Filter by CDR3
data_full_cdr3<-alignmentData[which(alignmentData$nSeqCDR3!=""),] #
data_full_cdr3$CDR3_length<-nchar(as.character(data_full_cdr3$nSeqCDR3)) 
data_full_cdr3<-data_full_cdr3[which(data_full_cdr3$CDR3_length!=3),] #3,591,002

##Obtain the ID for the clone call
data_full_cdr3$seqID<-seq(1,nrow(data_full_cdr3))
data_full_cdr3$SEQUENCE_ID<-paste(data_full_cdr3$Sample,data_full_cdr3$seqID,sep="_")

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
data_clonesInference_Ig<-data_full_cdr3_Ig[,c("SEQUENCE_ID","Sample","nSeqCDR3","CDR3_length","bestVGene","bestJGene","V_J_lenghCDR3")]
data_clonesInference_TCR<-data_full_cdr3_TCR[,c("SEQUENCE_ID","Sample","nSeqCDR3","CDR3_length","bestVGene","bestJGene","V_J_lenghCDR3")]
write.table(data_clonesInference_Ig,file="Data/GTEX/MIXCR/data_for_cloneInfered_Ig_GTEX.txt",row.names = F,sep="\t")
write.table(data_clonesInference_TCR,file="Data/GTEX/MIXCR/data_for_cloneInfered_TCR_GTEX.txt",row.names = F,sep="\t")

### After passing the nucleotides.py
##Read the clones and merge with the data
nucleotides_Ig<-read.csv("Data/GTEx/MIXCR/ClonesInfered_Ig_GTEX.csv")
nucleotides_TCR<-read.csv("Data/GTEx/MIXCR/ClonesInfered_TCR_GTEX.csv")
nucleotides<-rbind(nucleotides_Ig,nucleotides_TCR)
data_merge<-merge(data_full_cdr3,nucleotides[,c("SEQUENCE_ID","CloneId")],by=c("SEQUENCE_ID"))

###Add a column with the reads per chain for Ig and TRA
##Clones per chain
data_merge$V_J_lenghCDR3_CloneId = paste(data_merge$V_J_lenghCDR3,data_merge$CloneId,sep="_")

clones_count<- unique(data_merge[,c("Sample","V_J_lenghCDR3_CloneId","chainType")])
clones<-data.matrix(table(clones_count$Sample,clones_count$chainType))
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
  data_sample_unique<-data_merge[which(data_merge$Sample==sample[i]),]
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
  #write.delim(data.frame(table(table(clones_sample_TRA))),file=paste("clones_sample_TRA_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  #write.delim(data.frame(table(table(clones_sample_TRB))),file=paste("clones_sample_TRB_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  #write.delim(data.frame(table(table(clones_sample_TRD))),file=paste("clones_sample_TRD_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  #write.delim(data.frame(table(table(clones_sample_TRG))),file=paste("clones_sample_TRG_",sample[i],".txt",sep=""),sep="\t",col.names=F)
  
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
#recon<-read.table("Data/PAAD/RECON/test_D_number_table.txt",header=T)
#chain<-substr(recon$sample_name,61,63)
#sample<-substr(recon$sample_name,65,100)
##0.0D is species richness (Number of clones)
##Entropy is ln(1.0.D)
# diversity<-as.data.frame(diversity)
# chain_list<-unique(chain)
# for(i in chain_list){
#   recon_chain<-recon[which(chain==i),]
#   sample_chain<-substr(recon_chain$sample_name,65,100)
#   id<-match(rownames(diversity),sample_chain)
#   clone_chain<-ifelse(is.na(id)==F,recon_chain[id,"est_0.0D"],0)
#   assign(paste0("clones_recon_",i),clone_chain)
#   entroy_chain<-ifelse(is.na(id)==F,ln(recon_chain[id,"est_1.0D"]),0)
#   assign(paste0("entropy_recon_",i),entroy_chain)
# }
# diversity$clones_recon_IGH<-clones_recon_IGH
# diversity$clones_recon_IGK<-clones_recon_IGK
# diversity$clones_recon_IGL<-clones_recon_IGL
# diversity$clones_recon_TRA<-clones_recon_TRA
# diversity$clones_recon_TRB<-clones_recon_TRB
# diversity$clones_recon_TRD<-clones_recon_TRD
# diversity$clones_recon_TRG<-clones_recon_TRG
# 
# diversity$entropy_recon_IGH<-entropy_recon_IGH
# diversity$entropy_recon_IGK<-entropy_recon_IGK
# diversity$entropy_recon_IGL<-entropy_recon_IGL
# diversity$entropy_recon_TRA<-entropy_recon_TRA
# diversity$entropy_recon_TRB<-entropy_recon_TRB
# diversity$entropy_recon_TRD<-entropy_recon_TRD
# diversity$entropy_recon_TRG<-entropy_recon_TRG
#Simpson Index (1/2D) and BPI (1/∞D)

load("Data/GTEX/MIXCR/GTEX_RepertoireResults_diversity.Rdata")
id<-match(rownames(summaryMatrix_GTEX),rownames(diversity))
GTEX.repertoire.diversity<-cbind(summaryMatrix_GTEX,diversity[id,])
data_merge_GTEX<-data_merge
save(data_merge_GTEX,GTEX.repertoire.diversity,file="Data/GTEX/MIXCR/GTEX_FullData.Rdata")




