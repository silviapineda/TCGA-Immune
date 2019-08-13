library(pgirmess)
library(stringr)
library(SciViews)

load("Data/GTEx/Blood/MIXCR/GTEX_FullData.Rdata")
totalReads<-read.table("Data/GTEx/Blood/MIXCR/total_reads_GTEX.txt",sep=";")
id<-match(rownames(GTEX.repertoire.diversity),totalReads$V1)
GTEX.repertoire.diversity$TotalSeq<-totalReads[id,"V2"]
rownames(GTEX.repertoire.diversity)<-unlist(strsplit(as.character(rownames(GTEX.repertoire.diversity)), "\\}"))

####Annotation 
sra<-read.csv("Data/GTEx/Blood/SraRunTableBloodRNAseq.csv")
##Read gene reads to obtain the samples that have passed QC in GTEX
#Gene_reads<-read.table("Data/GTEx/TEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct",sep="\t",header=T)
#samples<-gsub(".", '-', as.character(colnames(Gene_reads)), fixed = T)
#samples<-samples[-c(1:2)]
#write.csv(samples,"Data/GTEx/samples_QC_GTEX.csv")
samples<-read.csv("Data/GTEx/samples_QC_GTEX.csv")
id<-match(samples$x, sra$Sample_Name)
sra<-sra[na.omit(id),]

##Match with GTEX blood diversity
id<-match(rownames(GTEX.repertoire.diversity),sra$Run)
GTEX.repertoire.diversity<-GTEX.repertoire.diversity[which(is.na(id)==F),]
GTEX.repertoire.diversity$SUBJID<-sra$submitted_subject_id[na.omit(id)]

annotation_gtex<-read.csv("Data/GTEx/GTEX_annotation_phenotypes.csv")
id<-match(GTEX.repertoire.diversity$SUBJID,annotation_gtex$SUBJID)
annotation_gtex_blood<-annotation_gtex[id,]

id<-match(data_merge_GTEX$Sample,rownames(GTEX.repertoire.diversity))
data_merge_GTEX<-data_merge_GTEX[which(is.na(id)!=F),]

##To obtain RECON files
sample<-unique(data_merge_GTEX$Sample)
for (i in 1:length(sample)){
  print(i)
  data_sample_unique<-data_merge_GTEX[which(data_merge_GTEX$Sample==sample[i]),]
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
}
