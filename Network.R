rm(list = ls(all = TRUE))
x <-date()
print(x)
###########################################################################################
### PROJECT: Immune Repertoire. Analysis B cells antibodies for pregnancy outcomes
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Network analysis
###         
###
### Author: Silvia Pineda
### Date: January, 2019
############################################################################################
library(ggplot2)
library("igraph")
library("ineq")
library("dplyr")
library("RColorBrewer")

setwd("~/TCGA-Immune/")
load("Data/PAAD/PAAD_FullData.Rdata")

##1.Obtain the vertex and edges
#Receptor
Obtain_vertex_edges<-function(data,receptor){
  data<-data[which(data$receptor==receptor),]
  data$CloneId_CDR3<-paste0(data[,c("V_J_lenghCDR3_CloneId")],data[,c("nSeqCDR3")])
  sample<-unique(data$sample)
  for (i in sample){
    print(i)
    data_sample<-data[which(data$sample==i),]
    df_sample<-data_sample[,c("CloneId_CDR3","V_J_lenghCDR3_CloneId")]
    groups <- group_by(df_sample,CloneId_CDR3,V_J_lenghCDR3_CloneId)
    assign(paste0("edges",i),unique(data.frame(groups)))
    df_vertex<-data.frame(table(data_sample$CloneId_CDR3))
    assign(paste0("vertex",i),df_vertex[which(df_vertex$Freq!=0),])
    write.table(get(paste0("edges",i)),paste0("Results/Network/edges_",receptor,"_",i,".txt"),sep="\t",row.names = F)
    write.table(get(paste0("vertex",i)),paste0("Results/Network/vertex_",receptor,"_",i,".txt"),sep="\t",row.names = F)
  }
}

#chainType
Obtain_vertex_edges<-function(data,chainType){
  data<-data[which(data$chainType==chainType),]
  data$CloneId_CDR3<-paste0(data[,c("V_J_lenghCDR3_CloneId")],data[,c("nSeqCDR3")])
  sample<-unique(data$sample)
  for (i in sample){
    print(i)
    data_sample<-data[which(data$sample==i),]
    df_sample<-data_sample[,c("CloneId_CDR3","V_J_lenghCDR3_CloneId")]
    groups <- group_by(df_sample,CloneId_CDR3,V_J_lenghCDR3_CloneId)
    assign(paste0("edges",i),unique(data.frame(groups)))
    df_vertex<-data.frame(table(data_sample$CloneId_CDR3))
    assign(paste0("vertex",i),df_vertex[which(df_vertex$Freq!=0),])
    write.table(get(paste0("edges",i)),paste0("Results/Network/edges_",chainType,"_",i,".txt"),sep="\t",row.names = F)
    write.table(get(paste0("vertex",i)),paste0("Results/Network/vertex_",chainType,"_",i,".txt"),sep="\t",row.names = F)
  }
}

network_IGHV<-Obtain_vertex_edges(data_merge,"IGHV")
network_IGKV<-Obtain_vertex_edges(data_merge,"IGKV")
network_IGLV<-Obtain_vertex_edges(data_merge,"IGLV")
#network_TRAV<-Obtain_vertex_edges(data_merge,"TRAV")
#network_TRBV<-Obtain_vertex_edges(data_merge,"TRBV")
#network_TRDV<-Obtain_vertex_edges(data_merge,"TRDV")
#network_TRGV<-Obtain_vertex_edges(data_merge,"TRGV")

data_merge$receptor<-ifelse(data_merge$chainType=="IGHV" | data_merge$chainType=="IGKV" | data_merge$chainType=="IGLV","IG","TCR")
network_IG<-Obtain_vertex_edges(data_merge,"IG")


##2.Apply the nucleotides-assembly-1.0.jar made by Mikel using the Network.sh 

##3.Plot the network
sample<-unique(data_merge$sample) 
#grep("cab8e91a-ca41-4c62-af53-3aa4057d68d5",sample) #IGH
#sample<-sample[-grep("cab8e91a-ca41-4c62-af53-3aa4057d68d5",sample)] #IGH

id<-grep("c174b41a-84c4-4a33-9c21-48dee5029ddb",sample) #IGKV
sample<-sample[-id]


receptor="IGLV"
for(i in sample) {
  print(i)
  edges <- read.delim(paste("Results/Network/",receptor,"/edges_",receptor,"_",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("Results/Network/",receptor,"/vertex_",receptor,"_",i,".txt",sep = ""))
  if(length(edges$edge1)!=0){
    net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
    V(net)$size <- V(net)$Freq/100
    V(net)$color <- c("#c9a47f")
    net <- simplify(net, remove.multiple = F, remove.loops = T) 
    E(net)$arrow.mode <- 0
    E(net)$width <- 0.4
    E(net)$color <- c("black")
    tiff(paste("Results/Network/",receptor,"/network_",receptor,"_",i,".tiff",sep=""),res=300,h=3000,w=3000)
    plot(net,vertex.label=NA,layout=layout_with_graphopt(net,niter=800))
    dev.off()
  }
}


##4. Obtain the vertex and cluster Gini index
#Receptor
Obtain_gini_index<-function(data,receptor,PAAD.repertoire.diversity){
  sample<-rownames(PAAD.repertoire.diversity)
  #grep("cab8e91a-ca41-4c62-af53-3aa4057d68d5",sample)
  #sample<-sample[-grep("cab8e91a-ca41-4c62-af53-3aa4057d68d5",sample)]
  
  vertex_max<-NULL
  vertex_gini<-NULL
  cluster_max<-NULL
  cluster_gini<-NULL
  num_reads_max_cluster<-NULL
  clusters<-NULL
  j<-1
  
  for (i in sample){
    assign(paste0("edges",i),read.delim(paste0("Results/Network/",receptor,"/edges_",receptor,"_",i,".txt")))
    assign(paste0("vertex",i),read.delim(paste0("Results/Network/",receptor,"/vertex_",receptor,"_",i,".txt")))
    vertex_max[j]<-max(get(paste0("vertex",i))$Freq)
    vertex_gini[j]<-Gini(get(paste0("vertex",i))$Freq)
    cluster_max[j]<-max(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId))
    clusters[j]<-sum(table(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId)))
    num_reads_max_cluster[j]<-tail(table(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId)),1)
    cluster_gini[j]<-Gini(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId))
    j=j+1
  }
  
  #clonal_expansion<-(num_reads_max_cluster/summary_data[,isotype])*100
  results<-cbind(cluster_gini,vertex_gini,vertex_max,cluster_max,num_reads_max_cluster,clusters)
  
  rownames(results)<-sample
  
  return(results)
}

#chainType
Obtain_gini_index<-function(data,chainType,PAAD.repertoire.diversity){
  sample<-rownames(PAAD.repertoire.diversity)
  id<-grep("c174b41a-84c4-4a33-9c21-48dee5029ddb",sample) #IGKV
  sample<-sample[-id]
  
  vertex_max<-NULL
  vertex_gini<-NULL
  cluster_max<-NULL
  cluster_gini<-NULL
  num_reads_max_cluster<-NULL
  clusters<-NULL
  j<-1
  
  for (i in sample){
    print(i)
    assign(paste0("edges",i),read.delim(paste0("Results/Network/",chainType,"/edges_",chainType,"_",i,".txt")))
    assign(paste0("vertex",i),read.delim(paste0("Results/Network/",chainType,"/vertex_",chainType,"_",i,".txt")))
    vertex_max[j]<-max(get(paste0("vertex",i))$Freq)
    vertex_gini[j]<-Gini(get(paste0("vertex",i))$Freq)
    cluster_max[j]<-max(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId))
    clusters[j]<-sum(table(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId)))
    num_reads_max_cluster[j]<-tail(table(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId)),1)
    cluster_gini[j]<-Gini(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId))
    j=j+1
  }
  
  #clonal_expansion<-(num_reads_max_cluster/summary_data[,isotype])*100
  results<-cbind(cluster_gini,vertex_gini,vertex_max,cluster_max,num_reads_max_cluster,clusters)
  
  rownames(results)<-sample
  
  return(results)
}

###Plot Gini boxplot
receptor = "IG"
assign(paste0("cluster_gini_",receptor),data.frame(Obtain_gini_index(data_merge,receptor,PAAD.repertoire.diversity)))
assign(paste0("PAAD.repertoire.diversity.gini_",receptor),merge(PAAD.repertoire.diversity,cluster_gini_IG,by="row.names"))
write.csv(get(paste0("PAAD.repertoire.diversity.gini_",receptor)),file=paste0("Results/Network/network",receptor,".csv"))
PAAD.repertoire.diversity.gini_IG_tumor<-PAAD.repertoire.diversity.gini_IG[which(PAAD.repertoire.diversity.gini_IG$Tumor_type=="Tumor_pancreas"),]
write.csv(get(paste0("PAAD.repertoire.diversity.gini_",receptor,"_tumor")),file=paste0("Results/Network/network_",receptor,"_tumor.csv"))

tiff(paste0("Results/Network/network_vertex_cluster_gini_",receptor,".tiff"),h=2000,w=2000,res=300)
plot(get(paste0("PAAD.repertoire.diversity.gini_",receptor,"_tumor"))[,"cluster_gini"], 
     get(paste0("PAAD.repertoire.diversity.gini_",receptor,"_tumor"))[,"vertex_gini"],pch=20,ylab = "Gini (Vextex)",xlab = "Gini (Cluster)")
dev.off()


chainType= "IGLV"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity)))
assign(paste0("PAAD.repertoire.diversity.gini_",chainType),merge(PAAD.repertoire.diversity,get(paste0("cluster_gini_",chainType)),by="row.names"))
#write.csv(get(paste0("PAAD.repertoire.diversity.gini_",chainType)),file=paste0("Results/Network/network",chainType,".csv"))
assign(paste0("PAAD.repertoire.diversity.gini_",chainType,"_tumor"),
       get(paste0("PAAD.repertoire.diversity.gini_",chainType))[which(get(paste0("PAAD.repertoire.diversity.gini_",chainType))$Tumor_type=="Tumor_pancreas"),])
write.csv(get(paste0("PAAD.repertoire.diversity.gini_",chainType,"_tumor")),file=paste0("Results/Network/network_",chainType,"_tumor.csv"))
tiff(paste0("Results/Network/network_vertex_cluster_gini_",chainType,".tiff"),h=2000,w=2000,res=300)
#par(fig=c(0,0.8,0,0.8))
plot(get(paste0("PAAD.repertoire.diversity.gini_",chainType,"_tumor"))[,"cluster_gini"], 
     get(paste0("PAAD.repertoire.diversity.gini_",chainType,"_tumor"))[,"vertex_gini"],pch=20,ylab = "Gini (Vextex)",xlab = "Gini (Cluster)")
dev.off()


###Info for the three samples with very high gini(vertex) and gini(cluster)
samples_high<-c("5e9e81e2-0e8b-4aca-aced-6ce451fa3262", "08c4c14a-94a8-4063-bbc3-916579960078", "49170bdc-8725-43a0-b2d7-ccb52fa415e6")
samples_high_TCGA<-substr(as.character(PAAD.repertoire.diversity.gini_IG[match(samples_high,PAAD.repertoire.diversity.gini_IG$Row.names),"TCGA_sample"]),1,12)

clinical.patient[match(samples_high_TCGA,clinical.patient$bcr_patient_barcode),]
