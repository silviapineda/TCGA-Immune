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
load("Data/Pancreas_Validation/Pancreas_Validation_FullData.Rdata")
##################################
##1.Obtain the vertex and edges
##################################
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
    write.table(get(paste0("edges",i)),paste0("Results/Validation/Tumor//Network/edges_",chainType,"_",i,".txt"),sep="\t",row.names = F)
    write.table(get(paste0("vertex",i)),paste0("Results/Validation/Tumor/Network/vertex_",chainType,"_",i,".txt"),sep="\t",row.names = F)
  }
}

network_IGH<-Obtain_vertex_edges(data_merge,"IGH")
network_IGK<-Obtain_vertex_edges(data_merge,"IGK")
network_IGL<-Obtain_vertex_edges(data_merge,"IGL")
network_TRA<-Obtain_vertex_edges(data_merge,"TRA")
network_TRB<-Obtain_vertex_edges(data_merge,"TRB")
network_TRD<-Obtain_vertex_edges(data_merge,"TRD")
network_TRG<-Obtain_vertex_edges(data_merge,"TRG")

#data_merge$receptor<-ifelse(data_merge$chainType=="IGH" | data_merge$chainType=="IGK" | data_merge$chainType=="IGL","IG","TCR")
#network_IG<-Obtain_vertex_edges(data_merge,"IG")

########
##2.Apply the nucleotides-assembly-1.0.jar made by Mikel using the Network.sh 
########

###############################################
##3. Obtain the vertex and cluster Gini index
###############################################

#chainType
Obtain_gini_index<-function(data,chainType,Pancreas.Validation.repertoire.diversity){
  sample<-rownames(Pancreas.Validation.repertoire.diversity)
  
  vertex_max<-NULL
  vertex_gini<-NULL
  cluster_max<-NULL
  cluster_gini<-NULL
  num_reads_max_cluster<-NULL
  clusters<-NULL
  j<-1
  sample2<-NULL
  
  for (i in sample){
    print(i)
    res<-try(assign(paste0("edges",i),read.delim(paste0("Results/Validation/Tumor/Network/edges_",chainType,"_",i,".txt"))))
    if(class(res) != "try-error"){ 
    assign(paste0("edges",i),read.delim(paste0("Results/Validation/Tumor/Network/edges_",chainType,"_",i,".txt")))
    assign(paste0("vertex",i),read.delim(paste0("Results/Validation/Tumor/Network/vertex_",chainType,"_",i,".txt")))
    vertex_max[j]<-max(get(paste0("vertex",i))$Freq)
    vertex_gini[j]<-Gini(get(paste0("vertex",i))$Freq)
    cluster_max[j]<-max(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId))
    clusters[j]<-sum(table(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId)))
    num_reads_max_cluster[j]<-tail(table(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId)),1)
    cluster_gini[j]<-Gini(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId))
    j=j+1
    sample2<-c(sample2,i)
    }
  }
  
  #clonal_expansion<-(num_reads_max_cluster/summary_data[,chainType])*100
  results<-cbind(cluster_gini,vertex_gini,vertex_max,cluster_max,num_reads_max_cluster,clusters)
  
  rownames(results)<-sample2
  
  return(results)
}


chainType= "IGH"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,Pancreas.Validation.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(Pancreas.Validation.repertoire.diversity))
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]
chainType= "IGK"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,Pancreas.Validation.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(Pancreas.Validation.repertoire.diversity))
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]
chainType= "IGL"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,Pancreas.Validation.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(Pancreas.Validation.repertoire.diversity))
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]
chainType= "TRA"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,Pancreas.Validation.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(Pancreas.Validation.repertoire.diversity))
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]
chainType= "TRB"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,Pancreas.Validation.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(Pancreas.Validation.repertoire.diversity))
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]
chainType= "TRD"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,Pancreas.Validation.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(Pancreas.Validation.repertoire.diversity))
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]
chainType= "TRG"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,Pancreas.Validation.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(Pancreas.Validation.repertoire.diversity))
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
Pancreas.Validation.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
Pancreas.Validation.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]

##To save the varibles with vertex and cluster
save(data_merge,Pancreas.Validation.repertoire.diversity,file="Data/Pancreas_Validation/Pancreas_Validation_FullData.Rdata")


########################
##4.Plot the network
########################
#Tumor

chainType="TRA"
Pancreas.Validation.repertoire.diversity.TRA<-Pancreas.Validation.repertoire.diversity[which(is.na(Pancreas.Validation.repertoire.diversity$cluster_gini_TRA)==F),]
sample_tumor<-rownames(Pancreas.Validation.repertoire.diversity.TRA)
for(i in sample_tumor) {
  print(i)
  edges <- read.delim(paste("Results/Validation/Tumor/Network/edges_",chainType,"_",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("Results/Validation/Tumor/Network/vertex_",chainType,"_",i,".txt",sep = ""))
  if(length(edges$edge1)!=0){
    net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
    V(net)$size <- V(net)$Freq
    V(net)$color <- c("#B3CDE3")
    net <- simplify(net, remove.multiple = F, remove.loops = T) 
    E(net)$arrow.mode <- 0
    E(net)$width <- 0.4
    E(net)$color <- c("black")
    tiff(paste("Results/Validation/Tumor//Network/network_",chainType,"_",i,".tiff",sep=""),res=300,h=3000,w=3000)
    plot(net,vertex.label=NA,layout=layout_with_graphopt(net,niter=800,charge=0.01))
    dev.off()
  }
}

#Normal
# sample_normal<-Pancreas.Validation.repertoire.diversity$sample[which(Pancreas.Validation.repertoire.diversity$tissue=="normal pancreas")]
# chainType="IGH"
# for(i in sample_normal) {
#   print(i)
#   edges <- read.delim(paste("Results/Validation/Tumor//Network/edges_",chainType,"_",i,".txt.outcome.txt",sep = ""))
#   vertex <- read.delim(paste("Results/Validation/Tumor//Network/vertex_",chainType,"_",i,".txt",sep = ""))
#   if(length(edges$edge1)!=0){
#     net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
#     V(net)$size <- V(net)$Freq/100
#     V(net)$color <- c("#F4CAE4")
#     net <- simplify(net, remove.multiple = F, remove.loops = T) 
#     E(net)$arrow.mode <- 0
#     E(net)$width <- 0.4
#     E(net)$color <- c("black")
#     tiff(paste("Results/Validation/Network/network_",chainType,"_",i,".tiff",sep=""),res=300,h=3000,w=3000)
#     plot(net,vertex.label=NA,layout=layout_with_graphopt(net,niter=800,charge=0.01))
#     dev.off()
#   }
# }
