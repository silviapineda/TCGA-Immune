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
load("Data/GTEx/Pancreas/GTEx_FullData.Rdata")
#load("Data/GTEx/Blood/GTEx_blood_FullData.Rdata")
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
    write.table(get(paste0("edges",i)),paste0("Results/GTEx//Network/edges_",chainType,"_",i,".txt"),sep="\t",row.names = F)
    write.table(get(paste0("vertex",i)),paste0("Results/GTEx//Network/vertex_",chainType,"_",i,".txt"),sep="\t",row.names = F)
  }
}

network_IGH<-Obtain_vertex_edges(data_merge_pancreas,"IGH")
network_IGK<-Obtain_vertex_edges(data_merge,"IGK")
network_IGL<-Obtain_vertex_edges(data_merge,"IGL")
#network_TRAV<-Obtain_vertex_edges(data_merge,"TRAV")
#network_TRBV<-Obtain_vertex_edges(data_merge,"TRBV")
#network_TRDV<-Obtain_vertex_edges(data_merge,"TRDV")
#network_TRGV<-Obtain_vertex_edges(data_merge,"TRGV")

#data_merge$receptor<-ifelse(data_merge$chainType=="IGH" | data_merge$chainType=="IGK" | data_merge$chainType=="IGL","IG","TCR")
#network_IG<-Obtain_vertex_edges(data_merge,"IG")

########
##2.Apply the nucleotides-assembly-1.0.jar made by Mikel using the Network.sh 
########

###############################################
##3. Obtain the vertex and cluster Gini index
###############################################

#chainType
Obtain_gini_index<-function(data,chainType,PAAD.GTEx.repertoire.diversity){
  sample<-rownames(PAAD.GTEx.repertoire.diversity)
  #id<-grep("cab8e91a-ca41-4c62-af53-3aa4057d68d5",sample) #IGH
  #sample<-sample[-id]
  id<-grep("SRR1089537",sample) #IGH
  sample<-sample[-id]
  id<-grep("SRR1095479",sample) #IGK
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
    assign(paste0("edges",i),read.delim(paste0("Results/GTEx/Network/edges_",chainType,"_",i,".txt")))
    assign(paste0("vertex",i),read.delim(paste0("Results/GTEx/Network/vertex_",chainType,"_",i,".txt")))
    vertex_max[j]<-max(get(paste0("vertex",i))$Freq)
    vertex_gini[j]<-Gini(get(paste0("vertex",i))$Freq)
    cluster_max[j]<-max(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId))
    clusters[j]<-sum(table(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId)))
    num_reads_max_cluster[j]<-tail(table(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId)),1)
    cluster_gini[j]<-Gini(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId))
    j=j+1
  }
  
  #clonal_expansion<-(num_reads_max_cluster/summary_data[,chainType])*100
  results<-cbind(cluster_gini,vertex_gini,vertex_max,cluster_max,num_reads_max_cluster,clusters)
  
  rownames(results)<-sample
  
  return(results)
}


chainType= "IGH"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge_pancreas,chainType,Pancreas.repertoire.diversity)))
id<-match(rownames(cluster_gini_IGH),rownames(Pancreas.repertoire.diversity))
Pancreas.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
Pancreas.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
Pancreas.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
Pancreas.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]

##To save the varibles with vertex and cluster
save(data_merge_pancreas,annotation_gtex_pancreas,Pancreas.repertoire.diversity,file="Data/GTEx/GTEx_FullData.Rdata")

Pancreas.repertoire.diversity.filter<-
  Pancreas.repertoire.diversity[which(Pancreas.repertoire.diversity[,paste0("cluster_gini_",chainType)]!=0 &
                                        Pancreas.repertoire.diversity[,paste0("vertex_gini_",chainType)]!=0 ),]

tiff(paste0("Results/Network/network_vertex_cluster_gini_",chainType,".tiff"),h=2000,w=2000,res=300)
brewer.pal(n = 5, name = "Accent")
cols=c( "#7FC97F","#386CB0","#FDC086", "#BEAED4")

par(fig=c(0,0.8,0,0.8))
plot(PAAD.GTEx.repertoire.diversity.filter[,paste0("cluster_gini_",chainType)], 
     PAAD.GTEx.repertoire.diversity.filter[,paste0("vertex_gini_",chainType)],
     col = cols[factor(PAAD.GTEx.repertoire.diversity.filter$outcome)],
     pch=20,ylab = c("Gini (Vextex)"),xlab = c("Gini (Cluster)"))
    legend("bottomright",legend=c("normal-pancreas (GTEx)", "normal-pancreas (TCGA)",
                                  "pseudonormal-pancreas (TCGA)","tumor-pancreas (TCGA)"), 
           col=cols,pch=20,cex=0.8)

par(fig=c(0,0.8,0.55,1), new=TRUE)
summary(glm(PAAD.GTEx.repertoire.diversity.filter[,paste0("cluster_gini_",chainType)]~PAAD.GTEx.repertoire.diversity.filter$outcome))
boxplot(PAAD.GTEx.repertoire.diversity.filter[,paste0("cluster_gini_",chainType)]~PAAD.GTEx.repertoire.diversity.filter$outcome,
        col=cols, horizontal=TRUE, axes=FALSE)

par(fig=c(0.65,1,0,0.8),new=TRUE)
summary(glm(PAAD.GTEx.repertoire.diversity.filter[,paste0("vertex_gini_",chainType)]~PAAD.GTEx.repertoire.diversity.filter$outcome))
boxplot(PAAD.GTEx.repertoire.diversity.filter[,paste0("vertex_gini_",chainType)]~PAAD.GTEx.repertoire.diversity.filter$outcome,
        col=cols, axes=FALSE)
dev.off()

########################
##4.Plot the network
########################

#Normal
sample_normal<-rownames(Pancreas.repertoire.diversity.filter)
chainType="IGH"
for(i in sample_normal) {
  print(i)
  edges <- read.delim(paste("Results/GTEx/Network/edges_",chainType,"_",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("Results/GTEx/Network/vertex_",chainType,"_",i,".txt",sep = ""))
  if(length(edges$edge1)!=0){
    net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
    V(net)$size <- V(net)$Freq/200
    V(net)$color <- c("#386CB0")
    net <- simplify(net, remove.multiple = F, remove.loops = T) 
    E(net)$arrow.mode <- 0
    E(net)$width <- 0.4
    E(net)$color <- c("black")
    tiff(paste("Results/GTEx/Network/network_",chainType,"_",i,".tiff",sep=""),res=300,h=3000,w=3000)
    plot(net,vertex.label=NA,layout=layout_with_graphopt(net,niter=800,charge=0.01))
    dev.off()
  }
}

