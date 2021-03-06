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
    write.table(get(paste0("edges",i)),paste0("Data/PAAD/Network/edges_",chainType,"_",i,".txt"),sep="\t",row.names = F)
    write.table(get(paste0("vertex",i)),paste0("Data/PAAD/Network/vertex_",chainType,"_",i,".txt"),sep="\t",row.names = F)
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
#Receptor
Obtain_gini_index<-function(data,receptor,PAAD.repertoire.diversity){
  sample<-rownames(PAAD.repertoire.diversity)
  grep("cab8e91a-ca41-4c62-af53-3aa4057d68d5",sample)
  sample<-sample[-grep("cab8e91a-ca41-4c62-af53-3aa4057d68d5",sample)]
  
  vertex_max<-NULL
  vertex_gini<-NULL
  cluster_max<-NULL
  cluster_gini<-NULL
  num_reads_max_cluster<-NULL
  clusters<-NULL
  j<-1
  
  for (i in sample){
    assign(paste0("edges",i),read.delim(paste0("Data/PAAD/edges_",receptor,"_",i,".txt")))
    assign(paste0("vertex",i),read.delim(paste0("Data/PAAD/vertex_",receptor,"_",i,".txt")))
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
  vertex_max<-NULL
  vertex_gini<-NULL
  cluster_max<-NULL
  cluster_gini<-NULL
  num_reads_max_cluster<-NULL
  clusters<-NULL
  clonal_expansion<-NULL
  j<-1
  CDR3_length<-NULL
  vertex_size<-NULL
  sample2<-NULL
  for (i in sample){
    print(i)
    res<-try(assign(paste0("edges",i),read.delim(paste0("Data/PAAD/Network/edges_",chainType,"_",i,".txt"))))
    if(class(res) != "try-error"){ 
      assign(paste0("edges",i),read.delim(paste0("Data/PAAD/Network/edges_",chainType,"_",i,".txt")))
      assign(paste0("vertex",i),read.delim(paste0("Data/PAAD/Network/vertex_",chainType,"_",i,".txt")))
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
  
  
  results<-cbind(cluster_gini,vertex_gini,vertex_max,cluster_max,num_reads_max_cluster,clusters)
  rownames(results)<-sample2
  
  return(results)
}


chainType= "IGH"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(PAAD.repertoire.diversity))
PAAD.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
PAAD.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]
chainType= "IGK"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(PAAD.repertoire.diversity))
PAAD.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
PAAD.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]
chainType= "IGL"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(PAAD.repertoire.diversity))
PAAD.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
PAAD.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]
chainType= "TRA"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(PAAD.repertoire.diversity))
PAAD.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
PAAD.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]
chainType= "TRB"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(PAAD.repertoire.diversity))
PAAD.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
PAAD.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]
chainType= "TRD"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(PAAD.repertoire.diversity))
PAAD.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
PAAD.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]
chainType= "TRG"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity)))
id<-match(rownames(get(paste0("cluster_gini_",chainType))),rownames(PAAD.repertoire.diversity))
PAAD.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
PAAD.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
PAAD.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]



##To save the varibles with vertex and cluster
save(data_merge,PAAD.repertoire.diversity,xCell.data.PAAD,xCell.pvalue.PAAD,clinical.drug,clinical.patient,clinical.radiation,clinical.new_tumor_event,clinical.folow_up,biospecimen.slide,annotation
     ,file="Data/PAAD/PAAD_FullData.Rdata")

####
load("Data/PAAD/PAAD_FullData.Rdata")
chainType= "IGH"
PAAD.repertoire.diversity.filter<-
  PAAD.repertoire.diversity[which(PAAD.repertoire.diversity[,paste0("cluster_gini_",chainType)]!=0 &
                                         PAAD.repertoire.diversity[,paste0("vertex_gini_",chainType)]!=0 ),]

tiff(paste0("Results/PAAD/Network/network_vertex_cluster_gini_",chainType,".tiff"),h=2000,w=2000,res=300)
brewer.pal(3,name = "Accent")
brewer.pal(3,name = "Pastel1")
cols=c( "#7FC97F", "#FBB4AE","#BEAED4", "#FDC086")

par(fig=c(0,0.8,0,0.8))
plot(PAAD.repertoire.diversity.filter[,paste0("cluster_gini_",chainType)], 
     PAAD.repertoire.diversity.filter[,paste0("vertex_gini_",chainType)],
     col = cols[factor(PAAD.repertoire.diversity.filter$Tumor_type_4categ)],
     pch=20,ylab = c("Gini (Vextex)"),xlab = c("Gini (Cluster)"))
    legend("bottomright",legend=levels(PAAD.repertoire.diversity.filter$Tumor_type_4categ), 
           col=cols,pch=20,cex=0.8)

par(fig=c(0,0.8,0.55,1), new=TRUE)
summary(glm(PAAD.repertoire.diversity.filter[,paste0("cluster_gini_",chainType)]~PAAD.repertoire.diversity.filter$Tumor_type_4categ))
boxplot(PAAD.repertoire.diversity.filter[,paste0("cluster_gini_",chainType)]~PAAD.repertoire.diversity.filter$Tumor_type_4categ,
        col=cols, horizontal=TRUE, axes=FALSE)

par(fig=c(0.65,1,0,0.8),new=TRUE)
summary(glm(PAAD.repertoire.diversity.filter[,paste0("vertex_gini_",chainType)]~PAAD.repertoire.diversity.filter$Tumor_type_4categ))
boxplot(PAAD.repertoire.diversity.filter[,paste0("vertex_gini_",chainType)]~PAAD.repertoire.diversity.filter$Tumor_type_4categ,
        col=cols, axes=FALSE)
dev.off()

########################
##4.Plot the network
########################

#Tumor
chainType="IGK"
PAAD.repertoire.diversity.IGK<-PAAD.repertoire.diversity[which(is.na(PAAD.repertoire.diversity$cluster_gini_IGK)==F),]
sample_tumor<-rownames(PAAD.repertoire.diversity.IGK)[which(PAAD.repertoire.diversity.IGK$Tumor_type_4categ=="PDAC")]


for(i in sample_tumor) {
  print(i)
  edges <- read.delim(paste("Data/PAAD/Network/edges_",chainType,"_",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("Data/PAAD/Network/vertex_",chainType,"_",i,".txt",sep = ""))
  if(length(edges$edge1)!=0){
    net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
    V(net)$size <- V(net)$Freq/100
    V(net)$color <- c("#BEAED4")
    net <- simplify(net, remove.multiple = F, remove.loops = T) 
    E(net)$arrow.mode <- 0
    E(net)$width <- 0.4
    E(net)$color <- c("black")
    tiff(paste("Results/PAAD/Network/network_",chainType,"_",i,".tiff",sep=""),res=300,h=3000,w=3000)
    plot(net,vertex.label=NA,layout=layout_with_graphopt(net,niter=800,charge=0.01,mass=30))
    dev.off()
  }
}

#Normal
sample_normal<-rownames(PAAD.repertoire.diversity)[which(PAAD.repertoire.diversity$Tumor_type_3categ=="normal_pancreas")]
for(i in sample_normal) {
  print(i)
  edges <- read.delim(paste("Data/PAAD/Network/edges_",chainType,"_",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("Data/PAAD/Network/vertex_",chainType,"_",i,".txt",sep = ""))
  if(length(edges$edge1)!=0){
    net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
    V(net)$size <- V(net)$Freq/100
    V(net)$color <- c("#FDC086")
    net <- simplify(net, remove.multiple = F, remove.loops = T) 
    E(net)$arrow.mode <- 0
    E(net)$width <- 0.4
    E(net)$color <- c("black")
    tiff(paste("Results/PAAD/Network/network_",chainType,"_normal_",i,".tiff",sep=""),res=300,h=3000,w=3000)
    plot(net,vertex.label=NA,layout=layout_with_graphopt(net,niter=800,charge=0.01))
    dev.off()
  }
}

#PseudoNormal
sample_normal<-rownames(PAAD.repertoire.diversity)[which(PAAD.repertoire.diversity$Tumor_type_3categ=="pseudonormal_pancreas")]
for(i in sample_normal) {
  print(i)
  edges <- read.delim(paste("Data/PAAD/Network/edges_",chainType,"_",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("Data/PAAD/Network/vertex_",chainType,"_",i,".txt",sep = ""))
  if(length(edges$edge1)!=0){
    net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
    V(net)$size <- V(net)$Freq/100
    V(net)$color <- c("#B3CDE3")
    net <- simplify(net, remove.multiple = F, remove.loops = T) 
    E(net)$arrow.mode <- 0
    E(net)$width <- 0.4
    E(net)$color <- c("black")
    tiff(paste("Results/PAAD/Network/network_",chainType,"_pseudonormal_",i,".tiff",sep=""),res=300,h=3000,w=3000)
    plot(net,vertex.label=NA,layout=layout_with_graphopt(net,niter=800,charge=0.01))
    dev.off()
  }
}

########
####Clonal expansion
######
tiff("Results/PAAD/Network/network_clonal_expansion_IGH.tiff",h=2000,w=2000,res=300)
ggboxplot(PAAD.repertoire.diversity.filter, x = "Tumor_type_4categ", y = "clonal_expansion_IGH",color = "Tumor_type_4categ",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=clonal_expansion_IGH, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))
dev.off()

tiff("Results/PAAD/Network/network_clonal_expansion_IGK.tiff",h=2000,w=2000,res=300)
ggboxplot(PAAD.repertoire.diversity.filter, x = "Tumor_type_4categ", y = "clonal_expansion_IGK",color = "Tumor_type_4categ",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=clonal_expansion_IGK, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))
dev.off()


tiff("Results/PAAD/Network/network_clonal_expansion_IGL.tiff",h=2000,w=2000,res=300)
ggboxplot(PAAD.repertoire.diversity.filter, x = "Tumor_type_4categ", y = "clonal_expansion_IGL",color = "Tumor_type_4categ",ggtheme = theme_bw()) +
  rotate_x_text() +
  geom_point(aes(x=Tumor_type_4categ, y=clonal_expansion_IGL, color=Tumor_type_4categ), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("normal pancreas","PAC-Other","PDAC", "pseudonormal pancreas")) +
  stat_compare_means(
    comparisons =list(c("normal_pancreas","PDAC"),c("normal_pancreas","pseudonormal_pancreas"),c("pseudonormal_pancreas","PDAC"),c("PDAC","PAC-Other")))
dev.off()
