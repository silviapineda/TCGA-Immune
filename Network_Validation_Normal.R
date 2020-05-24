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
load("Data/Validation_Normal_pancreas/Pancreas_Normal_Validation_FullData.Rdata")
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
    write.table(get(paste0("edges",i)),paste0("Results/Validation/Normal/Network/edges_",chainType,"_",i,".txt"),sep="\t",row.names = F)
    write.table(get(paste0("vertex",i)),paste0("Results/Validation/Normal/Network/vertex_",chainType,"_",i,".txt"),sep="\t",row.names = F)
  }
}

network_IGH<-Obtain_vertex_edges(data_merge,"IGH")
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
Obtain_gini_index<-function(data,chainType,Pancreas.Normal.Validation.repertoire.diversity){
  sample<-rownames(Pancreas.Normal.Validation.repertoire.diversity)
  
  vertex_max<-NULL
  vertex_gini<-NULL
  cluster_max<-NULL
  cluster_gini<-NULL
  num_reads_max_cluster<-NULL
  clusters<-NULL
  j<-1
  
  for (i in sample){
    print(i)
    assign(paste0("edges",i),read.delim(paste0("Results/Validation/Normal/Network/edges_",chainType,"_",i,".txt")))
    assign(paste0("vertex",i),read.delim(paste0("Results/Validation/Normal/Network/vertex_",chainType,"_",i,".txt")))
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


chainType= "IGL"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,Pancreas.Normal.Validation.repertoire.diversity)))
id<-match(rownames(cluster_gini_IGL),rownames(Pancreas.Normal.Validation.repertoire.diversity))
Pancreas.Normal.Validation.repertoire.diversity[id,paste0("cluster_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_gini"]
Pancreas.Normal.Validation.repertoire.diversity[id,paste0("vertex_gini_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_gini"]
Pancreas.Normal.Validation.repertoire.diversity[id,paste0("vertex_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"vertex_max"]
Pancreas.Normal.Validation.repertoire.diversity[id,paste0("cluster_max_",chainType)]<-get(paste0("cluster_gini_",chainType))[,"cluster_max"]

##To save the varibles with vertex and cluster
save(data_merge,Pancreas.Normal.Validation.repertoire.diversity,file="Data/Validation_Normal_pancreas//Pancreas_Normal_Validation_FullData.Rdata")


########################
##4.Plot the network
########################
#Normal
sample_normal<-Pancreas.Normal.Validation.repertoire.diversity$sample
chainType="IGH"
for(i in sample_normal) {
  print(i)
  edges <- read.delim(paste("Results/Validation/Normal/Network/edges_",chainType,"_",i,".txt.outcome.txt",sep = ""))
  vertex <- read.delim(paste("Results/Validation/Normal/Network/vertex_",chainType,"_",i,".txt",sep = ""))
  #if(length(edges$edge1)!=0){
    net<-graph_from_data_frame(d=edges,vertices = vertex,directed=F)
    V(net)$size <- V(net)$Freq/2
    V(net)$color <- c( "#FDC086")
    net <- simplify(net, remove.multiple = F, remove.loops = T) 
    E(net)$arrow.mode <- 0
    E(net)$width <- 0.4
    E(net)$color <- c("black")
    tiff(paste("Results/Validation/Normal/Network/network_",chainType,"_",i,".tiff",sep=""),res=300,h=3000,w=3000)
    plot(net,vertex.label=NA,layout=layout_with_graphopt(net,niter=800,charge=0.01))
    dev.off()
  
  #}
}

###Plot all the samples together
load("Data/PAAD/PAAD_FullData.Rdata")
load("Data/GTEx/Pancreas/GTEx_FullData.Rdata")
load("Data/Pancreas_Validation/Pancreas_Validation_FullData.Rdata")

Cluster_vertex_gini_distribution<-rbind(PAAD.repertoire.diversity[,c("cluster_gini_IGH","vertex_gini_IGH")],
                                        Pancreas.repertoire.diversity[,c("cluster_gini_IGH","vertex_gini_IGH")],
                                        Pancreas.Validation.repertoire.diversity[,c("cluster_gini_IGH","vertex_gini_IGH")],
                                        Pancreas.Normal.Validation.repertoire.diversity[,c("cluster_gini_IGH","vertex_gini_IGH")])
Cluster_vertex_gini_distribution$outcome<-c(as.character(PAAD.repertoire.diversity$Tumor_type_4categ),rep("GTEX_normal_pancreas",nrow(Pancreas.repertoire.diversity)),
                                            as.character(Pancreas.Validation.repertoire.diversity$tissue),rep("Validation_normal_pancreas",nrow(Pancreas.Normal.Validation.repertoire.diversity)))

Cluster_vertex_gini_distribution<-Cluster_vertex_gini_distribution[which(Cluster_vertex_gini_distribution$outcome!="PAC-Other" & 
                                                                           Cluster_vertex_gini_distribution$outcome!="pseudonormal_pancreas" &
                                                                           Cluster_vertex_gini_distribution$outcome!="normal pancreas" &
                                                                           Cluster_vertex_gini_distribution$outcome!="normal_pancreas"),]
Cluster_vertex_gini_distribution$outcome<-factor(Cluster_vertex_gini_distribution$outcome)
Cluster_vertex_gini_distribution.filter<-Cluster_vertex_gini_distribution
tiff(paste0("Results/network_vertex_cluster_gini_",chainType,"_ALL.tiff"),h=2000,w=2000,res=300)

#cols=brewer.pal(6,name = "Pastel1")
cols= c("#7FC97F","#FBB4AE", "#BEAED4" ,  "#FDC086")

par(fig=c(0,0.8,0,0.8))
plot(Cluster_vertex_gini_distribution.filter[,paste0("cluster_gini_",chainType)], 
     Cluster_vertex_gini_distribution.filter[,paste0("vertex_gini_",chainType)],
     col = cols[factor(Cluster_vertex_gini_distribution.filter$outcome)],
     pch=20,ylab = c("Gini (Vextex)"),xlab = c("Gini (Cluster)"))
legend("bottomright",legend=c("GTEX_normal_pancreas (n=180)",
                              "Validation_PDAC (n=10)",
                              "TCGA_PDAC (n=131)",
                              "Validation_normal_pancreas (n=4)"), 
       col=cols,pch=20,cex=0.8)

par(fig=c(0,0.8,0.55,1), new=TRUE)
summary(glm(Cluster_vertex_gini_distribution.filter[,paste0("cluster_gini_",chainType)]~Cluster_vertex_gini_distribution.filter$outcome))
boxplot(Cluster_vertex_gini_distribution.filter[,paste0("cluster_gini_",chainType)]~Cluster_vertex_gini_distribution.filter$outcome,
        col=cols, horizontal=TRUE, axes=FALSE)



par(fig=c(0.65,1,0,0.8),new=TRUE)
summary(glm(Cluster_vertex_gini_distribution.filter[,paste0("vertex_gini_",chainType)]~Cluster_vertex_gini_distribution.filter$outcome))
boxplot(Cluster_vertex_gini_distribution.filter[,paste0("vertex_gini_",chainType)]~Cluster_vertex_gini_distribution.filter$outcome,
        col=cols, axes=FALSE)

dev.off()

tiff(paste0("Results/network_cluster_gini_",chainType,"_ALL.tiff"),h=2000,w=2000,res=300)
ggboxplot(Cluster_vertex_gini_distribution.filter, x = "outcome" , y =  "cluster_gini_IGH",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y= cluster_gini_IGH, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX_normal_pancreas (n=180)",
                                                  "Validation_PDAC (n=10)",
                                                  "TCGA_PDAC (n=131)",
                                                  "Validation_normal_pancreas (n=4)")) +
  
  stat_compare_means(
    comparisons =list(c("GTEX_normal_pancreas","Validation_normal_pancreas"),c("GTEX_normal_pancreas","PDAC"),
                      c("PDAC","pancreas tumor"),c("Validation_normal_pancreas","pancreas tumor"),c("PDAC","Validation_normal_pancreas")))

dev.off()

tiff(paste0("Results/network_vertex_gini_",chainType,"_ALL.tiff"),h=2000,w=2000,res=300)
ggboxplot(Cluster_vertex_gini_distribution.filter, x = "outcome" , y =  "vertex_gini_IGH",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y= vertex_gini_IGH, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX_normal_pancreas (n=180)",
                                                  "Validation_PDAC (n=10)",
                                                  "TCGA_PDAC (n=131)",
                                                  "Validation_normal_pancreas (n=4)")) +
  
  stat_compare_means(
    comparisons =list(c("GTEX_normal_pancreas","Validation_normal_pancreas"),c("GTEX_normal_pancreas","PDAC"),
                      c("PDAC","pancreas tumor"),c("Validation_normal_pancreas","pancreas tumor"),c("PDAC","Validation_normal_pancreas")))
dev.off()
