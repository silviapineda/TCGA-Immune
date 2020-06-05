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
library(reshape)

setwd("~/TCGA-Immune/")

##################################
##1.Obtain the vertex and edges
#################################
#chainType
Obtain_vertex_edges<-function(data,chainType,sample_order){
  sample<-levels(factor(data$sample))
  data<-data[which(data$chainType==chainType),]
  data$CloneId_CDR3<-paste0(data[,c("V_J_lenghCDR3_CloneId")],data[,c("nSeqCDR3")])
  #sample<-unique(data$sample)
  for (i in sample){
    print(i)
    data_sample<-data[which(data$sample==i),]
    ##% of subsample 
    data_subsample<-data_sample[sample(dim(data_sample)[1], round(dim(data_sample)[1]*sample_order)),]
    df_sample<-data_subsample[,c("CloneId_CDR3","V_J_lenghCDR3_CloneId")]
    groups <- group_by(df_sample,CloneId_CDR3,V_J_lenghCDR3_CloneId)
    assign(paste0("edges",i),unique(data.frame(groups)))
    df_vertex<-data.frame(table(data_subsample$CloneId_CDR3))
    assign(paste0("vertex",i),df_vertex[which(df_vertex$Freq!=0),])
    write.table(get(paste0("edges",i)),paste0("Results/Subsample/edges_",chainType,"_",i,"_",sample_order,".txt"),sep="\t",row.names = F)
    write.table(get(paste0("vertex",i)),paste0("Results/Subsample/vertex_",chainType,"_",i,"_",sample_order,".txt"),sep="\t",row.names = F)
  }
}


###############################################
##3. Obtain the vertex and cluster Gini index
###############################################
#chainType
Obtain_gini_index<-function(data,chainType,PAAD.repertoire.diversity,sample_order){
  sample<-rownames(PAAD.repertoire.diversity)
  
  vertex_max<-NULL
  vertex_gini<-NULL
  cluster_max<-NULL
  cluster_gini<-NULL
  num_reads_max_cluster<-NULL
  clusters<-NULL
  j<-1
  for (i in sample){
    print(i)
    assign(paste0("edges",i),read.delim(paste0("Results/Subsample/edges_",chainType,"_",i,"_",sample_order,".txt")))
    assign(paste0("vertex",i),read.delim(paste0("Results/Subsample/vertex_",chainType,"_",i,"_",sample_order,".txt")))
    if(dim(get(paste0("edges",i)))[1]>0){
      vertex_max[j]<-max(get(paste0("vertex",i))$Freq)
      vertex_gini[j]<-Gini(get(paste0("vertex",i))$Freq)
      cluster_max[j]<-max(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId))
      clusters[j]<-sum(table(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId)))
      num_reads_max_cluster[j]<-tail(table(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId)),1)
      cluster_gini[j]<-Gini(table(get(paste0("edges",i))$V_J_lenghCDR3_CloneId))
    }
    j=j+1
  }
  results<-cbind(cluster_gini,vertex_gini,vertex_max,cluster_max,num_reads_max_cluster,clusters)
  
  rownames(results)<-sample
  
  return(results)
}

##########################
### Subsample PAAD TCGA 
##########################
load("Data/PAAD/PAAD_FullData.Rdata")
##20
chainType= "IGH"
for(i in 1:5){
  network_IGH<-Obtain_vertex_edges(data_merge,"IGH",0.2)
  assign(paste0("cluster_gini_",i,"_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity,0.2)))
}

id<-match(rownames(get(paste0("cluster_gini_",i,"_",chainType))),rownames(PAAD.repertoire.diversity))
PAAD.repertoire.diversity[id,paste0("cluster_gini_20_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"cluster_gini"],
                                                                         get(paste0("cluster_gini_2_",chainType))[,"cluster_gini"],
                                                                         get(paste0("cluster_gini_3_",chainType))[,"cluster_gini"],
                                                                         get(paste0("cluster_gini_4_",chainType))[,"cluster_gini"],
                                                                         get(paste0("cluster_gini_5_",chainType))[,"cluster_gini"]))
PAAD.repertoire.diversity[id,paste0("vertex_gini_20_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_2_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_3_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_4_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_5_",chainType))[,"vertex_gini"]))


##40
chainType= "IGH"
for(i in 1:5){
  network_IGH<-Obtain_vertex_edges(data_merge,"IGH",0.4)
  assign(paste0("cluster_gini_",i,"_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity,0.4)))
}

id<-match(rownames(get(paste0("cluster_gini_",i,"_",chainType))),rownames(PAAD.repertoire.diversity))
PAAD.repertoire.diversity[id,paste0("cluster_gini_40_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_2_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_3_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_4_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_5_",chainType))[,"cluster_gini"]))
PAAD.repertoire.diversity[id,paste0("vertex_gini_40_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_2_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_3_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_4_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_5_",chainType))[,"vertex_gini"]))

##60
chainType= "IGH"
for(i in 1:5){
  network_IGH<-Obtain_vertex_edges(data_merge,"IGH",0.6)
  assign(paste0("cluster_gini_",i,"_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity,0.6)))
}

id<-match(rownames(get(paste0("cluster_gini_",i,"_",chainType))),rownames(PAAD.repertoire.diversity))
PAAD.repertoire.diversity[id,paste0("cluster_gini_60_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_2_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_3_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_4_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_5_",chainType))[,"cluster_gini"]))
PAAD.repertoire.diversity[id,paste0("vertex_gini_60_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_2_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_3_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_4_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_5_",chainType))[,"vertex_gini"]))

##80
chainType= "IGH"
for(i in 1:5){
  network_IGH<-Obtain_vertex_edges(data_merge,"IGH",0.8)
  assign(paste0("cluster_gini_",i,"_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity,0.8)))
}

id<-match(rownames(get(paste0("cluster_gini_",i,"_",chainType))),rownames(PAAD.repertoire.diversity))
PAAD.repertoire.diversity[id,paste0("cluster_gini_80_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_2_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_3_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_4_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_5_",chainType))[,"cluster_gini"]))
PAAD.repertoire.diversity[id,paste0("vertex_gini_80_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_2_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_3_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_4_",chainType))[,"vertex_gini"],
                                                                                  get(paste0("cluster_gini_5_",chainType))[,"vertex_gini"]))

##100
chainType= "IGH"
for(i in 1:5){
  network_IGH<-Obtain_vertex_edges(data_merge,"IGH",1)
  assign(paste0("cluster_gini_",i,"_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity,1)))
}

id<-match(rownames(get(paste0("cluster_gini_",i,"_",chainType))),rownames(PAAD.repertoire.diversity))
PAAD.repertoire.diversity[id,paste0("cluster_gini_100_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_2_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_3_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_4_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_5_",chainType))[,"cluster_gini"]))
PAAD.repertoire.diversity[id,paste0("vertex_gini_100_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"vertex_gini"],
                                                                                   get(paste0("cluster_gini_2_",chainType))[,"vertex_gini"],
                                                                                   get(paste0("cluster_gini_3_",chainType))[,"vertex_gini"],
                                                                                   get(paste0("cluster_gini_4_",chainType))[,"vertex_gini"],
                                                                                   get(paste0("cluster_gini_5_",chainType))[,"vertex_gini"]))

##To save the varibles with vertex and cluster
save(PAAD.repertoire.diversity,file="Data/PAAD/PAAD_Subsample.Rdata")

##########################
### Subsample Pancreas GTEX 
##########################
load("Data/GTEx/Pancreas/GTEx_FullData.Rdata")
##20
chainType= "IGH"
for(i in 1:5){
  network_IGH<-Obtain_vertex_edges(data_merge_pancreas,"IGH",0.2)
  assign(paste0("cluster_gini_",i,"_",chainType),data.frame(Obtain_gini_index(data_merge_pancreas,chainType,Pancreas.repertoire.diversity,0.2)))
}

id<-match(rownames(get(paste0("cluster_gini_",i,"_",chainType))),rownames(Pancreas.repertoire.diversity))
Pancreas.repertoire.diversity[id,paste0("cluster_gini_20_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_2_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_3_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_4_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_5_",chainType))[,"cluster_gini"]))
Pancreas.repertoire.diversity[id,paste0("vertex_gini_20_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_2_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_3_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_4_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_5_",chainType))[,"vertex_gini"]))

##40
chainType= "IGH"
for(i in 1:5){
  network_IGH<-Obtain_vertex_edges(data_merge_pancreas,"IGH",0.4)
  assign(paste0("cluster_gini_",i,"_",chainType),data.frame(Obtain_gini_index(data_merge_pancreas,chainType,Pancreas.repertoire.diversity,0.4)))
}

id<-match(rownames(get(paste0("cluster_gini_",i,"_",chainType))),rownames(Pancreas.repertoire.diversity))
Pancreas.repertoire.diversity[id,paste0("cluster_gini_40_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_2_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_3_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_4_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_5_",chainType))[,"cluster_gini"]))
Pancreas.repertoire.diversity[id,paste0("vertex_gini_40_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_2_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_3_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_4_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_5_",chainType))[,"vertex_gini"]))

##60
chainType= "IGH"
for(i in 1:5){
  network_IGH<-Obtain_vertex_edges(data_merge_pancreas,"IGH",0.6)
  assign(paste0("cluster_gini_",i,"_",chainType),data.frame(Obtain_gini_index(data_merge_pancreas,chainType,Pancreas.repertoire.diversity,0.6)))
}

id<-match(rownames(get(paste0("cluster_gini_",i,"_",chainType))),rownames(Pancreas.repertoire.diversity))
Pancreas.repertoire.diversity[id,paste0("cluster_gini_60_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_2_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_3_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_4_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_5_",chainType))[,"cluster_gini"]))
Pancreas.repertoire.diversity[id,paste0("vertex_gini_60_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_2_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_3_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_4_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_5_",chainType))[,"vertex_gini"]))
##80
chainType= "IGH"
for(i in 1:5){
  network_IGH<-Obtain_vertex_edges(data_merge_pancreas,"IGH",0.8)
  assign(paste0("cluster_gini_",i,"_",chainType),data.frame(Obtain_gini_index(data_merge_pancreas,chainType,Pancreas.repertoire.diversity,0.8)))
}

id<-match(rownames(get(paste0("cluster_gini_",i,"_",chainType))),rownames(Pancreas.repertoire.diversity))
Pancreas.repertoire.diversity[id,paste0("cluster_gini_80_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_2_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_3_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_4_",chainType))[,"cluster_gini"],
                                                                                   get(paste0("cluster_gini_5_",chainType))[,"cluster_gini"]))
Pancreas.repertoire.diversity[id,paste0("vertex_gini_80_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_2_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_3_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_4_",chainType))[,"vertex_gini"],
                                                                                      get(paste0("cluster_gini_5_",chainType))[,"vertex_gini"]))

##100
chainType= "IGH"
for(i in 1:5){
  network_IGH<-Obtain_vertex_edges(data_merge_pancreas,"IGH",1)
  assign(paste0("cluster_gini_",i,"_",chainType),data.frame(Obtain_gini_index(data_merge_pancreas,chainType,Pancreas.repertoire.diversity,1)))
}

id<-match(rownames(get(paste0("cluster_gini_",i,"_",chainType))),rownames(Pancreas.repertoire.diversity))
Pancreas.repertoire.diversity[id,paste0("cluster_gini_100_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"cluster_gini"],
                                                                                    get(paste0("cluster_gini_2_",chainType))[,"cluster_gini"],
                                                                                    get(paste0("cluster_gini_3_",chainType))[,"cluster_gini"],
                                                                                    get(paste0("cluster_gini_4_",chainType))[,"cluster_gini"],
                                                                                    get(paste0("cluster_gini_5_",chainType))[,"cluster_gini"]))
Pancreas.repertoire.diversity[id,paste0("vertex_gini_100_",chainType)]<-rowMeans(cbind(get(paste0("cluster_gini_1_",chainType))[,"vertex_gini"],
                                                                                       get(paste0("cluster_gini_2_",chainType))[,"vertex_gini"],
                                                                                       get(paste0("cluster_gini_3_",chainType))[,"vertex_gini"],
                                                                                       get(paste0("cluster_gini_4_",chainType))[,"vertex_gini"],
                                                                                       get(paste0("cluster_gini_5_",chainType))[,"vertex_gini"]))

##To save the varibles with vertex and cluster
save(Pancreas.repertoire.diversity,file="Data/GTEx/Pancreas/GTEX_pancreas_Subsample.Rdata")


###Plot the results
#Cluster
cluster_gini_PAAD<-melt(PAAD.repertoire.diversity[,c("TCGA_sample","cluster_gini_20_IGH","cluster_gini_40_IGH","cluster_gini_60_IGH","cluster_gini_80_IGH","cluster_gini_100_IGH")])
cluster_gini_PAAD$percentage<-c(rep(20,dim(PAAD.repertoire.diversity)[1]),rep(40,dim(PAAD.repertoire.diversity)[1]),rep(60,dim(PAAD.repertoire.diversity)[1]),
                           rep(80,dim(PAAD.repertoire.diversity)[1]),rep(100,dim(PAAD.repertoire.diversity)[1]))
colnames(cluster_gini_PAAD)[1]<-"sample"

cluster_gini_GTEX<-melt(Pancreas.repertoire.diversity[,c("SUBJID","cluster_gini_20_IGH","cluster_gini_40_IGH","cluster_gini_60_IGH","cluster_gini_80_IGH","cluster_gini_100_IGH")])
cluster_gini_GTEX$percentage<-c(rep(20,dim(Pancreas.repertoire.diversity)[1]),rep(40,dim(Pancreas.repertoire.diversity)[1]),rep(60,dim(Pancreas.repertoire.diversity)[1]),
                                rep(80,dim(Pancreas.repertoire.diversity)[1]),rep(100,dim(Pancreas.repertoire.diversity)[1]))

colnames(cluster_gini_GTEX)[1]<-"sample"

cluster_gini<-rbind(cluster_gini_PAAD,cluster_gini_GTEX)
tiff("Results/ImmuneRep/Network/cluster_gini_subsampled.tiff",res=300,h=2000,w=2000)
ggplot(data = cluster_gini, aes(x=percentage,y=value,colour=sample))  + ylim(0,0.6) + geom_line() + geom_point() +
  theme(legend.position = "none") + scale_color_manual(values=c(rep("#BEAED4",dim(PAAD.repertoire.diversity)[1]),rep("#7FC97F",dim(Pancreas.repertoire.diversity)[1]))) +
  ylab("cluster_gini_index") + xlab("proportion of reads samples")
dev.off()

#Vertex
vertex_gini_PAAD<-melt(PAAD.repertoire.diversity[,c("TCGA_sample","vertex_gini_20_IGH","vertex_gini_40_IGH","vertex_gini_60_IGH","vertex_gini_80_IGH","vertex_gini_100_IGH")])
vertex_gini_PAAD$percentage<-c(rep(20,dim(PAAD.repertoire.diversity)[1]),rep(40,dim(PAAD.repertoire.diversity)[1]),rep(60,dim(PAAD.repertoire.diversity)[1]),
                                rep(80,dim(PAAD.repertoire.diversity)[1]),rep(100,dim(PAAD.repertoire.diversity)[1]))
colnames(vertex_gini_PAAD)[1]<-"sample"

vertex_gini_GTEX<-melt(Pancreas.repertoire.diversity[,c("SUBJID","vertex_gini_20_IGH","vertex_gini_40_IGH","vertex_gini_60_IGH","vertex_gini_80_IGH","vertex_gini_100_IGH")])
vertex_gini_GTEX$percentage<-c(rep(20,dim(Pancreas.repertoire.diversity)[1]),rep(40,dim(Pancreas.repertoire.diversity)[1]),rep(60,dim(Pancreas.repertoire.diversity)[1]),
                                rep(80,dim(Pancreas.repertoire.diversity)[1]),rep(100,dim(Pancreas.repertoire.diversity)[1]))

colnames(vertex_gini_GTEX)[1]<-"sample"

vertex_gini<-rbind(vertex_gini_PAAD,vertex_gini_GTEX)
tiff("Results/ImmuneRep/Network/vertex_gini_subsampled.tiff",res=300,h=2000,w=2000)
ggplot(data = vertex_gini, aes(x=percentage,y=value,colour=sample))  + ylim(0,1) + geom_line() + geom_point() +
  theme(legend.position = "none") + scale_color_manual(values=c(rep("#BEAED4",dim(PAAD.repertoire.diversity)[1]),rep("#7FC97F",dim(Pancreas.repertoire.diversity)[1]))) +
  ylab("vertex_gini_index") + xlab("proportion of reads samples")
dev.off()
