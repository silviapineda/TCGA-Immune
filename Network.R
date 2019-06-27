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

PAAD.repertoire.diversity$Tumor_type_2categ<-ifelse(PAAD.repertoire.diversity$Tumor_type=="Tumor_pancreas","Tumor_pancres",
                                                    ifelse(PAAD.repertoire.diversity$Tumor_type=="Solid_tissue_normal","normal_pseudonormal_pancreas",
                                                           ifelse(PAAD.repertoire.diversity$Tumor_type=="Adjacent_normal_pancreas","normal_pseudonormal_pancreas",
                                                                  ifelse(PAAD.repertoire.diversity$Tumor_type=="Pseudonormal (<1% neoplastic cellularity)","normal_pseudonormal_pancreas",NA))))
PAAD.repertoire.diversity$Tumor_type_2categ<-as.factor(PAAD.repertoire.diversity$Tumor_type_2categ)
PAAD.repertoire.diversity<-PAAD.repertoire.diversity[which(is.na(PAAD.repertoire.diversity$Tumor_type_2categ)==F),]


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
    V(net)$color <- c("#FDC086")
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
  #id<-grep("c174b41a-84c4-4a33-9c21-48dee5029ddb",sample) #IGKV
  #sample<-sample[-id]
  id<-grep("cab8e91a-ca41-4c62-af53-3aa4057d68d5",sample) #IGHV
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
# PAAD.repertoire.diversity.gini_IG_tumor<-PAAD.repertoire.diversity.gini_IG[which(PAAD.repertoire.diversity.gini_IG$Tumor_type=="Tumor_pancreas"),]
# write.csv(get(paste0("PAAD.repertoire.diversity.gini_",receptor,"_tumor")),file=paste0("Results/Network/network_",receptor,"_tumor.csv"))
brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")

tiff(paste0("Results/Network/network_vertex_cluster_gini_",receptor,".tiff"),h=2000,w=2000,res=300)
plot(get(paste0("PAAD.repertoire.diversity.gini_",receptor))[,"cluster_gini"], 
     get(paste0("PAAD.repertoire.diversity.gini_",receptor))[,"vertex_gini"],
     col = COLOR[get(paste0("PAAD.repertoire.diversity.gini_",chainType))[,"Tumor_type_2categ"]],pch=20,ylab = "Gini (Vextex)",xlab = "Gini (Cluster)")
legend("topleft",legend=c("Normal_pseudonormal_pancreas","Tumor_pancreas"), col=COLOR,pch=19,cex=0.8)
dev.off()


chainType= "IGHV"
assign(paste0("cluster_gini_",chainType),data.frame(Obtain_gini_index(data_merge,chainType,PAAD.repertoire.diversity)))
assign(paste0("PAAD.repertoire.diversity.gini_",chainType),merge(PAAD.repertoire.diversity,get(paste0("cluster_gini_",chainType)),by="row.names"))
#write.csv(get(paste0("PAAD.repertoire.diversity.gini_",chainType)),file=paste0("Results/Network/network",chainType,".csv"))
###only for tumor
#assign(paste0("PAAD.repertoire.diversity.gini_",chainType,"_tumor"),
#       get(paste0("PAAD.repertoire.diversity.gini_",chainType))[which(get(paste0("PAAD.repertoire.diversity.gini_",chainType))$Tumor_type=="Tumor_pancreas"),])
#write.csv(get(paste0("PAAD.repertoire.diversity.gini_",chainType,"_tumor")),file=paste0("Results/Network/network_",chainType,"_tumor.csv"))
tiff(paste0("Results/Network/network_vertex_cluster_gini_",chainType,".tiff"),h=2000,w=2000,res=300)
brewer.pal(n = 3, name = "Accent")
COLOR=c("#BEAED4","#7FC97F")
plot(get(paste0("PAAD.repertoire.diversity.gini_",chainType))[,"cluster_gini"], 
     get(paste0("PAAD.repertoire.diversity.gini_",chainType))[,"vertex_gini"],
     col = COLOR[get(paste0("PAAD.repertoire.diversity.gini_",chainType))[,"Tumor_type_2categ"]],pch=19,ylab = c("Gini (Vextex)"),xlab = c("Gini (Cluster)"))
    legend("topleft",legend=c("Normal_pseudonormal_pancreas","Tumor_pancreas"), col=COLOR,pch=19,cex=0.8)
dev.off()


###Info for the three samples with very high gini(vertex) and gini(cluster)
samples_high<-c("5e9e81e2-0e8b-4aca-aced-6ce451fa3262", "08c4c14a-94a8-4063-bbc3-916579960078", "49170bdc-8725-43a0-b2d7-ccb52fa415e6")
samples_high_TCGA<-substr(as.character(PAAD.repertoire.diversity.gini_IG[match(samples_high,PAAD.repertoire.diversity.gini_IG$Row.names),"TCGA_sample"]),1,12)

clinical.patient[match(samples_high_TCGA,clinical.patient$bcr_patient_barcode),]

#################################################
#### Analysis of the GINI vertex and cluster ####
#################################################
PAAD.repertoire.diversity.gini_IGHV<-read.csv("Results/Network/network_IGHV_tumor.csv")
##Delete outliers
PAAD.repertoire.diversity.gini_IGHV<-PAAD.repertoire.diversity.gini_IGHV[which(PAAD.repertoire.diversity.gini_IGHV$vertex_gini<0.8),]
clinical.Gini<-clinical.patient[match(substr(PAAD.repertoire.diversity.gini_IGHV$TCGA_sample,1,12),clinical.patient$bcr_patient_barcode),]

##Function to run the association between clinical outcome and BCR/TCR
association.test.immuneRep<- function (clinical.patient.tumor,clinical.var,PAAD.repertoire.tumor,markers){
  if (class(clinical.patient.tumor[,clinical.var])=="factor"){
    p.markers<-NULL
    for(i in 1:length(markers)){
      p.markers[i]<-kruskal.test(PAAD.repertoire.tumor[which(clinical.patient.tumor[,clinical.var]!=""),markers[i]]~
                                   clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var])$p.value
    }
    dfplot <- data.frame(PAAD.repertoire.tumor[which(clinical.patient.tumor[,clinical.var]!=""),markers[which(p.markers<0.05)]],
                         clinical.var = clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var])
    colnames(dfplot)[ncol(dfplot)]<- clinical.var
    
    if (dim(dfplot)[2]>1){
      for(i in 1:length(markers[which(p.markers<0.05)])){
        print(i)
        dfplot$marker<-PAAD.repertoire.tumor[which(clinical.patient.tumor[,clinical.var]!=""),markers[which(p.markers<0.05)][i]]
        tiff(paste0("Results/boxplot_",clinical.var,"_",markers[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
        print(ggplot(dfplot,aes(y=marker, fill=clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var], 
                                x=clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var])) + geom_boxplot() 
              + scale_y_continuous(name= markers[which(p.markers<0.05)][i]) + stat_compare_means()  + scale_x_discrete(name=clinical.var)
              + scale_fill_discrete(name=clinical.var))
        dev.off()
      }
    }
  } else {
    p.markers<-NULL
    for(i in 1:length(markers)){
      p.markers[i]<-coef(summary(glm(PAAD.repertoire.tumor[which(clinical.patient.tumor[,clinical.var]!=""),markers[i]]~
                                       clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var])))[2,4]
    }
    dfplot <- data.frame(marker =PAAD.repertoire.tumor[which(clinical.patient.tumor[,clinical.var]!=""),markers[which(p.markers<0.05)]],
                         clinical.var = clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var])
    colnames(dfplot)[ncol(dfplot)]<- clinical.var
    
    if (dim(dfplot)[2]>1){
      for(i in 1:length(markers[which(p.markers<0.05)])){
        print(i)
        dfplot$marker<-PAAD.repertoire.tumor[which(clinical.patient.tumor[,clinical.var]!=""),markers[which(p.markers<0.05)][i]]
        tiff(paste0("Results/boxplot_",clinical.var,"_",markers[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
        print(ggplot(dfplot,aes(y=marker, fill=clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var], 
                                x=clinical.patient.tumor[which(clinical.patient.tumor[,clinical.var]!=""),clinical.var])) + geom_point() + geom_smooth(method='lm')
              + scale_y_continuous(name= markers[which(p.markers<0.05)][i]) + stat_cor(method = "pearson") + scale_x_continuous(name=clinical.var)
              + scale_fill_continuous(name=clinical.var))
        dev.off()
      }
    }
  }
  return("Done")
}

markers<-c("cluster_gini","vertex_gini")
##Histological type
clinical.Gini$histological_type_2cat<-factor(ifelse(clinical.Gini$histological_type=="Pancreas-Adenocarcinoma Ductal Type","PDAC","Other"))
clinical.Gini$histological_type_2cat<-factor(clinical.Gini$histological_type_2cat)
association.test.immuneRep(clinical.Gini,"histological_type_2cat",PAAD.repertoire.diversity.gini_IGHV,markers)

##anatomic_neoplasm_subdivision
clinical.Gini$anatomic_neoplasm_subdivision<-factor(clinical.Gini$anatomic_neoplasm_subdivision)
association.test.immuneRep(clinical.Gini,"anatomic_neoplasm_subdivision",PAAD.repertoire.diversity.gini_IGHV,markers)

##gender
clinical.Gini$gender<-factor(clinical.Gini$gender)
association.test.immuneRep(clinical.Gini,"gender",PAAD.repertoire.diversity.gini_IGHV,markers)

##race_list
clinical.Gini$race_list<-factor(clinical.Gini$race_list)
association.test.immuneRep(clinical.Gini,"race_list",PAAD.repertoire.diversity.gini_IGHV,markers)

##History of Prior Malignancy
clinical.Gini$other_dx<-factor(clinical.Gini$other_dx)
association.test.immuneRep(clinical.Gini,"other_dx",PAAD.repertoire.diversity.gini_IGHV,markers)

##lymph_node_examined_count
association.test.immuneRep(clinical.Gini,"lymph_node_examined_count",PAAD.repertoire.diversity.gini_IGHV,markers)

##neoplasm_histologic_grade
clinical.Gini$neoplasm_histologic_grade_3cat<-factor(ifelse(clinical.Gini$neoplasm_histologic_grade=="G1","G1",
                                                                     ifelse(clinical.Gini$neoplasm_histologic_grade=="G2","G2",
                                                                            ifelse(clinical.Gini$neoplasm_histologic_grade=="G3","G3",""))))
association.test.immuneRep(clinical.Gini,"neoplasm_histologic_grade_3cat",PAAD.repertoire.diversity.gini_IGHV,markers)

##Age 
association.test.immuneRep(clinical.Gini,"age_at_initial_pathologic_diagnosis",PAAD.repertoire.diversity.gini_IGHV,markers)


##Smoking
clinical.Gini$smoking<-factor(ifelse(clinical.Gini$tobacco_smoking_history_master=="Current smoker (includes daily smokers and non-daily smokers or occasional smokers)","Current",
                                              ifelse(clinical.Gini$tobacco_smoking_history_master=="Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime)","Non-smoker","Former")))
association.test.immuneRep(clinical.Gini,"smoking",PAAD.repertoire.diversity.gini_IGHV,markers)

##number_pack_years_smoked
association.test.immuneRep(clinical.Gini,"number_pack_years_smoked",PAAD.repertoire.diversity.gini_IGHV,markers)

##Alcohol
clinical.Gini$alcohol_history_documented<-factor(clinical.Gini$alcohol_history_documented)
association.test.immuneRep(clinical.Gini,"alcohol_history_documented",PAAD.repertoire.diversity.gini_IGHV,markers)

##Alcohol category
clinical.Gini$alcoholic_exposure_category2<-ifelse(clinical.Gini$alcohol_history_documented=="NO","No-drinker",
                                                            ifelse(clinical.Gini$alcohol_history_documented=="YES" & clinical.Gini$alcoholic_exposure_category=="",NA,
                                                                   ifelse(clinical.Gini$alcoholic_exposure_category=="None","None-Drinker",
                                                                          ifelse(clinical.Gini$alcoholic_exposure_category=="Occasional Drinker","Occasional-Drinker",
                                                                                 ifelse(clinical.Gini$alcoholic_exposure_category=="Daily Drinker","Daily-Drinker",
                                                                                        ifelse(clinical.Gini$alcoholic_exposure_category=="Social Drinker","Social-Drinker",
                                                                                               ifelse(clinical.Gini$alcoholic_exposure_category=="Weekly Drinker","Weekly-Drinker",NA)))))))
clinical.Gini$alcoholic_exposure_category2<-factor(clinical.Gini$alcoholic_exposure_category2)
association.test.immuneRep(clinical.Gini,"alcoholic_exposure_category2",PAAD.repertoire.diversity.gini_IGHV,markers)

##family history
clinical.Gini$family_history_of_cancer<-factor(clinical.Gini$family_history_of_cancer)
association.test.immuneRep(clinical.Gini,"family_history_of_cancer",PAAD.repertoire.diversity.gini_IGHV,markers)

##radiation_therapy
clinical.Gini$radiation_therapy<-factor(clinical.Gini$radiation_therapy)
association.test.immuneRep(clinical.Gini,"radiation_therapy",PAAD.repertoire.diversity.gini_IGHV,markers)

##primary_therapy_outcome_success
clinical.Gini$primary_therapy_outcome_success<-factor(clinical.Gini$primary_therapy_outcome_success)
association.test.immuneRep(clinical.Gini,"primary_therapy_outcome_success",PAAD.repertoire.diversity.gini_IGHV,markers)

##history_chronic_pancreatitis
clinical.Gini$history_of_chronic_pancreatitis<-factor(clinical.Gini$history_of_chronic_pancreatitis)
association.test.immuneRep(clinical.Gini,"history_of_chronic_pancreatitis",PAAD.repertoire.diversity.gini_IGHV,markers)


##################################
#######Clinical follow-up########
#################################
clinical.Gini.follow_up<-clinical.folow_up[match(substr(PAAD.repertoire.diversity.gini_IGHV$TCGA_sample,1,12),clinical.folow_up$bcr_patient_barcode),]
##vital_status
clinical.Gini.follow_up$vital_status<-factor(clinical.Gini.follow_up$vital_status)
association.test.immuneRep(clinical.Gini.follow_up,"vital_status",PAAD.repertoire.diversity.gini_IGHV,markers)

##new tumor event
clinical.Gini.follow_up$new_tumor_event_type<-replace(clinical.Gini.follow_up$new_tumor_event_type,clinical.Gini.follow_up$new_tumor_event_type=="#N/A",NA)
clinical.Gini.follow_up$new_tumor_event_type<-replace(clinical.Gini.follow_up$new_tumor_event_type,
                                                       clinical.Gini.follow_up$new_tumor_event_type=="Locoregional Recurrence|Distant Metastasis" | 
                                                         clinical.Gini.follow_up$new_tumor_event_type=="New Primary Tumor",NA)
clinical.Gini.follow_up$new_tumor_event_type<-factor(clinical.Gini.follow_up$new_tumor_event_type)
association.test.immuneRep(clinical.Gini.follow_up,"new_tumor_event_type",PAAD.repertoire.diversity.gini_IGHV,markers)

#treatment_outcome_first_course
clinical.Gini.follow_up$treatment_outcome_first_course<-replace(clinical.Gini.follow_up$treatment_outcome_first_course,
                                                                 clinical.Gini.follow_up$treatment_outcome_first_course=="[Discrepancy]" |
                                                                   clinical.Gini.follow_up$treatment_outcome_first_course=="[Not Applicable]" | 
                                                                   clinical.Gini.follow_up$treatment_outcome_first_course=="[Not Available]" |
                                                                   clinical.Gini.follow_up$treatment_outcome_first_course=="[Unknown]",NA)
clinical.Gini.follow_up$treatment_outcome_first_course<-factor(clinical.Gini.follow_up$treatment_outcome_first_course)
association.test.immuneRep(clinical.Gini.follow_up,"treatment_outcome_first_course",PAAD.repertoire.diversity.gini_IGHV,markers)


##################################
### Survival Analysis ############
##################################
library(survival)
library(survminer)
library(survMisc)

##OS
surv_object <- Surv(time = clinical.Gini.follow_up$OS.time, event = clinical.Gini.follow_up$OS)
res.cox <- coxph(surv_object~PAAD.repertoire.diversity.gini_IGHV$KappaLambda_ratio_expression)
summary(res.cox)
##Categorical
KL_mean<-mean(PAAD.repertoire.diversity.gini_IGHV$KappaLambda_ratio_expression)
PAAD.repertoire.diversity.gini_IGHV$KL_ratio_2cat<-ifelse(PAAD.repertoire.diversity.gini_IGHV$KappaLambda_ratio_expression<=KL_mean,1,2)
fit1 <- survfit(surv_object ~ PAAD.repertoire.diversity.gini_IGHV$KL_ratio_2cat)
fit1
ggsurvplot(fit1, data = PAAD.repertoire.diversity.gini_IGHV)
comp(ten(fit1))$tests$lrTests

#DSS
clinical.Gini.follow_up2<-clinical.Gini.follow_up[which(clinical.Gini.follow_up$DSS!="#N/A"),]
clinical.Gini.follow_up2$DSS<-as.integer(as.character(clinical.Gini.follow_up2$DSS))
PAAD.repertoire.diversity.gini_IGHV2<-PAAD.repertoire.diversity.gini_IGHV[which(clinical.Gini.follow_up$DSS!="#N/A"),]
surv_object <- Surv(time = clinical.Gini.follow_up2$DSS.time, event = clinical.Gini.follow_up2$DSS)
res.cox <- coxph(surv_object~PAAD.repertoire.diversity.gini_IGHV2$KappaLambda_ratio_expression)
summary(res.cox)
#Categorical
KL_mean<-mean(PAAD.repertoire.diversity.gini_IGHV2$KappaLambda_ratio_expression)
PAAD.repertoire.diversity.gini_IGHV2$KL_ratio_2cat<-ifelse(PAAD.repertoire.diversity.gini_IGHV2$KappaLambda_ratio_expression<=KL_mean,1,2)
fit1 <- survfit(surv_object ~ PAAD.repertoire.diversity.gini_IGHV2$KL_ratio_2cat)
fit1
ggsurvplot(fit1, data = PAAD.repertoire.diversity.gini_IGHV2)
comp(ten(fit1))$tests$lrTests

##PFI
surv_object <- Surv(time = clinical.Gini.follow_up$PFI.time, event = clinical.Gini.follow_up$PFI)
res.cox <- coxph(surv_object~PAAD.repertoire.diversity.gini_IGHV$KappaLambda_ratio_expression)
summary(res.cox)
