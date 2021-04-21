rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. Descriptive analysis 
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: Downsampling 
###         
###
### Author: Silvia Pineda
### Date: April, 2019
############################################################################################


setwd("~/TCGA-Immune/")

load("Data/PAAD/PAAD_FullData.Rdata")
load("Data/GTEx/Pancreas/GTEx_FullData.Rdata")

##Diversity measures
Obtain_diversity<-function(data,diversity,chainType,sample_order){
  sample<-rownames(diversity)
  entropy<-NULL
  data<-data[which(data$chainType==chainType),]
  for (i in 1:length(sample)){
    print(i)
    
    data_sample_unique<-data[which(data$sample==sample[i]),]
    data_subsample<-data_sample_unique[sample(dim(data_sample_unique)[1], round(dim(data_sample_unique)[1]*sample_order)),]
    clones_sample<-data_subsample[,"V_J_lenghCDR3_CloneId"]
    
    fi<-as.numeric(table(clones_sample))/length(clones_sample)
    hi<-fi*log2(fi)
    entropy[i]=-sum(hi)
    
  }
  return(entropy)
}

chainType="IGH"
for (i in 1:5){
  assign(paste0("entropy_",i,"_",chainType),Obtain_diversity(data_merge,PAAD.repertoire.diversity,chainType,0.2))
}
PAAD.repertoire.diversity[,paste0("entropy_20_",chainType)]<-rowMeans(cbind(get(paste0("entropy_1_",chainType)),
                                                                            get(paste0("entropy_2_",chainType)),
                                                                            get(paste0("entropy_3_",chainType)),
                                                                            get(paste0("entropy_4_",chainType)),
                                                                            get(paste0("entropy_5_",chainType))))

for (i in 1:5){
  assign(paste0("entropy_",i,"_",chainType),Obtain_diversity(data_merge,PAAD.repertoire.diversity,chainType,0.4))
}
PAAD.repertoire.diversity[,paste0("entropy_40_",chainType)]<-rowMeans(cbind(get(paste0("entropy_1_",chainType)),
                                                                            get(paste0("entropy_2_",chainType)),
                                                                            get(paste0("entropy_3_",chainType)),
                                                                            get(paste0("entropy_4_",chainType)),
                                                                            get(paste0("entropy_5_",chainType))))

for (i in 1:5){
  assign(paste0("entropy_",i,"_",chainType),Obtain_diversity(data_merge,PAAD.repertoire.diversity,chainType,0.6))
}
PAAD.repertoire.diversity[,paste0("entropy_60_",chainType)]<-rowMeans(cbind(get(paste0("entropy_1_",chainType)),
                                                                            get(paste0("entropy_2_",chainType)),
                                                                            get(paste0("entropy_3_",chainType)),
                                                                            get(paste0("entropy_4_",chainType)),
                                                                            get(paste0("entropy_5_",chainType))))


for (i in 1:5){
  assign(paste0("entropy_",i,"_",chainType),Obtain_diversity(data_merge,PAAD.repertoire.diversity,chainType,0.8))
}
PAAD.repertoire.diversity[,paste0("entropy_80_",chainType)]<-rowMeans(cbind(get(paste0("entropy_1_",chainType)),
                                                                            get(paste0("entropy_2_",chainType)),
                                                                            get(paste0("entropy_3_",chainType)),
                                                                            get(paste0("entropy_4_",chainType)),
                                                                            get(paste0("entropy_5_",chainType))))


Obtain_diversity<-function(data,diversity,chainType,sample_order){
  sample<-rownames(diversity)
  entropy<-NULL
  data<-data[which(data$chainType==chainType),]
  for (i in 1:length(sample)){
    print(i)
    
    data_sample_unique<-data[which(data$sample==sample[i]),]
    #data_subsample<-data_sample_unique[sample(dim(data_sample_unique)[1], round(dim(data_sample_unique)[1]*sample_order)),]
    if( dim(data_sample_unique)[1]>10){
      data_subsample<-data_sample_unique[sample(dim(data_sample_unique)[1],10),]
      clones_sample<-data_subsample[,"V_J_lenghCDR3_CloneId"]
      
      fi<-as.numeric(table(clones_sample))/length(clones_sample)
      hi<-fi*log2(fi)
      entropy[i]=-sum(hi)
    }else {
      entropy[i]=NA
    }
  }
  return(entropy)
}

chainType="IGH"
for (i in 1:5){
  assign(paste0("entropy_",i,"_",chainType),Obtain_diversity(data_merge,PAAD.repertoire.tumor.filter,chainType,0.2))
}
PAAD.repertoire.tumor.filter[,paste0("entropy_50_",chainType)]<-rowMeans(cbind(get(paste0("entropy_1_",chainType)),
                                                                            get(paste0("entropy_2_",chainType)),
                                                                            get(paste0("entropy_3_",chainType)),
                                                                            get(paste0("entropy_4_",chainType)),
                                                                            get(paste0("entropy_5_",chainType))))

for (i in 1:5){
  assign(paste0("entropy_",i,"_",chainType),Obtain_diversity(data_merge_pancreas,Pancreas.repertoire.diversity,chainType,0.2))
}
Pancreas.repertoire.diversity[,paste0("entropy_50_",chainType)]<-rowMeans(cbind(get(paste0("entropy_1_",chainType)),
                                                                               get(paste0("entropy_2_",chainType)),
                                                                               get(paste0("entropy_3_",chainType)),
                                                                               get(paste0("entropy_4_",chainType)),
                                                                               get(paste0("entropy_5_",chainType))))

boxplot(PAAD.repertoire.tumor.filter$entropy_50_IGH,Pancreas.repertoire.diversity$entropy_50_IGH)

 