


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
### DESCRIP: Step 2: Association analysis with SNps genotypes
###         
###
### Author: Silvia Pineda
### Date: January, 2019
############################################################################################
library(ggplot2)
library("dplyr")
library("RColorBrewer")
library(stringr)
library(glmnet)

setwd("~/TCGA-Immune/")
load("Data/PAAD/PAAD_FullData.Rdata")

PAAD_genotypes<-read.csv("Data/PAAD/callSNPsPAAD.csv")
rownames(PAAD_genotypes)<-PAAD_genotypes$SNP_name
PAAD_genotypes<-PAAD_genotypes[,-1]
TCGA_names<-substr(colnames(PAAD_genotypes),1,12)
TCGA_names<-str_replace(TCGA_names,"\\.","-")
TCGA_names<-str_replace(TCGA_names,"\\.","-")
colnames(PAAD_genotypes)<-TCGA_names

PAAD.repertoire.diversity$TCGA_sample<-substr(PAAD.repertoire.diversity$TCGA_sample,1,12)
id<-match(colnames(PAAD_genotypes),PAAD.repertoire.diversity$TCGA_sample)
PAAD.repertoire.diversity.genotypes<-PAAD.repertoire.diversity[na.omit(id),]
id<-match(PAAD.repertoire.diversity.genotypes$TCGA_sample,colnames(PAAD_genotypes))
PAAD.genotypes<-t(PAAD_genotypes[,id])


###### ENET 
alphalist<-seq(0.1,0.9,by=0.01)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(PAAD.genotypes,log10(PAAD.repertoire.diversity.genotypes$IG_expression),family="gaussian"
                                                        ,standardize=TRUE,alpha=a,nfolds=5))})
xx<-rep(NA,length(alphalist))
yy<-rep(NA,length(alphalist))
for (j in 1:length(alphalist)) {
  #print(j)
  if(class(elasticnet[[j]]) != "try-error"){
    xx[j]<-elasticnet[[j]]$lambda.min
    id.cv.opt<-grep(elasticnet[[j]]$lambda.min,elasticnet[[j]]$lambda,fixed=TRUE)
    yy[j]<-elasticnet[[j]]$cvm[id.cv.opt]
  }
}
id.min<-which(yy==min(yy,na.rm=TRUE))
lambda<-xx[id.min]
alpha<-alphalist[id.min]

enet<-glmnet(PAAD.genotypes,log10(PAAD.repertoire.diversity.genotypes$IG_expression),family="gaussian",standardize=TRUE,alpha=alpha,lambda=lambda)
snps<-rownames(enet$beta)[which(enet$beta!=0)]
coef<-enet$beta[which(enet$beta!=0)]

significant_snps<-PAAD_genotypes[match(snps,rownames(PAAD_genotypes)),] #16

#"SNP_A-1813059" "SNP_A-1825881" "SNP_A-1978600" "SNP_A-2010608" "SNP_A-2085054" "SNP_A-2181368" "SNP_A-2285564" "SNP_A-2309881"
#"SNP_A-2310951" "SNP_A-4292944" "SNP_A-8344170" "SNP_A-8378841" "SNP_A-8423745" "SNP_A-8433670" "SNP_A-8442577" "SNP_A-8448026"

snp_annotation<-read.csv("Data/PAAD/GenomeWideSNP_6.na35.annot.csv",sep="//")


library(caret)
set.seed(3456)
trainIndex <- createDataPartition(log10(PAAD.repertoire.diversity.genotypes$IG_expression), p = .8, 
                                  list = FALSE, 
                                  times = 1)
head(trainIndex)