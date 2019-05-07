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
### DESCRIP: Preparing Immune and Clinical for PAAD 
###         
###
### Author: Silvia Pineda
### Date: April, 2019
############################################################################################
library(ggplot2)
library(Rtsne)
library(glmnet)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(viridis)
setwd("~/TCGA-Immune/")

load("Data/PAAD/PAAD_FullData.Rdata")
PAAD.repertoire.diversity$Tumor_type_2categ<-ifelse(PAAD.repertoire.diversity$Tumor_type=="Tumor_pancreas","Tumor_pancres",
                                             ifelse(PAAD.repertoire.diversity$Tumor_type=="Solid_tissue_normal","Adjacent_normal_pancreas",
                                             ifelse(PAAD.repertoire.diversity$Tumor_type=="Adjacent_normal_pancreas","Adjacent_normal_pancreas",NA)))
PAAD.repertoire.diversity$Tumor_type_2categ<-as.factor(PAAD.repertoire.diversity$Tumor_type_2categ)
##################
####Descriptive analysis to see if there are differences by tumor and adjacent_normal
##################
###T markers
PAAD.repertoire.diversity_treads<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$T_Reads>100),]
T_markers<-c("TRA_expression","TRB_expression","TRD_expression","TRG_expression","Alpha_Beta_ratio_expression",
             "clones_recon_TRA","clones_recon_TRB","clones_recon_TRD","clones_recon_TRG","entropy_recon_TRA",
             "entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG")
p.T_markers=NULL
for(i in 1:length(T_markers)){
  p.T_markers[i]<-coef(summary(glm(PAAD.repertoire.diversity_treads[,T_markers[i]]~PAAD.repertoire.diversity_treads$Tumor_type_2categ)))[2,4]
  PAAD.repertoire.diversity_treads$Marker<-PAAD.repertoire.diversity_treads[,T_markers[i]]
  tiff(paste0("Results/boxplot_",T_markers[i],".tiff"),res=300,h=2000,w=2000)
  print(ggplot(PAAD.repertoire.diversity_treads[which(is.na(PAAD.repertoire.diversity_treads$Tumor_type_2categ)==F),], 
         aes(y=Marker, fill=Tumor_type_2categ, x=Tumor_type_2categ)) + geom_boxplot() + scale_y_continuous(name=T_markers[i])
        + stat_compare_means())
  dev.off()
  
}
T_markers[which(p.T_markers<0.05)] ## "TRD_expression" "TRG_expression"

##Ig markers
PAAD.repertoire.diversity_Igreads<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$IG_Reads>100),]
Ig_markers<-c("IGH_expression","IGK_expression","IGL_expression","KappaLambda_ratio_expression",
             "clones_recon_IGH","clones_recon_IGK","clones_recon_IGL","entropy_recon_IGH",
             "entropy_recon_IGK","entropy_recon_IGL")
p.Ig_markers=NULL
for(i in 1:length(Ig_markers)){
  p.Ig_markers[i]<-coef(summary(glm(PAAD.repertoire.diversity_treads[,Ig_markers[i]]~PAAD.repertoire.diversity_treads$Tumor_type_2categ)))[2,4]
  PAAD.repertoire.diversity_treads$Marker<-PAAD.repertoire.diversity_treads[,Ig_markers[i]]
  tiff(paste0("Results/boxplot_",Ig_markers[i],".tiff"),res=300,h=2000,w=2000)
  print(ggplot(PAAD.repertoire.diversity_treads[which(is.na(PAAD.repertoire.diversity_treads$Tumor_type_2categ)==F),], 
               aes(y=Marker, fill=Tumor_type_2categ, x=Tumor_type_2categ)) + geom_boxplot() + scale_y_continuous(name=Ig_markers[i])
        + stat_compare_means())
  dev.off()
}
Ig_markers[which(p.Ig_markers<0.05)] ## 0

###################################
##### Only tumors #################
###################################
PAAD.repertoire.tumor<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type=="Tumor_pancreas"),] #14
PAAD.repertoire.tumor$TCGA_sample<-substr(PAAD.repertoire.tumor$TCGA_sample,1,15)
xCell.data.tumor<-t(xCell.data.PAAD[,match(PAAD.repertoire.tumor$TCGA_sample,colnames(xCell.data.PAAD))])
xCell.pvalue.tumor<-t(xCell.pvalue.PAAD[,match(PAAD.repertoire.tumor$TCGA_sample,colnames(xCell.pvalue.PAAD))])

####Filter by p-value
xCell.pvalue.tumor.filter<-t(apply(xCell.pvalue.tumor,1,function (x) replace(x,x>=0.2,1)))
xcell.data.tumor.filter <-  xCell.data.tumor[,colSums(xCell.pvalue.tumor.filter==1) <= 148*0.8] ##(148*0.8)  45 cells

#Delete the scores
mat<-t(xcell.data.tumor.filter[ , -which(colnames(xcell.data.tumor.filter) %in% c("ImmuneScore","StromaScore","MicroenvironmentScore"))]) ##Delete the scores
tiff("Results/xCell_CompareSamples.tiff",res=300,w=3500,h=3000)
pheatmap(t(xcell.data.tumor.filter),scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(9,name="PuOr"))(12))
dev.off()

###Applied ENET to find cells associated with T or B expression
##1. Filter those cells that have more than 80% of 0's (148*0.2)
##############
##T_reads
##############
PAAD.repertoire.diversity_treads<-PAAD.repertoire.tumor[which(PAAD.repertoire.tumor$T_Reads>100),] #139
xcell.data.tumor.filter_treads<-xcell.data.tumor.filter[which(PAAD.repertoire.tumor$T_Reads>100),] #139
xcell.data.tumor.filter_treads_mat<-xcell.data.tumor.filter_treads[ , -which(colnames(xcell.data.tumor.filter_treads) %in% c("ImmuneScore","StromaScore","MicroenvironmentScore"))]

alphalist<-seq(0.01,0.99,by=0.01)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(xcell.data.tumor.filter_treads_mat,PAAD.repertoire.diversity_treads$TRA_expression,family="gaussian"
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

enet<-glmnet(xcell.data.tumor.filter_treads_mat,PAAD.repertoire.diversity_treads$TRA_expression,family="gaussian",standardize=TRUE,alpha=alpha,lambda=lambda)
cells<-rownames(enet$beta)[which(enet$beta!=0)]
coef<-enet$beta[which(enet$beta!=0)]

significant_cells<-xcell.data.tumor.filter_treads_mat[,match(cells,colnames(xcell.data.tumor.filter_treads_mat))] #12

## Heatmap ####
annotation_col = data.frame(
  T_expression = PAAD.repertoire.diversity_treads$TRA_expression)

rownames(annotation_col)<-rownames(significant_cells)
ann_colors = list (T_expression = brewer.pal(9,"Reds"))
tiff("Results/T_expression.tiff",res=300,w=3500,h=3000)
pheatmap(t(significant_cells),scale="row",annotation_col = annotation_col,annotation_colors = ann_colors,show_colnames = F,border_color=F,
         color = colorRampPalette(brewer.pal(9,name="PuOr"))(12))
dev.off()

##############
##IG_reads
#############
PAAD.repertoire.diversity_Igreads<-PAAD.repertoire.tumor[which(PAAD.repertoire.tumor$IGH_Reads>100),] #147
xcell.data.tumor.filter_Igreads<-xcell.data.tumor.filter[which(PAAD.repertoire.tumor$IGH_Reads>100),] #147
xcell.data.tumor.filter_Igreads_mat<-xcell.data.tumor.filter_Igreads[ , -which(colnames(xcell.data.tumor.filter_Igreads) %in% c("ImmuneScore","StromaScore","MicroenvironmentScore"))]

alphalist<-seq(0.01,0.99,by=0.01)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(xcell.data.tumor.filter_Igreads_mat,PAAD.repertoire.diversity_Igreads$IG_expression,family="gaussian"
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

enet<-glmnet(xcell.data.tumor.filter_Igreads_mat,PAAD.repertoire.diversity_Igreads$IG_expression,family="gaussian",standardize=TRUE,alpha=alpha,lambda=lambda)
cells<-rownames(enet$beta)[which(enet$beta!=0)]
coef<-enet$beta[which(enet$beta!=0)]

significant_cells<-xcell.data.tumor.filter_Igreads_mat[,match(cells,colnames(xcell.data.tumor.filter_Igreads_mat))] #15

## Heatmap ####
annotation_col = data.frame(
  IG_expression = PAAD.repertoire.diversity_Igreads$IG_expression)

rownames(annotation_col)<-rownames(significant_cells)
ann_colors = list (IG_expression = brewer.pal(6,"Greens"))
tiff("Results/IG_expression.tiff",res=300,w=3500,h=3000)
pheatmap(t(significant_cells),scale="row",annotation_col = annotation_col,annotation_colors = ann_colors,show_colnames = F,border_color=F,
         color = colorRampPalette(brewer.pal(9,name="PuOr"))(12))
dev.off()


###############################
## Merge with Clinical data ###
###############################
clinical.patient.tumor<-clinical.patient[match(substr(PAAD.repertoire.tumor$TCGA_sample,1,12),clinical.patient$bcr_patient_barcode),]
biospecimen.slide.tumor<-biospecimen.slide[match(substr(PAAD.repertoire.tumor$TCGA_sample,1,12),biospecimen.slide$bcr_patient_barcode),]

markers<-c("IG_expression","IGH_expression","IGK_expression","IGL_expression","T_expression","TRA_expression","TRB_expression","TRD_expression",
  "TRG_expression","Alpha_Beta_ratio_expression", "KappaLambda_ratio_expression","clones_recon_IGH","clones_recon_IGK","clones_recon_IGL",
  "clones_recon_TRA","clones_recon_TRB","clones_recon_TRD","clones_recon_TRG","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL",           
  "entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG")

##Histological type
clinical.patient.tumor$histological_type_2cat<-factor(ifelse(clinical.patient.tumor$histological_type=="Pancreas-Adenocarcinoma Ductal Type","PDAC","Other"))
p.markers<-NULL
for(i in 1:length(markers)){
  p.markers[i]<-coef(summary(glm(PAAD.repertoire.tumor[,markers[i]]~clinical.patient.tumor$histological_type_2cat)))[2,4]
}
dfplot <- data.frame(marker=PAAD.repertoire.tumor[,markers[which(p.markers<0.05)]],
            histological_type_2cat=clinical.patient.tumor$histological_type_2cat)
tiff(paste0("Results/boxplot_histological_type_2cat_",markers[which(p.markers<0.05)],".tiff"),res=300,h=2000,w=2000)
ggplot(dfplot,aes(y=marker, fill=histological_type_2cat, x=histological_type_2cat)) + geom_boxplot() + scale_y_continuous(name= markers[which(p.markers<0.05)]) + stat_compare_means()
dev.off()

##anatomic_neoplasm_subdivision
p.markers<-NULL
for(i in 1:length(markers)){
  p.markers[i]<-kruskal.test(PAAD.repertoire.tumor[,markers[i]]~clinical.patient.tumor$anatomic_neoplasm_subdivision)$p.value
}
dfplot <- data.frame(PAAD.repertoire.tumor[,markers[which(p.markers<0.05)]],
                     anatomic_neoplasm_subdivision=clinical.patient.tumor$anatomic_neoplasm_subdivision)

for(i in 1:length(markers[which(p.markers<0.05)])){
  print(i)
  dfplot$marker<-PAAD.repertoire.tumor[,markers[which(p.markers<0.05)][i]]
  tiff(paste0("Results/boxplot_anatomic_neoplasm_subdivision_",markers[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
   print(ggplot(dfplot,aes(y=marker, fill=anatomic_neoplasm_subdivision, x=anatomic_neoplasm_subdivision)) + geom_boxplot() 
        + scale_y_continuous(name= markers[which(p.markers<0.05)][i]) + stat_compare_means())
  dev.off()
}

##gender
p.markers<-NULL
for(i in 1:length(markers)){
  p.markers[i]<-kruskal.test(PAAD.repertoire.tumor[,markers[i]]~clinical.patient.tumor$gender)$p.value
}
dfplot <- data.frame(PAAD.repertoire.tumor[,markers[which(p.markers<0.05)]],
                     gender=clinical.patient.tumor$gender)

for(i in 1:length(markers[which(p.markers<0.05)])){
  print(i)
  dfplot$marker<-PAAD.repertoire.tumor[,markers[which(p.markers<0.05)][i]]
  tiff(paste0("Results/boxplot_gender_",markers[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
  print(ggplot(dfplot,aes(y=marker, fill=gender, x=gender)) + geom_boxplot() 
        + scale_y_continuous(name= markers[which(p.markers<0.05)][i]) + stat_compare_means())
  dev.off()
}

##race_list
p.markers<-NULL
for(i in 1:length(markers)){
  p.markers[i]<-kruskal.test(PAAD.repertoire.tumor[,markers[i]]~clinical.patient.tumor$race_list)$p.value
}
dfplot <- data.frame(marker =PAAD.repertoire.tumor[,markers[which(p.markers<0.05)]],
                     race_list=clinical.patient.tumor$race_list)

for(i in 1:length(markers[which(p.markers<0.05)])){
  print(i)
   tiff(paste0("Results/boxplot_race_",markers[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
  print(ggplot(dfplot,aes(y=marker, fill=race_list, x=race_list)) + geom_boxplot() 
        + scale_y_continuous(name= markers[which(p.markers<0.05)][i]) + stat_compare_means())
  dev.off()
}

##lymph_node_examined_count
p.markers<-NULL
for(i in 1:length(markers)){
  p.markers[i]<-coef(summary(glm(PAAD.repertoire.tumor[,markers[i]]~clinical.patient.tumor$lymph_node_examined_count)))[2,4]
}
dfplot <- data.frame(marker =PAAD.repertoire.tumor[,markers[which(p.markers<0.05)]],
                     lymph_node=clinical.patient.tumor$lymph_node_examined_count)

for(i in 1:length(markers[which(p.markers<0.05)])){
  print(i)
  tiff(paste0("Results/boxplot_lymph_node_",markers[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
  print(ggplot(dfplot,aes(y=marker, fill=lymph_node, x=lymph_node)) + geom_point() + geom_smooth(method='lm')
        + scale_y_continuous(name= markers[which(p.markers<0.05)][i]) + stat_cor(method = "pearson"))
  dev.off()
}

##neoplasm_histologic_grade
clinical.patient.tumor$neoplasm_histologic_grade_3cat<-factor(ifelse(clinical.patient.tumor$neoplasm_histologic_grade=="G1","G1",
                                                                     ifelse(clinical.patient.tumor$neoplasm_histologic_grade=="G2","G2",
                                                                            ifelse(clinical.patient.tumor$neoplasm_histologic_grade=="G3","G3",""))))

##vital_status
p.markers<-NULL
for(i in 1:length(markers)){
  p.markers[i]<-kruskal.test(PAAD.repertoire.tumor[,markers[i]]~clinical.patient.tumor$vital_status)$p.value
}
dfplot <- data.frame(PAAD.repertoire.tumor[,markers[which(p.markers<0.05)]],
                     vital_status=clinical.patient.tumor$vital_status)

for(i in 1:length(markers[which(p.markers<0.05)])){
  print(i)
  dfplot$marker<-PAAD.repertoire.tumor[,markers[which(p.markers<0.05)][i]]
  tiff(paste0("Results/boxplot_vital_status_",markers[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
  print(ggplot(dfplot,aes(y=marker, fill=vital_status, x=vital_status)) + geom_boxplot() 
        + scale_y_continuous(name= markers[which(p.markers<0.05)][i]) + stat_compare_means())
  dev.off()
}

##Smoking
clinical.patient.tumor$smoking<-factor(ifelse(clinical.patient.tumor$tobacco_smoking_history_master=="Current smoker (includes daily smokers and non-daily smokers or occasional smokers)","Current",
                                       ifelse(clinical.patient.tumor$tobacco_smoking_history_master=="Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime)","Non-smoker","Former")))

##number_pack_years_smoked
p.markers<-NULL
for(i in 1:length(markers)){
  p.markers[i]<-coef(summary(glm(PAAD.repertoire.tumor[,markers[i]]~clinical.patient.tumor$number_pack_years_smoked)))[2,4]
}
dfplot <- data.frame(marker =PAAD.repertoire.tumor[,markers[which(p.markers<0.05)]],
                     number_pack_years_smoked=clinical.patient.tumor$number_pack_years_smoked)

for(i in 1:length(markers[which(p.markers<0.05)])){
  print(i)
  dfplot$marker<-PAAD.repertoire.tumor[,markers[which(p.markers<0.05)][i]]
  tiff(paste0("Results/boxplot_number_pack_years_smoked_",markers[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
  print(ggplot(dfplot,aes(y=marker, fill=number_pack_years_smoked, x=number_pack_years_smoked)) + geom_point() + geom_smooth(method='lm')
        + scale_y_continuous(name= markers[which(p.markers<0.05)][i]) + stat_cor(method = "pearson"))
  dev.off()
}


##history_diabetes
clinical.patient.tumor$history_of_diabetes<-factor(clinical.patient.tumor$history_of_diabetes)
p.markers<-NULL
for(i in 1:length(markers)){
  p.markers[i]<-kruskal.test(PAAD.repertoire.tumor[which(clinical.patient.tumor$history_of_diabetes!=""),markers[i]]~clinical.patient.tumor$history_of_diabetes[which(clinical.patient.tumor$history_of_diabetes!="")])$p.value
}
dfplot <- data.frame(PAAD.repertoire.tumor[which(clinical.patient.tumor$history_of_diabetes!=""),markers[which(p.markers<0.05)]],
                     history_of_diabetes=clinical.patient.tumor$history_of_diabetes[which(clinical.patient.tumor$history_of_diabetes!="")])

for(i in 1:length(markers[which(p.markers<0.05)])){
  print(i)
  dfplot$marker<-PAAD.repertoire.tumor[which(clinical.patient.tumor$history_of_diabetes!=""),markers[which(p.markers<0.05)][i]]
  tiff(paste0("Results/boxplot_history_of_diabetes_",markers[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
  print(ggplot(dfplot,aes(y=marker, fill=history_of_diabetes, x=history_of_diabetes)) + geom_boxplot() 
        + scale_y_continuous(name= markers[which(p.markers<0.05)][i]) + stat_compare_means())
  dev.off()
}

##history_chronic_pancreatitis
clinical.patient.tumor$history_of_chronic_pancreatitis<-factor(clinical.patient.tumor$history_of_chronic_pancreatitis)
p.markers<-NULL
for(i in 1:length(markers)){
  p.markers[i]<-kruskal.test(PAAD.repertoire.tumor[which(clinical.patient.tumor$history_of_chronic_pancreatitis!=""),markers[i]]~clinical.patient.tumor$history_of_chronic_pancreatitis
                             [which(clinical.patient.tumor$history_of_chronic_pancreatitis!="")])$p.value
}
dfplot <- data.frame(PAAD.repertoire.tumor[which(clinical.patient.tumor$history_of_chronic_pancreatitis!=""),markers[which(p.markers<0.05)]],
                     history_of_chronic_pancreatitis=clinical.patient.tumor$history_of_chronic_pancreatitis[which(clinical.patient.tumor$history_of_chronic_pancreatitis!="")])

for(i in 1:length(markers[which(p.markers<0.05)])){
  print(i)
  dfplot$marker<-PAAD.repertoire.tumor[which(clinical.patient.tumor$history_of_chronic_pancreatitis!=""),markers[which(p.markers<0.05)][i]]
  tiff(paste0("Results/boxplot_history_of_chronic_pancreatitis_",markers[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
  print(ggplot(dfplot,aes(y=marker, fill=history_of_chronic_pancreatitis, x=history_of_chronic_pancreatitis)) + geom_boxplot() 
        + scale_y_continuous(name= markers[which(p.markers<0.05)][i]) + stat_compare_means())
  dev.off()
}

##Gender with xcell
clinical.patient.xcell<-clinical.patient[match(substr(rownames(xcell.data.tumor.filter),1,12),clinical.patient$bcr_patient_barcode),]

p.markers<-NULL
for(i in 1:ncol(xcell.data.tumor.filter)){
  p.markers[i]<-kruskal.test(xcell.data.tumor.filter[,i]~clinical.patient.xcell$gender)$p.value
}
dfplot <- data.frame(xcell.data.tumor.filter[,which(p.markers<0.05)],
                     gender=clinical.patient.xcell$gender)

for(i in 1:length(which(p.markers<0.05))){
  print(i)
  dfplot$marker<-xcell.data.tumor.filter[,which(p.markers<0.05)[i]]
  tiff(paste0("Results/boxplot_gender_",colnames(xcell.data.tumor.filter)[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
  print(ggplot(dfplot,aes(y=marker, fill=gender, x=gender)) + geom_boxplot() 
        + scale_y_continuous(name= colnames(xcell.data.tumor.filter)[which(p.markers<0.05)][i]) + stat_compare_means())
  dev.off()
}
