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

##################################
## Heatmap for the BCR and TCR ###
#################################
#IG
PAAD.repertoire.diversity_Igreads<-PAAD.repertoire.tumor[which(PAAD.repertoire.tumor$Ig_Reads>100),]
Ig_markers<-c("IGH_expression","IGK_expression","IGL_expression")
mat<-PAAD.repertoire.diversity_Igreads[,Ig_markers]
rownames(mat)<-substr(PAAD.repertoire.diversity_Igreads$TCGA_sample,1,12)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))

#TR
PAAD.repertoire.diversity_treads<-PAAD.repertoire.tumor[which(PAAD.repertoire.tumor$T_Reads>100),]
T_markers<-c("TRA_expression","TRB_expression","TRD_expression","TRG_expression")
mat<-PAAD.repertoire.diversity_treads[,T_markers]
rownames(mat)<-substr(PAAD.repertoire.diversity_treads$TCGA_sample,1,12)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))

#All
ALL_markers<-c("IGH_expression","IGK_expression","IGL_expression","TRA_expression","TRB_expression","TRD_expression","TRG_expression")
mat<-PAAD.repertoire.tumor[,ALL_markers]
rownames(mat)<-substr(PAAD.repertoire.tumor$TCGA_sample,1,12)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))

###Applied ENET to find cells associated with T or B expression
##############
##T_reads
##############
PAAD.repertoire.diversity_treads<-PAAD.repertoire.tumor[which(PAAD.repertoire.tumor$T_Reads>100),] #139
xcell.data.tumor.filter_treads<-xcell.data.tumor.filter[which(PAAD.repertoire.tumor$T_Reads>100),] #139
xcell.data.tumor.filter_treads_mat<-xcell.data.tumor.filter_treads[ , -which(colnames(xcell.data.tumor.filter_treads) %in% c("ImmuneScore","StromaScore","MicroenvironmentScore"))]

library(dplyr)
PAAD.repertoire.diversity_treads$T_expression_tertiles<-ntile(PAAD.repertoire.diversity_treads$T_expression,3)

alphalist<-seq(0.1,0.9,by=0.1)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(xcell.data.tumor.filter_treads_mat,PAAD.repertoire.diversity_treads$T_expression,family="gaussian"
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

enet<-glmnet(xcell.data.tumor.filter_treads_mat,PAAD.repertoire.diversity_treads$T_expression,family="gaussian",standardize=TRUE,alpha=alpha,lambda=lambda)
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
PAAD.repertoire.diversity_Igreads<-PAAD.repertoire.tumor[which(PAAD.repertoire.tumor$Ig_Reads>100),] #147
xcell.data.tumor.filter_Igreads<-xcell.data.tumor.filter[which(PAAD.repertoire.tumor$Ig_Reads>100),] #147
xcell.data.tumor.filter_Igreads_mat<-xcell.data.tumor.filter_Igreads[ , -which(colnames(xcell.data.tumor.filter_Igreads) %in% c("ImmuneScore","StromaScore","MicroenvironmentScore"))]

alphalist<-seq(0.1,0.9,by=0.01)
set.seed(54)
elasticnet<-lapply(alphalist, function(a){try(cv.glmnet(xcell.data.tumor.filter_Igreads_mat,log10(PAAD.repertoire.diversity_Igreads$IG_expression),family="gaussian"
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

enet<-glmnet(xcell.data.tumor.filter_Igreads_mat,log10(PAAD.repertoire.diversity_Igreads$IG_expression),family="gaussian",standardize=TRUE,alpha=alpha,lambda=lambda)
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
PAAD.repertoire.tumor<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type_2categ=="Tumor_pancreas"),]
clinical.patient.tumor<-clinical.patient[match(substr(PAAD.repertoire.tumor$TCGA_sample,1,12),clinical.patient$bcr_patient_barcode),]
PAAD.repertoire.tumor.clinical.patient<-cbind(PAAD.repertoire.tumor,clinical.patient.tumor)

##Function to run the association between clinical outcome and BCR/TCR
association.test.immuneRep<- function (PAAD.repertoire.tumor.clinical.patient,clinical.var){
  Ig_expr<-melt(PAAD.repertoire.tumor.clinical.patient[,c("TCGA_sample","IGH_expression","IGK_expression","IGL_expression",clinical.var)])
  Ig_expr$value<-log10(Ig_expr$value)
  Ig_expr<-na.omit(Ig_expr)
  tiff(paste0("Results/boxplot_Ig_expression_TCGA_",clinical.var,".tiff"),res=300,h=2500,w=3500)
  if(length(table(Ig_expr[,clinical.var]))<10){  
  print(ggboxplot(Ig_expr, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
    rotate_x_text() +
    geom_point(aes(x=Ig_expr[,clinical.var], y=value,color=Ig_expr[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
    stat_compare_means(label = "p.format"))
  } else {  
    print(ggplot(Ig_expr,aes(x = as.numeric(Ig_expr[,clinical.var]), y = value)) + geom_point() + 
      facet_grid(.~Ig_expr$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()

  Ig_entropy<-melt(PAAD.repertoire.tumor.clinical.patient[,c("TCGA_sample","entropy_recon_IGH","entropy_recon_IGK","entropy_recon_IGL",clinical.var)])
  Ig_entropy<-Ig_entropy[which(Ig_entropy$value!=0),]
  Ig_entropy<-na.omit(Ig_entropy)
  tiff(paste0("Results/boxplot_Ig_entropy_TCGA_",clinical.var,".tiff"),res=300,h=2500,w=3500)
  if(length(table(Ig_entropy[,clinical.var]))<10){  
    print(ggboxplot(Ig_entropy, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
    rotate_x_text() +
    geom_point(aes(x=Ig_entropy[,clinical.var], y=value,color=Ig_entropy[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
    stat_compare_means(label = "p.format"))
  } else {  
    print(ggplot(Ig_entropy,aes(x = as.numeric(Ig_entropy[,clinical.var]), y = value)) + geom_point() + 
          facet_grid(.~Ig_entropy$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()

  T_expr<-melt(PAAD.repertoire.tumor.clinical.patient[,c("TCGA_sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression",clinical.var)])
  T_expr$value<-log10(T_expr$value)
  T_expr<-na.omit(T_expr)
  tiff(paste0("Results/boxplot_T_expr_TCGA_",clinical.var,".tiff"),res=300,h=2500,w=3500)
  if(length(table(T_expr[,clinical.var]))<10){  
    print(ggboxplot(T_expr, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
    rotate_x_text() +
    geom_point(aes(x=T_expr[,clinical.var], y=value,color=T_expr[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
    stat_compare_means(label = "p.format"))
  } else {  
    print(ggplot(T_expr,aes(x = as.numeric(T_expr[,clinical.var]), y = value)) + geom_point() + 
            facet_grid(.~T_expr$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()

  T_entropy<-melt(PAAD.repertoire.tumor.clinical.patient[,c("TCGA_sample","entropy_recon_TRA","entropy_recon_TRB","entropy_recon_TRD","entropy_recon_TRG",clinical.var)])
  T_entropy<-T_entropy[which(T_entropy$value!=0),]
  T_entropy<-na.omit(T_entropy)
  tiff(paste0("Results/boxplot_T_entropy_TCGA_",clinical.var,".tiff"),res=300,h=2500,w=3500)
  if(length(table(T_entropy[,clinical.var]))<10){  
    print(ggboxplot(T_entropy, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
    rotate_x_text() +
    geom_point(aes(x=T_entropy[,clinical.var], y=value,color=T_entropy[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
    stat_compare_means(label = "p.format"))
  } else {  
    print(ggplot(T_entropy,aes(x = as.numeric(T_entropy[,clinical.var]), y = value)) + geom_point() + 
            facet_grid(.~T_entropy$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()
  

  Alpha_Beta_ratio_expression<-melt(PAAD.repertoire.tumor.clinical.patient[,c("TCGA_sample","Alpha_Beta_ratio_expression",clinical.var)])
  Alpha_Beta_ratio_expression<-na.omit(Alpha_Beta_ratio_expression)
  tiff(paste0("Results/boxplot_Alpha_Beta_ratio_expression_TCGA_",clinical.var,".tiff"),res=300,h=2500,w=3500)
  if(length(table(Alpha_Beta_ratio_expression[,clinical.var]))<10){  
    print(ggboxplot(Alpha_Beta_ratio_expression, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
    rotate_x_text() +
    geom_point(aes(x=Alpha_Beta_ratio_expression[,clinical.var], y=value,color=Alpha_Beta_ratio_expression[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
    stat_compare_means(label = "p.format"))
  } else {  
    print(ggplot(Alpha_Beta_ratio_expression,aes(x = as.numeric(Alpha_Beta_ratio_expression[,clinical.var]), y = value)) + geom_point() + 
            facet_grid(.~Alpha_Beta_ratio_expression$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson") + xlab(clinical.var))
  }
  dev.off()

  KappaLambda_ratio_expression<-melt(PAAD.repertoire.tumor.clinical.patient[,c("TCGA_sample","KappaLambda_ratio_expression",clinical.var)])
  KappaLambda_ratio_expression<-na.omit(KappaLambda_ratio_expression)
  tiff(paste0("Results/boxplot_KappaLambda_ratio_expression_TCGA_",clinical.var,".tiff"),res=300,h=2500,w=3500)
  if(length(table(KappaLambda_ratio_expression[,clinical.var]))<10){  
    print(ggboxplot(KappaLambda_ratio_expression, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
    rotate_x_text() +
    geom_point(aes(x=KappaLambda_ratio_expression[,clinical.var], y=value,color=KappaLambda_ratio_expression[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
    stat_compare_means(label = "p.format"))
  } else {  
    print(ggplot(KappaLambda_ratio_expression,aes(x = as.numeric(KappaLambda_ratio_expression[,clinical.var]), y = value)) + geom_point() + 
            facet_grid(.~KappaLambda_ratio_expression$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()
}

##Histological type
PAAD.repertoire.tumor.clinical.patient$histological_type_2cat<-factor(ifelse(PAAD.repertoire.tumor.clinical.patient$histological_type=="Pancreas-Adenocarcinoma Ductal Type","PDAC","PC-Other"))
PAAD.repertoire.tumor.clinical.patient$histological_type_2cat<-factor(PAAD.repertoire.tumor.clinical.patient$histological_type_2cat)
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"histological_type_2cat")

##anatomic_neoplasm_subdivision
PAAD.repertoire.tumor.clinical.patient$anatomic_neoplasm_subdivision<-factor(PAAD.repertoire.tumor.clinical.patient$anatomic_neoplasm_subdivision)
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"anatomic_neoplasm_subdivision")

##gender
PAAD.repertoire.tumor.clinical.patient$gender<-factor(PAAD.repertoire.tumor.clinical.patient$gender)
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"gender")

##race_list
PAAD.repertoire.tumor.clinical.patient$race_list<-ifelse(PAAD.repertoire.tumor.clinical.patient$race_list=="",NA,as.character(PAAD.repertoire.tumor.clinical.patient$race_list))
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"race_list")

##History of Prior Malignancy
PAAD.repertoire.tumor.clinical.patient$other_dx<-factor(PAAD.repertoire.tumor.clinical.patient$other_dx)
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"other_dx")


##neoplasm_histologic_grade
PAAD.repertoire.tumor.clinical.patient$neoplasm_histologic_grade_3cat<-factor(ifelse(PAAD.repertoire.tumor.clinical.patient$neoplasm_histologic_grade=="G1","G1",
                                                                     ifelse(PAAD.repertoire.tumor.clinical.patient$neoplasm_histologic_grade=="G2","G2",
                                                                            ifelse(PAAD.repertoire.tumor.clinical.patient$neoplasm_histologic_grade=="G3","G3",NA))))
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"neoplasm_histologic_grade_3cat")

##Age 
PAAD.repertoire.tumor.clinical.patient$age_at_initial_pathologic_diagnosis<-factor(PAAD.repertoire.tumor.clinical.patient$age_at_initial_pathologic_diagnosis)
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"age_at_initial_pathologic_diagnosis")


##Smoking
PAAD.repertoire.tumor.clinical.patient$smoking<-factor(ifelse(PAAD.repertoire.tumor.clinical.patient$tobacco_smoking_history_master=="Current smoker (includes daily smokers and non-daily smokers or occasional smokers)","Current",
                                       ifelse(PAAD.repertoire.tumor.clinical.patient$tobacco_smoking_history_master=="Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime)","Non-smoker","Former")))
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"smoking")

##number_pack_years_smoked
PAAD.repertoire.tumor.clinical.patient$number_pack_years_smoked<-factor(PAAD.repertoire.tumor.clinical.patient$number_pack_years_smoked)
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"number_pack_years_smoked")

##Alcohol
PAAD.repertoire.tumor.clinical.patient$alcohol_history_documented<-ifelse(PAAD.repertoire.tumor.clinical.patient$alcohol_history_documented=="",NA,
                                                                          as.character(PAAD.repertoire.tumor.clinical.patient$alcohol_history_documented))
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"alcohol_history_documented")

##Alcohol category
PAAD.repertoire.tumor.clinical.patient$alcoholic_exposure_category2<-ifelse(PAAD.repertoire.tumor.clinical.patient$alcohol_history_documented=="NO","No-drinker",
                                                     ifelse(PAAD.repertoire.tumor.clinical.patient$alcohol_history_documented=="YES" & PAAD.repertoire.tumor.clinical.patient$alcoholic_exposure_category=="",NA,
                                                     ifelse(PAAD.repertoire.tumor.clinical.patient$alcoholic_exposure_category=="None","None-Drinker",
                                                     ifelse(PAAD.repertoire.tumor.clinical.patient$alcoholic_exposure_category=="Occasional Drinker","Occasional-Drinker",
                                                     ifelse(PAAD.repertoire.tumor.clinical.patient$alcoholic_exposure_category=="Daily Drinker","Daily-Drinker",
                                                     ifelse(PAAD.repertoire.tumor.clinical.patient$alcoholic_exposure_category=="Social Drinker","Social-Drinker",
                                                     ifelse(PAAD.repertoire.tumor.clinical.patient$alcoholic_exposure_category=="Weekly Drinker","Weekly-Drinker",NA)))))))
PAAD.repertoire.tumor.clinical.patient$alcoholic_exposure_category2<-factor(PAAD.repertoire.tumor.clinical.patient$alcoholic_exposure_category2)
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"alcoholic_exposure_category2")

##family history
PAAD.repertoire.tumor.clinical.patient$family_history_of_cancer<-factor(ifelse(PAAD.repertoire.tumor.clinical.patient$family_history_of_cancer=="",NA,PAAD.repertoire.tumor.clinical.patient$family_history_of_cancer))
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"family_history_of_cancer")

##radiation_therapy
PAAD.repertoire.tumor.clinical.patient$radiation_therapy<-factor(ifelse(PAAD.repertoire.tumor.clinical.patient$radiation_therapy=="",NA,
                                                                        PAAD.repertoire.tumor.clinical.patient$radiation_therapy))
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"radiation_therapy")

##primary_therapy_outcome_success
PAAD.repertoire.tumor.clinical.patient$primary_therapy_outcome_success<-ifelse(PAAD.repertoire.tumor.clinical.patient$primary_therapy_outcome_success=="",NA,
                                                                               PAAD.repertoire.tumor.clinical.patient$primary_therapy_outcome_success)
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"primary_therapy_outcome_success")

##history_chronic_pancreatitis
PAAD.repertoire.tumor.clinical.patient$history_of_chronic_pancreatitis<-factor(ifelse(PAAD.repertoire.tumor.clinical.patient$history_of_chronic_pancreatitis=="",NA,
                                                                                      PAAD.repertoire.tumor.clinical.patient$history_of_chronic_pancreatitis))
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"history_of_chronic_pancreatitis")

#stage_event_tnm_categories
PAAD.repertoire.tumor.clinical.patient$pathologic_stage<-factor(ifelse(PAAD.repertoire.tumor.clinical.patient$stage_event_pathologic_stage=="Stage IA" | 
                                                                         PAAD.repertoire.tumor.clinical.patient$stage_event_pathologic_stage=="Stage IB", "Stage I",
                                        ifelse(PAAD.repertoire.tumor.clinical.patient$stage_event_pathologic_stage == "Stage IIA" |
                                                 PAAD.repertoire.tumor.clinical.patient$stage_event_pathologic_stage=="Stage IIB","Stage II",
                                         ifelse(PAAD.repertoire.tumor.clinical.patient$stage_event_pathologic_stage=="Stage III","Stage III",
                                         ifelse(PAAD.repertoire.tumor.clinical.patient$stage_event_pathologic_stage=="Stage IV", "Stage IV",NA)))))
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"pathologic_stage")

##history_diabetes
PAAD.repertoire.tumor.clinical.patient$history_of_diabetes<-ifelse(PAAD.repertoire.tumor.clinical.patient$history_of_diabetes=="",NA,
                                                                   as.character(PAAD.repertoire.tumor.clinical.patient$history_of_diabetes))
association.test.immuneRep(PAAD.repertoire.tumor.clinical.patient,"history_of_diabetes")



######################
### Random Forest ###
#####################
library("randomForest")
clinical_variablers<-c("histological_type_2cat","anatomic_neoplasm_subdivision","gender","race_list","other_dx",
                       "lymph_node_examined_count","neoplasm_histologic_grade_3cat","age_at_initial_pathologic_diagnosis",
                       "smoking","number_pack_years_smoked","alcohol_history_documented","alcoholic_exposure_category2",
                       "family_history_of_cancer","radiation_therapy","primary_therapy_outcome_success",
                       "history_of_chronic_pancreatitis","history_of_diabetes")

rf_output <- randomForest(PAAD.repertoire.tumor$clones_recon_TRA~.,data=clinical.patient.tumor[,clinical_variablers],
                          importance=T,proximity=TRUE, keep.forest=T,na.action=na.omit)

## Look at variable importance:
round(importance(rf_output), 2)
varImpPlot(rf_output)


##################################
#######Clinical follow-up########
#################################
PAAD.repertoire.tumor<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type_2categ=="Tumor_pancreas"),]
clinical.follow_up.tumor<-clinical.folow_up[match(substr(PAAD.repertoire.tumor$TCGA_sample,1,12),clinical.folow_up$bcr_patient_barcode),]
PAAD.repertoire.tumor.followuop<-cbind(PAAD.repertoire.tumor,clinical.follow_up.tumor)

##vital_status
PAAD.repertoire.tumor.followuop$vital_status<-factor(PAAD.repertoire.tumor.followuop$vital_status)
association.test.immuneRep(PAAD.repertoire.tumor.followuop,"vital_status")

##new tumor event
PAAD.repertoire.tumor.followuop$new_tumor_event_type<-replace(PAAD.repertoire.tumor.followuop$new_tumor_event_type,PAAD.repertoire.tumor.followuop$new_tumor_event_type=="#N/A",NA)
PAAD.repertoire.tumor.followuop$new_tumor_event_type<-replace(PAAD.repertoire.tumor.followuop$new_tumor_event_type,
                                                              PAAD.repertoire.tumor.followuop$new_tumor_event_type=="Locoregional Recurrence|Distant Metastasis" | 
                                                                PAAD.repertoire.tumor.followuop$new_tumor_event_type=="New Primary Tumor",NA)
PAAD.repertoire.tumor.followuop$new_tumor_event_type<-factor(PAAD.repertoire.tumor.followuop$new_tumor_event_type)
association.test.immuneRep(PAAD.repertoire.tumor.followuop,"new_tumor_event_type")

#treatment_outcome_first_course
PAAD.repertoire.tumor.followuop$treatment_outcome_first_course<-replace(PAAD.repertoire.tumor.followuop$treatment_outcome_first_course,
                                                                        PAAD.repertoire.tumor.followuop$treatment_outcome_first_course=="[Discrepancy]" |
                                                                          PAAD.repertoire.tumor.followuop$treatment_outcome_first_course=="[Not Applicable]" | 
                                                                          PAAD.repertoire.tumor.followuop$treatment_outcome_first_course=="[Not Available]" |
                                                                          PAAD.repertoire.tumor.followuop$treatment_outcome_first_course=="[Unknown]",NA)
PAAD.repertoire.tumor.followuop$treatment_outcome_first_course<-factor(PAAD.repertoire.tumor.followuop$treatment_outcome_first_course)
association.test.immuneRep(PAAD.repertoire.tumor.followuop,"treatment_outcome_first_course")

##############################################
### Survival Analysis  ####
##############################################
PAAD.repertoire.tumor.survival<-merge(PAAD.repertoire.tumor.clinical.patient,PAAD.repertoire.tumor.followuop,by="bcr_patient_barcode")
library(survival)
library(survminer)
library(survMisc)
##OS
surv_object <- Surv(time = PAAD.repertoire.tumor.survival$OS.time, event = PAAD.repertoire.tumor.survival$OS)
res.cox <- coxph(surv_object~PAAD.repertoire.tumor.survival$entropy_recon_IGH.x+PAAD.repertoire.tumor.survival$gender+PAAD.repertoire.tumor.survival$race_list
               +PAAD.repertoire.tumor.survival$age_at_initial_pathologic_diagnosis+PAAD.repertoire.tumor.survival$pathologic_stage)
summary(res.cox)
##Categorical
KL_mean<-mean(PAAD.repertoire.tumor$entropy_recon_IGH)
PAAD.repertoire.tumor$KL_ratio_2cat<-ifelse(PAAD.repertoire.tumor$entropy_recon_IGH<=KL_mean,1,2)
fit1 <- survfit(surv_object ~ PAAD.repertoire.tumor$KL_ratio_2cat)
fit1
ggsurvplot(fit1, data = PAAD.repertoire.tumor)
comp(ten(fit1))$tests$lrTests

#DSS
clinical.follow_up.tumor.DSS<-clinical.follow_up.tumor[which(clinical.follow_up.tumor$DSS!="#N/A"),]
clinical.follow_up.tumor.DSS$DSS<-as.integer(as.character(clinical.follow_up.tumor.DSS$DSS))
surv_object <- Surv(time = clinical.follow_up.tumor.DSS$DSS.time, event = clinical.follow_up.tumor.DSS$DSS)
res.cox <- coxph(surv_object~PAAD.repertoire.tumor$entropy_recon_TRB[which(clinical.follow_up.tumor$DSS!="#N/A")])
summary(res.cox)

#Categorical
KL_mean<-mean(PAAD.repertoire.diversity.gini_IGHV2$KappaLambda_ratio_expression)
PAAD.repertoire.diversity.gini_IGHV2$KL_ratio_2cat<-ifelse(PAAD.repertoire.diversity.gini_IGHV2$KappaLambda_ratio_expression<=KL_mean,1,2)
fit1 <- survfit(surv_object ~ PAAD.repertoire.diversity.gini_IGHV2$KL_ratio_2cat)
fit1
ggsurvplot(fit1, data = PAAD.repertoire.diversity.gini_IGHV2)
comp(ten(fit1))$tests$lrTests

##PFI
surv_object <- Surv(time = clinical.follow_up.tumor$PFI.time, event = clinical.follow_up.tumor$PFI)
res.cox <- coxph(surv_object~PAAD.repertoire.tumor$entropy_recon_TRB)
summary(res.cox)



################################
##clinical outcome with xcell###
################################
clinical.patient.xcell<-clinical.patient[match(substr(rownames(xcell.data.tumor.filter),1,12),clinical.patient$bcr_patient_barcode),]
clinical.follow.up.xcell<-clinical.folow_up[match(substr(rownames(xcell.data.tumor.filter),1,12),clinical.folow_up$bcr_patient_barcode),]

##Function to run the association between clinical outcome and xcell
association.test.xcell<- function (clinical.patient.xcell,clinical.var,xcell.data.tumor.filter){
  if (class(clinical.patient.xcell[,clinical.var])=="factor"){
    p.markers<-NULL
    for(i in 1: ncol(xcell.data.tumor.filter)){
      p.markers[i]<-kruskal.test(xcell.data.tumor.filter[which(clinical.patient.xcell[,clinical.var]!=""),i]~
                                   clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var])$p.value
    }
    dfplot <- data.frame(xcell.data.tumor.filter[which(clinical.patient.xcell[,clinical.var]!=""),which(p.markers<0.05)],
                         clinical.var = clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var])
    colnames(dfplot)[ncol(dfplot)]<- clinical.var
    
    if (dim(dfplot)[2]>1){
      for(i in 1:length(which(p.markers<0.05))){
        print(i)
        dfplot$marker<-xcell.data.tumor.filter[which(clinical.patient.xcell[,clinical.var]!=""),which(p.markers<0.05)[i]]
        tiff(paste0("Results/boxplot_",clinical.var,"_",colnames(xcell.data.tumor.filter)[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
        print(ggplot(dfplot,aes(y=marker, fill=clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var], 
                                x=clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var])) + geom_boxplot() 
              + scale_y_continuous(name=colnames(xcell.data.tumor.filter)[which(p.markers<0.05)][i]) + stat_compare_means()  + scale_x_discrete(name=clinical.var)
              + scale_fill_discrete(name=clinical.var))
        dev.off()
      }
    }
  } else {
    p.markers<-NULL
    for(i in 1:ncol(xcell.data.tumor.filter)){
      p.markers[i]<-coef(summary(glm(xcell.data.tumor.filter[which(clinical.patient.xcell[,clinical.var]!=""),i]~
                                       clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var])))[2,4]
    }
    dfplot <- data.frame(marker =xcell.data.tumor.filter[which(clinical.patient.xcell[,clinical.var]!=""),which(p.markers<0.05)],
                         clinical.var = clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var])
    colnames(dfplot)[ncol(dfplot)]<- clinical.var
    
    if (dim(dfplot)[2]>1){
      for(i in 1:length(which(p.markers<0.05))){
        print(i)
        dfplot$marker<-xcell.data.tumor.filter[which(clinical.patient.xcell[,clinical.var]!=""),which(p.markers<0.05)[i]]
        tiff(paste0("Results/boxplot_",clinical.var,"_",colnames(xcell.data.tumor.filter)[which(p.markers<0.05)][i],".tiff"),res=300,h=2000,w=3000)
        print(ggplot(dfplot,aes(y=marker, fill=clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var], 
                                x=clinical.patient.xcell[which(clinical.patient.xcell[,clinical.var]!=""),clinical.var])) + geom_point() + geom_smooth(method='lm')
              + scale_y_continuous(name= colnames(xcell.data.tumor.filter)[which(p.markers<0.05)][i]) + stat_cor(method = "pearson") + scale_x_continuous(name=clinical.var)
              + scale_fill_continuous(name=clinical.var))
        dev.off()
      }
    }
  }
  return("Done")
}

##Histological type
clinical.patient.xcell$histological_type_2cat<-factor(ifelse(clinical.patient.xcell$histological_type=="Pancreas-Adenocarcinoma Ductal Type","PDAC","Other"))
association.test.xcell(clinical.patient.xcell,"histological_type_2cat",xcell.data.tumor.filter)


##anatomic_neoplasm_subdivision
clinical.patient.xcell$anatomic_neoplasm_subdivision<-factor(clinical.patient.xcell$anatomic_neoplasm_subdivision)
association.test.xcell(clinical.patient.xcell,"anatomic_neoplasm_subdivision",xcell.data.tumor.filter)

##gender
clinical.patient.xcell$gender<-factor(clinical.patient.xcell$gender)
association.test.xcell(clinical.patient.xcell,"gender",xcell.data.tumor.filter)

##race_list
clinical.patient.xcell$race_list<-factor(clinical.patient.xcell$race_list)
association.test.xcell(clinical.patient.xcell,"race_list",xcell.data.tumor.filter)

##History of Prior Malignancy
clinical.patient.xcell$other_dx<-factor(clinical.patient.xcell$other_dx)
association.test.xcell(clinical.patient.xcell,"other_dx",xcell.data.tumor.filter)

##number_of_lymphnodes_positive_by_he
association.test.xcell(clinical.patient.xcell,"number_of_lymphnodes_positive_by_he",xcell.data.tumor.filter)

##neoplasm_histologic_grade
clinical.patient.xcell$neoplasm_histologic_grade_3cat<-factor(ifelse(clinical.patient.xcell$neoplasm_histologic_grade=="G1","G1",
                                                                     ifelse(clinical.patient.xcell$neoplasm_histologic_grade=="G2","G2",
                                                                            ifelse(clinical.patient.xcell$neoplasm_histologic_grade=="G3","G3",""))))
association.test.xcell(clinical.patient.xcell,"neoplasm_histologic_grade_3cat",xcell.data.tumor.filter)

##Age 
association.test.xcell(clinical.patient.xcell,"age_at_initial_pathologic_diagnosis",xcell.data.tumor.filter)

## Heatmap ####
annotation_col = data.frame(
  vital_status = factor(clinical.patient.xcell$vital_status))
rownames(annotation_col)<-rownames(xcell.data.tumor.filter)
xcell.data.tumor.filter_vital_status<-xcell.data.tumor.filter[,c("aDC","Adipocytes","CD4+ naive T-cells","CD8+ T-cells","cDC",
                                                                 "DC","Epithelial cells","Fibroblasts","HSC","ImmuneScore",
                                                                 "Keratinocytes","Melanocytes","MicroenvironmentScore","Monocytes",
                                                                 "NKT","StromaScore")]

ann_colors = list (vital_status = c("Alive" = brewer.pal(3,"Set2")[1], "Dead"= brewer.pal(3,"Set2")[2]))
pheatmap(t(xcell.data.tumor.filter_vital_status),scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12),
         annotation_col = annotation_col,annotation_colors = ann_colors)



##Smoking
clinical.patient.xcell$smoking<-factor(ifelse(clinical.patient.xcell$tobacco_smoking_history_master=="Current smoker (includes daily smokers and non-daily smokers or occasional smokers)","Current",
                                              ifelse(clinical.patient.xcell$tobacco_smoking_history_master=="Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime)","Non-smoker","Former")))
association.test.xcell(clinical.patient.xcell,"smoking",xcell.data.tumor.filter)

##number_pack_years_smoked
association.test.xcell(clinical.patient.xcell,"number_pack_years_smoked",xcell.data.tumor.filter)

##Alcohol
clinical.patient.xcell$alcohol_history_documented<-factor(clinical.patient.xcell$alcohol_history_documented)
association.test.xcell(clinical.patient.xcell,"alcohol_history_documented",xcell.data.tumor.filter)

##Alcohol category
clinical.patient.xcell$alcoholic_exposure_category2<-ifelse(clinical.patient.xcell$alcohol_history_documented=="NO","No-drinker",
                                                            ifelse(clinical.patient.xcell$alcohol_history_documented=="YES" & clinical.patient.xcell$alcoholic_exposure_category=="",NA,
                                                                   ifelse(clinical.patient.xcell$alcoholic_exposure_category=="None","None-Drinker",
                                                                          ifelse(clinical.patient.xcell$alcoholic_exposure_category=="Occasional Drinker","Occasional-Drinker",
                                                                                 ifelse(clinical.patient.xcell$alcoholic_exposure_category=="Daily Drinker","Daily-Drinker",
                                                                                        ifelse(clinical.patient.xcell$alcoholic_exposure_category=="Social Drinker","Social-Drinker",
                                                                                               ifelse(clinical.patient.xcell$alcoholic_exposure_category=="Weekly Drinker","Weekly-Drinker",NA)))))))
clinical.patient.xcell$alcoholic_exposure_category2<-factor(clinical.patient.xcell$alcoholic_exposure_category2)
association.test.xcell(clinical.patient.xcell,"alcoholic_exposure_category2",xcell.data.tumor.filter)

##family history
clinical.patient.xcell$family_history_of_cancer<-factor(clinical.patient.xcell$family_history_of_cancer)
association.test.xcell(clinical.patient.xcell,"family_history_of_cancer",xcell.data.tumor.filter)

##radiation_therapy
clinical.patient.xcell$radiation_therapy<-factor(clinical.patient.xcell$radiation_therapy)
association.test.xcell(clinical.patient.xcell,"radiation_therapy",xcell.data.tumor.filter)

##primary_therapy_outcome_success
clinical.patient.xcell$primary_therapy_outcome_success<-factor(clinical.patient.xcell$primary_therapy_outcome_success)
association.test.xcell(clinical.patient.xcell,"primary_therapy_outcome_success",xcell.data.tumor.filter)

##history_chronic_pancreatitis
clinical.patient.xcell$history_of_chronic_pancreatitis<-factor(clinical.patient.xcell$history_of_chronic_pancreatitis)
association.test.xcell(clinical.patient.xcell,"history_of_chronic_pancreatitis",xcell.data.tumor.filter)

#Pathologic_stage
clinical.patient.xcell$pathologic_stage<-factor(ifelse(clinical.patient.xcell$stage_event_pathologic_stage=="Stage IA" | 
                                                         clinical.patient.xcell$stage_event_pathologic_stage=="Stage IB", "Stage I",
                                                       ifelse(clinical.patient.xcell$stage_event_pathologic_stage == "Stage IIA" |
                                                                clinical.patient.xcell$stage_event_pathologic_stage=="Stage IIB","Stage II",
                                                              ifelse(clinical.patient.xcell$stage_event_pathologic_stage=="Stage III","Stage III",
                                                                     ifelse(clinical.patient.xcell$stage_event_pathologic_stage=="Stage IV", "Stage IV",NA)))))
association.test.xcell(clinical.patient.xcell,"pathologic_stage",xcell.data.tumor.filter)

#################################
#######Clinical follow-up########
#################################
##vital_status
clinical.follow.up.xcell$vital_status<-factor(clinical.follow.up.xcell$vital_status)
association.test.xcell(clinical.follow.up.xcell,"vital_status",xcell.data.tumor.filter)

#treatment_outcome_first_course
clinical.follow.up.xcell$treatment_outcome_first_course<-replace(clinical.follow.up.xcell$treatment_outcome_first_course,
                                                                 clinical.follow.up.xcell$treatment_outcome_first_course=="[Discrepancy]" |
                                                                   clinical.follow.up.xcell$treatment_outcome_first_course=="[Not Applicable]" | 
                                                                   clinical.follow.up.xcell$treatment_outcome_first_course=="[Not Available]" |
                                                                   clinical.follow.up.xcell$treatment_outcome_first_course=="[Unknown]",NA)
clinical.follow.up.xcell$treatment_outcome_first_course<-factor(clinical.follow.up.xcell$treatment_outcome_first_course)
association.test.xcell(clinical.follow.up.xcell,"treatment_outcome_first_course",xcell.data.tumor.filter)

##new tumor event
clinical.follow.up.xcell$new_tumor_event_type<-replace(clinical.follow.up.xcell$new_tumor_event_type,clinical.follow.up.xcell$new_tumor_event_type=="#N/A" |
                                                       clinical.follow.up.xcell$new_tumor_event_type=="Locoregional Recurrence|Distant Metastasis" | 
                                                         clinical.follow.up.xcell$new_tumor_event_type=="New Primary Tumor",NA)
clinical.follow.up.xcell$new_tumor_event_type<-factor(clinical.follow.up.xcell$new_tumor_event_type)
association.test.xcell(clinical.follow.up.xcell,"new_tumor_event_type",xcell.data.tumor.filter)

###########################
##### Biospecimen ########
##########################
biospecimen.slide.xcell<-biospecimen.slide[match(substr(rownames(xcell.data.tumor.filter),1,12),biospecimen.slide$bcr_patient_barcode),]

#percent_tumor_cells
association.test.xcell(biospecimen.slide.xcell,"percent_tumor_cells",xcell.data.tumor.filter)

#percent_tumor_nuclei
association.test.xcell(biospecimen.slide.xcell,"percent_tumor_nuclei",xcell.data.tumor.filter)

#percent_stromal_cells
association.test.xcell(biospecimen.slide.xcell,"percent_stromal_cells",xcell.data.tumor.filter)

#percent_lymphocyte_infiltration
association.test.xcell(biospecimen.slide.xcell,"percent_lymphocyte_infiltration",xcell.data.tumor.filter)



####Xcell
clinical.follow.up.xcell<-clinical.folow_up[match(substr(rownames(xcell.data.tumor.filter),1,12),clinical.folow_up$bcr_patient_barcode),]
##OS
surv_object <- Surv(time = clinical.follow.up.xcell$OS.time, event = clinical.follow.up.xcell$OS)
res.cox <- coxph(surv_object~xcell.data.tumor.filter[,"NKT"])
summary(res.cox)
##Categorical
KL_mean<-mean(xcell.data.tumor.filter[,"NKT"])
NKT_2cat<-ifelse(xcell.data.tumor.filter[,"NKT"]<=KL_mean,1,2)
fit1 <- survfit(surv_object ~ NKT_2cat)
fit1
ggsurvplot(fit1, data = as.data.frame(xcell.data.tumor.filter))
comp(ten(fit1))$tests$lrTests

#DSS
clinical.follow.up.xcell.DSS<-clinical.follow.up.xcell[which(clinical.follow.up.xcell$DSS!="#N/A"),]
clinical.follow.up.xcell.DSS$DSS<-as.integer(as.character(clinical.follow.up.xcell.DSS$DSS))
surv_object <- Surv(time = clinical.follow.up.xcell.DSS$DSS.time, event = clinical.follow.up.xcell.DSS$DSS)
res.cox <- coxph(surv_object~xcell.data.tumor.filter[,"NKT"][which(clinical.follow.up.xcell.DSS$DSS!="#N/A")])
summary(res.cox)
#Categorical
KL_mean<-mean(xcell.data.tumor.filter[,"NKT"])
NKT_2cat<-ifelse(xcell.data.tumor.filter[,"NKT"]<=KL_mean,1,2)
fit1 <- survfit(surv_object ~ NKT_2cat[which(clinical.follow.up.xcell.DSS$DSS!="#N/A")])
fit1
ggsurvplot(fit1, data = as.data.frame(xcell.data.tumor.filter[which(clinical.follow.up.xcell.DSS$DSS!="#N/A"),]))
comp(ten(fit1))$tests$lrTests


##PFI
surv_object <- Surv(time = clinical.follow_up.tumor$PFI.time, event = clinical.follow_up.tumor$PFI)
res.cox <- coxph(surv_object~PAAD.repertoire.tumor$KappaLambda_ratio_expression)
summary(res.cox)






###Cluster
##Cluster
#id<-match(clinical.follow_up.tumor$bcr_patient_barcode,substr(names(fit.k$cluster),1,12))
#clinical.follow_up.tumor$cluster<-factor(fit.k$cluster[id])

##OS
surv_object <- Surv(time = clinical.follow_up.tumor$OS.time, event = clinical.follow_up.tumor$OS)
fit1 <- survfit(surv_object ~ clinical.follow_up.tumor$cluster, data = clinical.follow_up.tumor)
fit1
tiff("Results/KM_cluster_Kmeans.tiff",res=300,h=2000,w=3000)
ggsurvplot(fit1, data = clinical.follow_up.tumor)
dev.off()
comp(ten(fit1))$tests$lrTests

#DSS
clinical.follow_up.tumor2<-clinical.follow_up.tumor[which(clinical.follow_up.tumor$DSS!="#N/A"),]
clinical.follow_up.tumor2$DSS<-as.integer(as.character(clinical.follow_up.tumor2$DSS))
surv_object <- Surv(time = clinical.follow_up.tumor2$DSS.time, event = clinical.follow_up.tumor2$DSS)
fit1 <- survfit(surv_object ~ clinical.follow_up.tumor2$cluster, data = clinical.follow_up.tumor2)
ggsurvplot(fit1, data = clinical.follow_up.tumor2)
comp(ten(fit1))$tests$lrTests

##PFI
surv_object <- Surv(time = clinical.follow_up.tumor$PFI.time, event = clinical.follow_up.tumor$PFI)
fit1 <- survfit(surv_object ~ clinical.follow_up.tumor$cluster, data = clinical.follow_up.tumor)
summary(fit1)
tiff("Results/KM_cluster_Kmeans.tiff",res=300,h=2000,w=3000)
ggsurvplot(fit1, data = clinical.follow_up.tumor)
dev.off()
comp(ten(fit1))$tests$lrTests


################################
### Compare data with GTEX ####
##############################
load("Data/GTEx/Blood/MIXCR/GTEX_FullData.Rdata")

repertoire.diversity<-rbind(PAAD.repertoire.diversity[,c("Total_Reads","IGH_Reads","IGK_Reads", "IGL_Reads","TRA_Reads","TRB_Reads","TRD_Reads","TRG_Reads","IG_Reads",
                                                         "T_Reads","IG_expression", "IGH_expression", "IGK_expression","IGL_expression","T_expression","TRA_expression",
                                                         "TRB_expression","TRD_expression","TRG_expression","Alpha_Beta_ratio_expression","KappaLambda_ratio_expression",
                                                         "clones_IGH","clones_IGK","clones_IGL", "clones_TRA", "clones_TRB","clones_TRD","clones_TRG","entropy_IGH",
                                                         "entropy_IGL","entropy_TRA","entropy_TRB","entropy_TRD","entropy_TRG")],
                            GTEX.repertoire.diversity[,c("Total_Reads","IGH_Reads","IGK_Reads", "IGL_Reads","TRA_Reads","TRB_Reads","TRD_Reads","TRG_Reads","IG_Reads",
                                                         "T_Reads","IG_expression", "IGH_expression", "IGK_expression","IGL_expression","T_expression","TRA_expression",
                                                         "TRB_expression","TRD_expression","TRG_expression","Alpha_Beta_ratio_expression","KappaLambda_ratio_expression",
                                                         "clones_IGH","clones_IGK","clones_IGL", "clones_TRA", "clones_TRB","clones_TRD","clones_TRG","entropy_IGH",
                                                         "entropy_IGL","entropy_TRA","entropy_TRB","entropy_TRD","entropy_TRG")])
repertoire.diversity$sample_type<-factor(c(rep("PAAD",160),rep("GTEX",441)))
repertoire.diversity$TotalSeqReads<-repertoire.diversity$IG_Reads/repertoire.diversity$IG_expression

ggplot(repertoire.diversity,aes(y=entropy_TRA, fill=sample_type, x=sample_type)) + geom_boxplot()  + stat_compare_means()


