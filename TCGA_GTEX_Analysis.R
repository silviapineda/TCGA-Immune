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
library(reshape)
library(stringr)

setwd("~/TCGA-Immune/")

load("Data/PAAD/PAAD_FullData.Rdata")
load("Data/GTEx/Pancreas/GTEx_FullData.Rdata")
load("Data/Pancreas_Validation/Pancreas_Validation_FullData.Rdata")
load("Data/Validation_Normal_pancreas/Pancreas_Normal_Validation_FullData.Rdata")


###################################################
#### Immune reperotire Analysis TCGA PDAC only ###
##################################################

#### TCGA - Only PDAC samples
PAAD.repertoire.tumor<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type_4categ=="PDAC"),] #144
PAAD.repertoire.tumor$TCGA_sample<-substr(PAAD.repertoire.tumor$TCGA_sample,1,15)


##############################
### Merge with subtypes #####
##############################
paad.subtype.tumor<-paad.subtype[match(substr(PAAD.repertoire.tumor$TCGA_sample,1,12),substr(paad.subtype$Tumor.Sample.ID,1,12)),]
PAAD.repertoire.tumor.subtype<-cbind(PAAD.repertoire.tumor,paad.subtype.tumor)

###############################
###Merged with mutated genes ##
###############################
paad.mutated.genes<-read.csv("Data/PAAD/Clinical/mutated.genes.csv")
gene_list<-NULL
for(i in 1:ncol(paad.mutated.genes)){
  gene_list<-c(gene_list,as.character(paad.mutated.genes[,i]))
}
gene_list<-unique(gene_list)
gene_list<-gene_list[order(gene_list)]
gene_list<-gene_list[-c(1,536)]

paad.mutated.genes.total<-matrix(NA,150,length(gene_list))
for(i in 1:150){
  id<-match(paad.mutated.genes[,i],gene_list)
  id.gene<-unique(na.omit(id))
  paad.mutated.genes.total[i,id.gene]<-1
  paad.mutated.genes.total[i,]<-replace(paad.mutated.genes.total[i,],is.na(paad.mutated.genes.total[i,])==T,0)
}
rownames(paad.mutated.genes.total)<-colnames(paad.mutated.genes)
colnames(paad.mutated.genes.total)<-gene_list
rownames(paad.mutated.genes.total)<-str_replace(rownames(paad.mutated.genes.total),"\\.","-")
rownames(paad.mutated.genes.total)<-str_replace(rownames(paad.mutated.genes.total),"\\.","-")
rownames(paad.mutated.genes.total)<-str_replace(rownames(paad.mutated.genes.total),"\\.","-")

paad.mutated.genes.tumor<-paad.mutated.genes.total[match(substr(PAAD.repertoire.tumor$TCGA_sample,1,12),substr(rownames(paad.mutated.genes.total),1,12)),]
PAAD.repertoire.tumor.subtype<-cbind(PAAD.repertoire.tumor,paad.subtype.tumor,paad.mutated.genes.tumor)
PAAD.repertoire.tumor.subtype$Mutated.Genes<-rowSums(paad.mutated.genes.tumor)

###############################
## Merge with Clinical data ###
###############################
clinical.patient.tumor<-clinical.patient[match(substr(PAAD.repertoire.tumor.subtype$TCGA_sample,1,12),clinical.patient$bcr_patient_barcode),]
PAAD.repertoire.tumor.clinical.patient<-cbind(PAAD.repertoire.tumor.subtype,clinical.patient.tumor)

#################################
#######Clinical follow-up########
#################################
clinical.follow_up.tumor<-clinical.folow_up[match(substr(PAAD.repertoire.tumor.clinical.patient$TCGA_sample,1,12),clinical.folow_up$bcr_patient_barcode),]
PAAD.repertoire.tumor.clinical.followuop<-cbind(PAAD.repertoire.tumor.clinical.patient,clinical.follow_up.tumor)

##############################################
## Merge with ImmuneFeaturesLandscapePAAD ###
#############################################
ImmuneFeaturesLandscapePAAD<-read.csv("Data/PAAD/ImmuneFeaturesLandscapePAAD.csv")
ImmuneFeaturesLandscape.tumor<-ImmuneFeaturesLandscapePAAD[match(substr(PAAD.repertoire.tumor.clinical.followuop$TCGA_sample,1,12),ImmuneFeaturesLandscapePAAD$TCGA.Participant.Barcode),]
PAAD.repertoire.tumor.ImmuneFeaturesLandscape<-cbind(PAAD.repertoire.tumor.clinical.followuop,ImmuneFeaturesLandscape.tumor)

#############################################
### Final data frame with all variables #####
############################################
PAAD.repertoire.tumor<-PAAD.repertoire.tumor.ImmuneFeaturesLandscape



#### Plot all the variables with the outlier 
tiff("Results/Clinical/Mutated_Outlier.tiff",res=300,h=2500,w=500)
par(mfrow = c(5,1))
col<-brewer.pal(9, "Set1")
barplot(PAAD.repertoire.tumor$IGH,col=col[1], main = "IgH reads")
barplot(PAAD.repertoire.tumor$Mutated.Genes,col=col[2],main = "Mutated Genes")
barplot(PAAD.repertoire.tumor$SNV.Neoantigens,col=col[3],main = "SNV neoantigens")
barplot(PAAD.repertoire.tumor$Nonsilent.Mutation.Rate,col=col[4],main = "NonSilent mutation rate")
barplot(PAAD.repertoire.tumor$Silent.Mutation.Rate,col=col[5],main = "Silent mutation rate")
dev.off()

### 5e9e81e2-0e8b-4aca-aced-6ce451fa3262
grep("5e9e81e2-0e8b-4aca-aced-6ce451fa3262",rownames(PAAD.repertoire.tumor))
PAAD.repertoire.tumor.filter<-PAAD.repertoire.tumor[-53,]

save(data_merge,PAAD.repertoire.diversity,PAAD.repertoire.tumor.filter,xCell.data.PAAD,xCell.pvalue.PAAD,paad.subtype,clinical.drug,clinical.patient,clinical.radiation,clinical.new_tumor_event,clinical.folow_up,biospecimen.slide,annotation,
     file="Data/PAAD/PAAD_FullData.Rdata")


##############################################################
# 1. Association analysis with subtypes and mutation genes ###
##############################################################

##Mutated.Genes
summary(glm(PAAD.repertoire.tumor.filter$Mutated.Genes~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC))
PAAD.repertoire.tumor.filter$Mutated.Genes<-factor(PAAD.repertoire.tumor.filter$Mutated.Genes)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Mutated.Genes")

##ABSOLUTE PURITY
plot(PAAD.repertoire.tumor.filter$ABSOLUTE.Purity~PAAD.repertoire.tumor.filter$Risk_Score_compositional)
summary(glm(PAAD.repertoire.tumor.filter$Mutated.Genes~factor(PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC)))

PAAD.repertoire.tumor.filter$ABSOLUTE.Purity<-factor(PAAD.repertoire.tumor.filter$ABSOLUTE.Purity)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"ABSOLUTE.Purity")

##Ploidy
PAAD.repertoire.tumor.filter$Ploidy<-factor(PAAD.repertoire.tumor.filter$Ploidy)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Ploidy") #No significant

##Purity.Class..high.or.low.
PAAD.repertoire.tumor.filter$Purity<-factor(PAAD.repertoire.tumor.filter$Purity.Class..high.or.low.)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Purity")

##mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical
PAAD.repertoire.tumor.filter$Moffitt.clusters..All.150.Samples..1basal..2classical<-factor(PAAD.repertoire.tumor.filter$mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical)
PAAD.repertoire.tumor.filter$subtypes_Moffit<-factor(ifelse(PAAD.repertoire.tumor.filter$Moffitt.clusters..All.150.Samples..1basal..2classical==1,"Basal",
                                                             ifelse(PAAD.repertoire.tumor.filter$Moffitt.clusters..All.150.Samples..1basal..2classical==2,"Classical",NA)))
association.test.immuneRep(PAAD.repertoire.tumor.filter,"subtypes_Moffit")

chisq.test(PAAD.repertoire.tumor.filter$subtypes_Moffit,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC)
summary(glm(PAAD.repertoire.tumor.filter$Risk_Score_compositional~factor(PAAD.repertoire.tumor.filter$subtypes_Moffit)))

##mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM
PAAD.repertoire.tumor.filter$mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM<-factor(PAAD.repertoire.tumor.filter$mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM)
PAAD.repertoire.tumor.filter$subtypes_Collisson<-factor(ifelse(PAAD.repertoire.tumor.filter$mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM==1,"Classical",
                                                                ifelse(PAAD.repertoire.tumor.filter$mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM==2,"Exocrine",
                                                                       ifelse(PAAD.repertoire.tumor.filter$mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM==3,"QM",NA))))

association.test.immuneRep(PAAD.repertoire.tumor.filter,"subtypes_Collisson")
chisq.test(PAAD.repertoire.tumor.filter$subtypes_Collisson,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
summary(glm(PAAD.repertoire.tumor.filter$Risk_Score_compositional~factor(PAAD.repertoire.tumor.filter$subtypes_Collisson)))
boxplot(PAAD.repertoire.tumor.filter$Risk_Score_compositional~factor(PAAD.repertoire.tumor.filter$subtypes_Collisson))

##mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX
PAAD.repertoire.tumor.filter$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX<-factor(PAAD.repertoire.tumor.filter$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX)
PAAD.repertoire.tumor.filter$subtypes_Bailey<-factor(ifelse(PAAD.repertoire.tumor.filter$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX==1,"Squamous",
                                                             ifelse(PAAD.repertoire.tumor.filter$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX==2,"Immunogenic",
                                                                    ifelse(PAAD.repertoire.tumor.filter$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX==3,"Progenitor",
                                                                           ifelse(PAAD.repertoire.tumor.filter$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX==4,"ADEX",NA)))))
association.test.immuneRep(PAAD.repertoire.tumor.filter,"subtypes_Bailey")
chisq.test(PAAD.repertoire.tumor.filter$subtypes_Bailey,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC)
summary(glm(PAAD.repertoire.tumor.filter$Risk_Score_compositional~factor(PAAD.repertoire.tumor.filter$subtypes_Bailey)))
boxplot(PAAD.repertoire.tumor.filter$Risk_Score_compositional~factor(PAAD.repertoire.tumor.filter$subtypes_Bailey))

##Copy.Number.Clusters..All.150.Samples.
PAAD.repertoire.tumor.filter$Copy.Number.Clusters..All.150.Samples.<-factor(PAAD.repertoire.tumor.filter$Copy.Number.Clusters..All.150.Samples.)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Copy.Number.Clusters..All.150.Samples.")

##KRAS.Mutated..1.or.0.
chisq.test(PAAD.repertoire.tumor.filter$KRAS.Mutated..1.or.0.,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC)
summary(glm(PAAD.repertoire.tumor.filter$Risk_Score_compositional~factor(PAAD.repertoire.tumor.filter$KRAS.Mutated..1.or.0.)))
PAAD.repertoire.tumor.filter$KRAS.Mutated..1.or.0.<-factor(PAAD.repertoire.tumor.filter$KRAS.Mutated..1.or.0.)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"KRAS.Mutated..1.or.0.")

##CDKN2A.Expression
PAAD.repertoire.tumor.filter$CDKN2A.Expression<-factor(PAAD.repertoire.tumor.filter$CDKN2A.Expression)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"CDKN2A.Expression")

##DNA.methylation.leukocyte.percent.estimate
PAAD.repertoire.tumor.filter$Leukocyte_DNAmethylation<-factor(PAAD.repertoire.tumor.filter$DNA.methylation.leukocyte.percent.estimate)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Leukocyte_DNAmethylation")

##DNA.hypermethylation.mode.purity
PAAD.repertoire.tumor.filter$DNA.hypermethylation.mode.purity<-factor(PAAD.repertoire.tumor.filter$DNA.hypermethylation.mode.purity)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"DNA.hypermethylation.mode.purity")

##miRNA.Clusters..All.150.Samples.
PAAD.repertoire.tumor.filter$miRNA.Clusters..All.150.Samples.<-factor(PAAD.repertoire.tumor.filter$miRNA.Clusters..All.150.Samples.)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"miRNA.Clusters..All.150.Samples.")

##miRNA.Clusters..All.150.Samples.
PAAD.repertoire.tumor.filter$lncRNA.Clusters..All.150.Samples.<-factor(PAAD.repertoire.tumor.filter$lncRNA.Clusters..All.150.Samples.)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"lncRNA.Clusters..All.150.Samples.")

##ZNF521 mutattion
chisq.test(PAAD.repertoire.tumor.filter$ZNF521.mutated,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC)
boxplot(PAAD.repertoire.tumor.filter$Risk_Score_compositional~factor(PAAD.repertoire.tumor.filter$ZNF521.mutated))
summary(glm(PAAD.repertoire.tumor.filter$Risk_Score_compositional~factor(PAAD.repertoire.tumor.filter$ZNF521.mutated)))
table(PAAD.repertoire.tumor.filter$ZNF521.mutated,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
PAAD.repertoire.tumor.filter$ZNF521.mutated<-factor(PAAD.repertoire.tumor.filter$ZNF521)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"ZNF521.mutated")

##############################################################
# 2. Association analysis with Immune Landscape ###
##############################################################
PAAD.repertoire.tumor.filter$Immune.Subtype<-factor(PAAD.repertoire.tumor.filter$Immune.Subtype)
chisq.test(PAAD.repertoire.tumor.filter$Immune.Subtype,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Immune.Subtype")

PAAD.repertoire.tumor.filter$Aneuploidy.Score<-factor(PAAD.repertoire.tumor.filter$Aneuploidy.Score)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Aneuploidy.Score")

boxplot(PAAD.repertoire.tumor.filter$Leukocyte.Fraction~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC)
plot(PAAD.repertoire.tumor.filter$Leukocyte.Fraction~PAAD.repertoire.tumor.filter$Risk_Score_compositional)
summary(glm(PAAD.repertoire.tumor.filter$Leukocyte.Fraction~PAAD.repertoire.tumor.filter$Risk_Score_compositional))

PAAD.repertoire.tumor.filter$Leukocyte.Fraction<-factor(PAAD.repertoire.tumor.filter$Leukocyte.Fraction)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Leukocyte.Fraction")

boxplot(PAAD.repertoire.tumor.filter$Stromal.Fraction~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
PAAD.repertoire.tumor.filter$Stromal.Fraction<-factor(PAAD.repertoire.tumor.filter$Stromal.Fraction)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Stromal.Fraction")

PAAD.repertoire.tumor.filter$Intratumor.Heterogeneity<-factor(PAAD.repertoire.tumor.filter$Intratumor.Heterogeneity)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Intratumor.Heterogeneity")

PAAD.repertoire.tumor.filter$TIL.Regional.Fraction<-factor(PAAD.repertoire.tumor.filter$TIL.Regional.Fraction)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"TIL.Regional.Fraction")

PAAD.repertoire.tumor.filter$Proliferation<-factor(PAAD.repertoire.tumor.filter$Proliferation)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Proliferation")

boxplot(PAAD.repertoire.tumor.filter$Wound.Healing~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
PAAD.repertoire.tumor.filter$Wound.Healing<-factor(PAAD.repertoire.tumor.filter$Wound.Healing)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Wound.Healing")

boxplot(PAAD.repertoire.tumor.filter$Macrophage.Regulation~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
PAAD.repertoire.tumor.filter$Macrophage.Regulation<-factor(PAAD.repertoire.tumor.filter$Macrophage.Regulation)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Macrophage.Regulation")

PAAD.repertoire.tumor.filter$Lymphocyte.Infiltration.Signature.Score<-factor(PAAD.repertoire.tumor.filter$Lymphocyte.Infiltration.Signature.Score)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Lymphocyte.Infiltration.Signature.Score")

boxplot(PAAD.repertoire.tumor.filter$IFN.gamma.Response~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
PAAD.repertoire.tumor.filter$IFN.gamma.Response<-factor(PAAD.repertoire.tumor.filter$IFN.gamma.Response)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"IFN.gamma.Response")

boxplot(PAAD.repertoire.tumor.filter$TGF.beta.Response~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
summary(glm(PAAD.repertoire.tumor.filter$TGF.beta.Response~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster))
PAAD.repertoire.tumor.filter$TGF.beta.Response<-factor(PAAD.repertoire.tumor.filter$TGF.beta.Response)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"TGF.beta.Response")

boxplot(PAAD.repertoire.tumor.filter$SNV.Neoantigens~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
plot(PAAD.repertoire.tumor.filter$SNV.Neoantigens~PAAD.repertoire.tumor.filter$Risk_Score_compositional)
summary(glm(PAAD.repertoire.tumor.filter$SNV.Neoantigens~PAAD.repertoire.tumor.filter$Risk_Score_compositional))

PAAD.repertoire.tumor.filter$SNV.Neoantigens<-factor(PAAD.repertoire.tumor.filter$SNV.Neoantigens)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"SNV.Neoantigens")

boxplot(PAAD.repertoire.tumor.filter$Indel.Neoantigens~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
PAAD.repertoire.tumor.filter$Indel.Neoantigens<-factor(PAAD.repertoire.tumor.filter$Indel.Neoantigens)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Indel.Neoantigens")

boxplot(PAAD.repertoire.tumor.filter$Silent.Mutation.Rate~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
PAAD.repertoire.tumor.filter$Silent.Mutation.Rate<-factor(PAAD.repertoire.tumor.filter$Silent.Mutation.Rate)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Silent.Mutation.Rate")

boxplot(PAAD.repertoire.tumor.filter$Nonsilent.Mutation.Rate~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
PAAD.repertoire.tumor.filter$Nonsilent.Mutation.Rate<-factor(PAAD.repertoire.tumor.filter$Nonsilent.Mutation.Rate)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Nonsilent.Mutation.Rate")

PAAD.repertoire.tumor.filter$BCR.Evenness<-factor(PAAD.repertoire.tumor.filter$BCR.Evenness)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"BCR.Evenness")

PAAD.repertoire.tumor.filter$BCR.Shannon<-factor(PAAD.repertoire.tumor.filter$BCR.Shannon)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"BCR.Shannon")

PAAD.repertoire.tumor.filter$BCR.Richness<-factor(PAAD.repertoire.tumor.filter$BCR.Richness)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"BCR.Richness")

PAAD.repertoire.tumor.filter$TCR.Evenness<-factor(PAAD.repertoire.tumor.filter$TCR.Evenness)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"TCR.Evenness")

PAAD.repertoire.tumor.filter$TCR.Shannon<-factor(PAAD.repertoire.tumor.filter$TCR.Shannon)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"TCR.Shannon")

PAAD.repertoire.tumor.filter$TCR.Richness<-factor(PAAD.repertoire.tumor.filter$TCR.Richness)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"TCR.Richness")

PAAD.repertoire.tumor.filter$TCR.Richness<-factor(PAAD.repertoire.tumor.filter$TCR.Richness)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"TCR.Richness")


boxplot(PAAD.repertoire.tumor.filter$Th1.Cells~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
plot(PAAD.repertoire.tumor.filter$Th1.Cells~PAAD.repertoire.tumor.filter$Risk_Score_compositional)
summary(glm(PAAD.repertoire.tumor.filter$Th17.Cells~PAAD.repertoire.tumor.filter$Risk_Score_compositional))

PAAD.repertoire.tumor.filter$Th1.Cells<-factor(PAAD.repertoire.tumor.filter$Th1.Cells)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Th1.Cells")

boxplot(PAAD.repertoire.tumor.filter$Th2.Cells~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
PAAD.repertoire.tumor.filter$Th2.Cells<-factor(PAAD.repertoire.tumor.filter$Th2.Cells)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Th2.Cells")


boxplot(PAAD.repertoire.tumor.filter$Th17.Cells~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
PAAD.repertoire.tumor.filter$Th17.Cells<-factor(PAAD.repertoire.tumor.filter$Th17.Cells)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"Th17.Cells")

##############################################################
# 3. Association analysis with Clinical outcomes ###
##############################################################

##Histological type
PAAD.repertoire.tumor.filter$histological_type<-factor(PAAD.repertoire.tumor.filter$histological_type)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"histological_type")

##anatomic_neoplasm_subdivision
PAAD.repertoire.tumor.filter$anatomic_neoplasm_subdivision<-factor(PAAD.repertoire.tumor.filter$anatomic_neoplasm_subdivision)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"anatomic_neoplasm_subdivision")

##gender
PAAD.repertoire.tumor.filter$gender<-factor(PAAD.repertoire.tumor.filter$gender)
chisq.test(PAAD.repertoire.tumor.filter$gender,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC)
summary(glm(PAAD.repertoire.tumor.filter$Risk_Score_compositional~PAAD.repertoire.tumor.filter$gender))

association.test.immuneRep(PAAD.repertoire.tumor.filter,"gender")

##race_list
PAAD.repertoire.tumor.filter$race_list<-ifelse(PAAD.repertoire.tumor.filter$race_list=="",NA,as.character(PAAD.repertoire.tumor.filter$race_list))
chisq.test(PAAD.repertoire.tumor.filter$race_list,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"race_list")

##History of Prior Malignancy
PAAD.repertoire.tumor.filter$other_dx<-factor(PAAD.repertoire.tumor.filter$other_dx)
chisq.test(PAAD.repertoire.tumor.filter$other_dx,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC)
summary(glm(PAAD.repertoire.tumor.filter$Risk_Score_compositional~PAAD.repertoire.tumor.filter$other_dx))

association.test.immuneRep(PAAD.repertoire.tumor.filter,"other_dx")

##neoplasm_histologic_grade
PAAD.repertoire.tumor.filter$neoplasm_histologic_grade_3cat<-factor(ifelse(PAAD.repertoire.tumor.filter$neoplasm_histologic_grade=="G1","G1",
                                                                                     ifelse(PAAD.repertoire.tumor.filter$neoplasm_histologic_grade=="G2","G2",
  
                                                                                                                                                                                      ifelse(PAAD.repertoire.tumor.filter$neoplasm_histologic_grade=="G3","G3",NA))))
chisq.test(PAAD.repertoire.tumor.filter$neoplasm_histologic_grade_3cat,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"neoplasm_histologic_grade_3cat")

##Age 
boxplot(PAAD.repertoire.tumor.filter$age_at_initial_pathologic_diagnosis~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
PAAD.repertoire.tumor.filter$age_at_initial_pathologic_diagnosis<-factor(PAAD.repertoire.tumor.filter$age_at_initial_pathologic_diagnosis)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"age_at_initial_pathologic_diagnosis")


##Smoking
PAAD.repertoire.tumor.filter$smoking<-factor(ifelse(PAAD.repertoire.tumor.filter$tobacco_smoking_history_master=="Current smoker (includes daily smokers and non-daily smokers or occasional smokers)","Current",
                                                              ifelse(PAAD.repertoire.tumor.filter$tobacco_smoking_history_master=="Lifelong Non-smoker (less than 100 cigarettes smoked in Lifetime)","Non-smoker",
                                                                     ifelse(PAAD.repertoire.tumor.filter$tobacco_smoking_history_master=="Current reformed smoker for > 15 years (greater than 15 years)","Former",
                                                                            ifelse(PAAD.repertoire.tumor.filter$tobacco_smoking_history_master=="Current reformed smoker for ≤15 years (less than or equal to 15 years)","Former",
                                                                                   ifelse(PAAD.repertoire.tumor.filter$tobacco_smoking_history_master=="Current reformed smoker, duration not specified","Former",NA))))))
chisq.test(PAAD.repertoire.tumor.filter$smoking,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC)
summary(glm(PAAD.repertoire.tumor.filter$Risk_Score_compositional~PAAD.repertoire.tumor.filter$smoking))

association.test.immuneRep(PAAD.repertoire.tumor.filter,"smoking")

##Smoking 2
PAAD.repertoire.tumor.filter$smoking2<-factor(ifelse(PAAD.repertoire.tumor.filter$smoking=="Current" |
                                                                 PAAD.repertoire.tumor.filter$smoking=="Former" ,"Ever-Smoker",
                                                               ifelse(PAAD.repertoire.tumor.filter$smoking=="Non-smoker","Non-smoker",NA)))
chisq.test(PAAD.repertoire.tumor.filter$smoking2,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"smoking2")

summary(glm(PAAD.repertoire.tumor.filter$Risk_Score_compositional~PAAD.repertoire.tumor.filter$smoking))

##number_pack_years_smoked
boxplot(PAAD.repertoire.tumor.filter$number_pack_years_smoked~PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
summary(glm(PAAD.repertoire.tumor.filter$Risk_Score_compositional~PAAD.repertoire.tumor.filter$number_pack_years_smoked))

PAAD.repertoire.tumor.filter$number_pack_years_smoked_all<-ifelse(PAAD.repertoire.tumor.filter$smoking2=="Non-smoker",0,PAAD.repertoire.tumor.filter$number_pack_years_smoked)
PAAD.repertoire.tumor.filter$number_pack_years_smoked_all<-factor(PAAD.repertoire.tumor.filter$number_pack_years_smoked_all)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"number_pack_years_smoked_all")

##Alcohol
PAAD.repertoire.tumor.filter$alcohol_history_documented<-ifelse(PAAD.repertoire.tumor.filter$alcohol_history_documented=="",NA,
                                                                          as.character(PAAD.repertoire.tumor.filter$alcohol_history_documented))
association.test.immuneRep(PAAD.repertoire.tumor.filter,"alcohol_history_documented")

##Alcohol category
PAAD.repertoire.tumor.filter$alcoholic_exposure_category2<-ifelse(PAAD.repertoire.tumor.filter$alcohol_history_documented=="NO","No-drinker",
                                                                            ifelse(PAAD.repertoire.tumor.filter$alcohol_history_documented=="YES" & PAAD.repertoire.tumor.filter$alcoholic_exposure_category=="",NA,
                                                                                   ifelse(PAAD.repertoire.tumor.filter$alcoholic_exposure_category=="None","None-Drinker",
                                                                                          ifelse(PAAD.repertoire.tumor.filter$alcoholic_exposure_category=="Occasional Drinker","Occasional-Drinker",
                                                                                                 ifelse(PAAD.repertoire.tumor.filter$alcoholic_exposure_category=="Daily Drinker","Daily-Drinker",
                                                                                                        ifelse(PAAD.repertoire.tumor.filter$alcoholic_exposure_category=="Social Drinker","Social-Drinker",
                                                                                                               ifelse(PAAD.repertoire.tumor.filter$alcoholic_exposure_category=="Weekly Drinker","Weekly-Drinker",NA)))))))
PAAD.repertoire.tumor.filter$alcoholic_exposure_category2<-factor(PAAD.repertoire.tumor.filter$alcoholic_exposure_category2)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"alcoholic_exposure_category2")

##family history
PAAD.repertoire.tumor.filter$family_history_of_cancer<-factor(ifelse(PAAD.repertoire.tumor.filter$family_history_of_cancer=="",NA,PAAD.repertoire.tumor.filter$family_history_of_cancer))
chisq.test(PAAD.repertoire.tumor.filter$family_history_of_cancer,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"family_history_of_cancer")

##radiation_therapy
PAAD.repertoire.tumor.filter$radiation_therapy<-factor(ifelse(PAAD.repertoire.tumor.filter$radiation_therapy=="",NA,
                                                                        PAAD.repertoire.tumor.filter$radiation_therapy))
chisq.test(PAAD.repertoire.tumor.filter$radiation_therapy,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)

association.test.immuneRep(PAAD.repertoire.tumor.filter,"radiation_therapy")

##primary_therapy_outcome_success
PAAD.repertoire.tumor.filter$primary_therapy_outcome_success<-ifelse(PAAD.repertoire.tumor.filter$primary_therapy_outcome_success=="",NA,
                                                                               PAAD.repertoire.tumor.filter$primary_therapy_outcome_success)
chisq.test(PAAD.repertoire.tumor.filter$primary_therapy_outcome_success,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"primary_therapy_outcome_success")

##history_chronic_pancreatitis
PAAD.repertoire.tumor.filter$history_of_chronic_pancreatitis<-factor(ifelse(PAAD.repertoire.tumor.filter$history_of_chronic_pancreatitis=="",NA,
                                                                                      PAAD.repertoire.tumor.filter$history_of_chronic_pancreatitis))
chisq.test(PAAD.repertoire.tumor.filter$history_of_chronic_pancreatitis,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"history_of_chronic_pancreatitis")

#stage_event_tnm_categories
PAAD.repertoire.tumor.filter$pathologic_stage<-factor(ifelse(PAAD.repertoire.tumor.filter$stage_event_pathologic_stage=="Stage IA" | 
                                                                         PAAD.repertoire.tumor.filter$stage_event_pathologic_stage=="Stage IB", "Stage I",
                                                                       ifelse(PAAD.repertoire.tumor.filter$stage_event_pathologic_stage == "Stage IIA" |
                                                                                PAAD.repertoire.tumor.filter$stage_event_pathologic_stage=="Stage IIB","Stage II",
                                                                              ifelse(PAAD.repertoire.tumor.filter$stage_event_pathologic_stage=="Stage III","Stage III-IV",
                                                                                     ifelse(PAAD.repertoire.tumor.filter$stage_event_pathologic_stage=="Stage IV", "Stage III-IV",NA)))))
chisq.test(PAAD.repertoire.tumor.filter$pathologic_stage,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"pathologic_stage")

##history_diabetes
PAAD.repertoire.tumor.filter$history_of_diabetes<-ifelse(PAAD.repertoire.tumor.filter$history_of_diabetes=="",NA,
                                                                   as.character(PAAD.repertoire.tumor.filter$history_of_diabetes))
chisq.test(PAAD.repertoire.tumor.filter$history_of_diabetes,PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"history_of_diabetes")



##vital_status
PAAD.repertoire.tumor.filter$vital_status<-factor(PAAD.repertoire.tumor.filter$vital_status)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"vital_status")

##new tumor event
PAAD.repertoire.tumor.filter$new_tumor_event_type<-replace(PAAD.repertoire.tumor.filter$new_tumor_event_type,PAAD.repertoire.tumor.filter$new_tumor_event_type=="#N/A",NA)
PAAD.repertoire.tumor.filter$new_tumor_event_type<-replace(PAAD.repertoire.tumor.filter$new_tumor_event_type,
                                                                       PAAD.repertoire.tumor.filter$new_tumor_event_type=="Locoregional Recurrence|Distant Metastasis" | 
                                                                         PAAD.repertoire.tumor.filter$new_tumor_event_type=="New Primary Tumor",NA)
PAAD.repertoire.tumor.filter$new_tumor_event_type<-factor(PAAD.repertoire.tumor.filter$new_tumor_event_type)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"new_tumor_event_type")

#treatment_outcome_first_course
PAAD.repertoire.tumor.filter$treatment_outcome_first_course<-replace(PAAD.repertoire.tumor.filter$treatment_outcome_first_course,
                                                                                 PAAD.repertoire.tumor.filter$treatment_outcome_first_course=="[Discrepancy]" |
                                                                                   PAAD.repertoire.tumor.filter$treatment_outcome_first_course=="[Not Applicable]" | 
                                                                                   PAAD.repertoire.tumor.filter$treatment_outcome_first_course=="[Not Available]" |
                                                                                   PAAD.repertoire.tumor.filter$treatment_outcome_first_course=="[Unknown]",NA)
PAAD.repertoire.tumor.filter$treatment_outcome_first_course<-factor(PAAD.repertoire.tumor.filter$treatment_outcome_first_course)
association.test.immuneRep(PAAD.repertoire.tumor.filter,"treatment_outcome_first_course")

#############################################
### 4. Survival Analysis  ####
##############################################
library(survival)
library(survminer)
library(survMisc)
library(metaviz)
##OS
#PAAD.repertoire.tumor.survival$number_pack_years_smoked<-as.numeric(as.character(PAAD.repertoire.tumor.survival$number_pack_years_smoked))
PAAD.repertoire.tumor.filter$IGH_expression_log<-log10(PAAD.repertoire.tumor.filter$IGH_expression)
PAAD.repertoire.tumor.filter$IGK_expression_log<-log10(PAAD.repertoire.tumor.filter$IGK_expression)
PAAD.repertoire.tumor.filter$IGL_expression_log<-log10(PAAD.repertoire.tumor.filter$IGL_expression)
PAAD.repertoire.tumor.filter$TRA_expression_log<-log10(PAAD.repertoire.tumor.filter$TRA_expression)
PAAD.repertoire.tumor.filter$TRB_expression_log<-log10(PAAD.repertoire.tumor.filter$TRB_expression)
surv_object <- Surv(time = PAAD.repertoire.tumor.filter$OS.time, event = PAAD.repertoire.tumor.filter$OS)
PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster<-factor(PAAD.repertoire.tumor.filter$IGK_clonotypes_cluster_PDAC)
PAAD.repertoire.tumor.filter.smoking<-PAAD.repertoire.tumor.filter[which(is.na(PAAD.repertoire.tumor.filter$number_pack_years_smoked_all)==F),]
res.cox <- coxph(Surv(time = OS.time, event = OS)~IGK_clonotypes_cluster+gender+ as.numeric(as.character(age_at_initial_pathologic_diagnosis))
                 +pathologic_stage,data=PAAD.repertoire.tumor.filter)
summary(res.cox)
ggforest(res.cox)

##Preparing plot
All_markers<-c("IGH_expression_log","IGK_expression_log","IGL_expression_log","TRA_expression_log","TRB_expression_log"
               ,"entropy_IGH", "entropy_IGK", "entropy_IGL",
               "entropy_TRA", "entropy_TRB","cluster_gini_IGH","vertex_gini_IGH",
               "cluster_gini_IGK","vertex_gini_IGK","cluster_gini_IGL","vertex_gini_IGL")
coef<-NULL
HR<-NULL
se<-NULL
p<-NULL
for(i in 1:length(All_markers)){
  res.cox <- coxph(Surv(time = OS.time, event = OS)~PAAD.repertoire.tumor.filter[,All_markers[i]]+gender
                   +as.numeric(as.character(age_at_initial_pathologic_diagnosis))+pathologic_stage+smoking
                   ,data=PAAD.repertoire.tumor.filter)
  coef[i]<-coef(summary(res.cox))[1,1]
  HR[i]<-coef(summary(res.cox))[1,2]
  se[i]<-coef(summary(res.cox))[1,3]
  p[i]<-coef(summary(res.cox))[1,5]
}

Survival_summary<-cbind(coef,HR,se,p)
rownames(Survival_summary)<-All_markers
tiff("Manuscript/Figures/SurvivalPlot.tiff",res=300,h=1700,w=2400)
viz_forest(x = Survival_summary[, c("coef", "se")], study_labels = rownames(Survival_summary),variant = "rain",
           type="study_only", x_trans_function = exp, xlab = "Hazard-Ratio",annotate_CI = TRUE)

dev.off()

##When restricted to number pack years smoked the sample size is too small and we loose significance
##Categorical
KL_mean<-mean(PAAD.repertoire.tumor.filter$entropy_IGH,na.rm=T)
PAAD.repertoire.tumor.filter$KL_ratio_2cat<-ifelse(as.numeric(PAAD.repertoire.tumor.filter$entropy_IGH)<=KL_mean,1,2)
fit1 <- survfit(surv_object ~ PAAD.repertoire.tumor.filter$smoking)
fit1
tiff("Manuscript/Figures/KM_BCR.tiff",res=300,h=1000,w=1000)
ggsurvplot(fit1, data = PAAD.repertoire.tumor.filter,legend="right",legend.title="BCR infiltration",legend.labs=c("low","medium","high"))
dev.off()
tiff("Results/Clinical/KM_tobacco.tiff",res=300,h=1000,w=1000)
ggsurvplot(fit1, data = PAAD.repertoire.tumor.filter,legend="right",legend.title="Smoking",legend.labs=c("Current","Former","Non-smoker"),pval = T)
dev.off()
comp(ten(fit1))$tests$lrTests


PAAD.repertoire.tumor.filter$number_pack_years_smoked<-as.numeric(as.character(PAAD.repertoire.tumor.filter$number_pack_years_smoked))
PAAD.repertoire.tumor.filter_smokers<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$number_pack_years_smoked>0),]

surv_object <- Surv(time = PAAD.repertoire.tumor.filter_smokers$OS.time, event = PAAD.repertoire.tumor.filter_smokers$OS)
res.cox <- coxph(Surv(time = OS.time, event = OS)~entropy_IGH,data=PAAD.repertoire.tumor.filter_smokers)
summary(res.cox)
ggforest(res.cox)

PAAD.repertoire.tumor.filter_smokers<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$smoking2=="Ever-Smoker"),]
surv_object <- Surv(time = PAAD.repertoire.tumor.filter_smokers$OS.time, event = PAAD.repertoire.tumor.filter_smokers$OS)
res.cox <- coxph(Surv(time = OS.time, event = OS)~entropy_IGH,data=PAAD.repertoire.tumor.filter_smokers)
summary(res.cox)
ggforest(res.cox)





############################################
####Compare tumor vs normal with its pair##
##########################################
normal_samples<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type_4categ=="normal_pancreas"),"TCGA_sample"]
summary_data<-NULL
for(i in 1:length(normal_samples)){
  id<-grep(substr(normal_samples[i],1,12),PAAD.repertoire.diversity$TCGA_sample)
  summary_data<-rbind(summary_data,PAAD.repertoire.diversity[id,c("IGH_expression","IGK_expression","IGL_expression",
                                                                  "entropy_IGH","entropy_IGK","entropy_IGL",
                                                                  "clones_IGH", "clones_IGK", "clones_IGL",
                                                                  "TRA_expression","TRB_expression","TRD_expression","TRG_expression",
                                                                  "entropy_TRA","entropy_TRB","entropy_TRD","entropy_TRG",
                                                                  "clones_TRA", "clones_TRB", "clones_TRD", "clones_TRG",
                                                                  "Tumor_type_4categ","TCGA_sample")])
}
summary_data<-summary_data[-3,]
marker<-c("IGH_expression","IGK_expression","IGL_expression",
          "entropy_IGH","entropy_IGK","entropy_IGL",
          "clones_IGH", "clones_IGK", "clones_IGL",
          "TRA_expression","TRB_expression","TRD_expression","TRG_expression",
          "entropy_TRA","entropy_TRB","entropy_TRD","entropy_TRG",
          "clones_TRA", "clones_TRB", "clones_TRD", "clones_TRG")
for(i in 1:length(marker)){
  summary_data$sample<-substr(summary_data$TCGA_sample,1,12)
  summary_data$Marker<-summary_data[,marker[i]]
  tiff(paste0("Results/ImmuneRep/Comparisons/tumor_adjacent_normal_",marker[i]))
  print(ggplot(summary_data, aes(x = Tumor_type_4categ, y = Marker)) + scale_y_continuous(name=marker[i]) + 
          scale_x_discrete("") +
    geom_boxplot() +
    geom_line(aes(group = sample)) +
    geom_point())
  dev.off()
}

for(i in 1:length(normal_samples)){
  #id_pair<-grep(substr(normal_samples[i],1,12),PAAD.repertoire.diversity$TCGA_sample)
  #PAAD.repertoire.diversity[id_pair,]
  print(ggline(summary_data, x = "Tumor_type_4categ", y = "IGH_expression", color = "Tumor_type_4categ",
               add = c("mean_se", "jitter"),  palette = c("#66C2A5", "#FC8D62")) + scale_y_continuous(name="IGH_expression") +
          stat_compare_means(aes(group = sample),method="anova",label = "p.format"))
 
   
}


####Ig/T expression only in PDAC samples
PAAD.repertoire.tumor.subtype.Ig<-PAAD.repertoire.tumor.subtype[which(PAAD.repertoire.tumor.subtype$Ig_Reads>1000),]
Ig_markers<-c("IGH_expression","IGK_expression", "IGL_expression","entropy_IGH","entropy_IGK", "entropy_IGL")

PAAD.repertoire.tumor.subtype.T<-PAAD.repertoire.tumor.subtype[which(PAAD.repertoire.tumor.subtype$T_Reads>100),]
T_markers<-c("TRA_expression","TRB_expression", "entropy_TRA","entropy_TRB")

##Correlation matrix
mat<-PAAD.repertoire.tumor.filter[,c(Ig_markers,T_markers)]
library(corrplot)
M <- cor(mat)
tiff("Results/ImmuneRep/Comparisons/Correlation.tiff",res=300,w=2200,h=2200)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
          sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

dev.off()

library("PerformanceAnalytics")
tiff("Results/ImmuneRep/Comparisons/Chat_Correlation.tiff",res=300,w=2300,h=2300)
chart.Correlation(mat, pch=19)
dev.off()

##### PCA plot
pca <- prcomp((mat), scale = TRUE)
pc <- c(1,2)
tiff(paste0("Results/PCA_clones_relative_abundance_IG.tiff"),width = 2000, height = 2000, res = 300)
plot(pca$x[,pc[1]], pca$x[,pc[2]], pch=20,xlab="PCA1",ylab="PCA2")
text(pca$x[,1], pca$x[,2], rownames(mat), pos= 2 )
dev.off()





#### Clustering ###
#res <- pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))
#mat.clust <- cbind(mat, cluster = cutree(res$tree_col, k = 4))
#mat.clust$cluster<-replace(mat.clust$cluster,mat.clust$cluster==4,3)

mat<-PAAD.repertoire.tumor.subtype.T[,T_markers]
annotation_col = data.frame(ABSOLUTE_Purity = 100*as.numeric(as.character(PAAD.repertoire.tumor.subtype.T$ABSOLUTE.Purity)),
                            Leukocyte_DNAmethylation = 100*as.numeric(as.character(PAAD.repertoire.tumor.subtype.T$Leukocyte_DNAmethylation)))
rownames(annotation_col)<-rownames(mat)

##Plot heatmap with cluster
ann_colors = list (ABSOLUTE_Purity = colorRampPalette(brewer.pal(9,name="Purples"))(120),
                   Leukocyte_DNAmethylation = colorRampPalette(brewer.pal(9,name="Greens"))(120))
tiff("Results/ImmuneRep/Comparisons/Cluster_TExpression_heatmap.tiff",res=300,w=2500,h=1000)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120),
         annotation_col = annotation_col,annotation_colors = ann_colors)

dev.off()




############################################ 
####Ig/T expression only in PDAC samples ##
##########################################
PAAD.repertoire.tumor.filter.Ig<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$Ig_Reads>1000),]
Ig_markers<-c("IGH_expression","IGK_expression", "IGL_expression")
Ig_markers<-c("entropy_IGH","entropy_IGK", "entropy_IGL")

PAAD.repertoire.tumor.filter.T<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$T_Reads>100),]
T_markers<-c("TRA_expression","TRB_expression", "TRD_expression", "TRG_expression")
T_markers<-c("entropy_TRA","entropy_TRB")

All_markers<-c("IGH_expression","IGK_expression", "IGL_expression","entropy_IGH","entropy_IGK", "entropy_IGL",
               "TRA_expression","TRB_expression", "TRD_expression", "TRG_expression","entropy_TRA","entropy_TRB")

##Correlation matrix
mat<-PAAD.repertoire.tumor.filter[,All_markers]
library(corrplot)
M <- cor(mat)
tiff("Results/ImmuneRep/Comparisons/Correlation_Allmarkers.tiff",res=300,w=2000,h=2000)
corrplot(M)
corrplot(M, order = "hclust", addrect = 4)
dev.off()

##### PCA plot
pca <- prcomp((mat), scale = TRUE)
pc <- c(1,2)
tiff(paste0("Results/PCA_clones_relative_abundance_IG.tiff"),width = 2000, height = 2000, res = 300)
plot(pca$x[,pc[1]], pca$x[,pc[2]], pch=20,xlab="PCA1",ylab="PCA2")
text(pca$x[,1], pca$x[,2], rownames(mat), pos= 2 )
dev.off()


#### Heatmaps ###
mat<-PAAD.repertoire.tumor.filter.Ig[,Ig_markers]
annotation_col = data.frame(ABSOLUTE_Purity = 100*as.numeric(as.character(PAAD.repertoire.tumor.filter.Ig$ABSOLUTE.Purity)),
                            Leukocyte_DNAmethylation = 100*as.numeric(as.character(PAAD.repertoire.tumor.filter.Ig$Leukocyte_DNAmethylation)))

mat<-PAAD.repertoire.tumor.filter.T[,T_markers]
annotation_col = data.frame(ABSOLUTE_Purity = 100*as.numeric(as.character(PAAD.repertoire.tumor.filter.T$ABSOLUTE.Purity)),
                            Leukocyte_DNAmethylation = 100*as.numeric(as.character(PAAD.repertoire.tumor.filter.T$Leukocyte_DNAmethylation)))

rownames(annotation_col)<-rownames(mat)

##Plot heatmap with cluster
#ann_colors = list (cluster = c("1" = brewer.pal(3,"Set1")[1], "2"= brewer.pal(3,"Set1")[3], "3" = brewer.pal(3,"Set1")[2]))
ann_colors = list (ABSOLUTE_Purity = colorRampPalette(brewer.pal(9,name="Purples"))(120),
                   Leukocyte_DNAmethylation = colorRampPalette(brewer.pal(9,name="Greens"))(120))
tiff("Results/ImmuneRep/Comparisons/Cluster_TExpression_heatmap.tiff",res=300,w=2500,h=1000)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120),
         annotation_col = annotation_col,annotation_colors = ann_colors)




#############################################################
#### Immune reperotire Analysis TCGA GTEx and Validation ###
############################################################

#### TCGA - Only PDAC samples
PAAD.repertoire.tumor.filter$TCGA_sample<-substr(PAAD.repertoire.tumor.filter$TCGA_sample,1,15)
#### GTEX - Pancreas
GTEX.repertoire.normal<-Pancreas.repertoire.diversity
##Tumor Validation
Validation.repertoire.tumor<-Pancreas.Validation.repertoire.diversity[which(Pancreas.Validation.repertoire.diversity$tissue=="pancreas tumor"),]
##Normal Validation
Validation.repertoire.normal<-Pancreas.Normal.Validation.repertoire.diversity

Repertoire.Diversity<-rbind(PAAD.repertoire.tumor.filter[,c("Ig_Reads","T_Reads","IGH_expression","IGK_expression","IGL_expression",
                                               "TRA_expression", "TRB_expression", "TRD_expression", "TRG_expression",
                                               "Alpha_Beta_ratio_expression", "KappaLambda_ratio_expression",
                                               "clones_IGH","clones_IGK","clones_IGL",
                                               "clones_TRA","clones_TRB","clones_TRD","clones_TRG",
                                               "entropy_IGH", "entropy_IGK", "entropy_IGL",
                                               "entropy_TRA", "entropy_TRB", "entropy_TRD", "entropy_TRG")],
                            GTEX.repertoire.normal[,c("Ig_Reads","T_Reads","IGH_expression","IGK_expression","IGL_expression",
                                                     "TRA_expression", "TRB_expression", "TRD_expression", "TRG_expression",
                                                     "Alpha_Beta_ratio_expression", "KappaLambda_ratio_expression",
                                                     "clones_IGH","clones_IGK","clones_IGL",
                                                     "clones_TRA","clones_TRB","clones_TRD","clones_TRG",
                                                     "entropy_IGH", "entropy_IGK", "entropy_IGL",
                                                     "entropy_TRA", "entropy_TRB", "entropy_TRD", "entropy_TRG")],
                            Validation.repertoire.tumor[,c("Ig_Reads","T_Reads","IGH_expression","IGK_expression","IGL_expression",
                                                     "TRA_expression", "TRB_expression", "TRD_expression", "TRG_expression",
                                                     "Alpha_Beta_ratio_expression", "KappaLambda_ratio_expression",
                                                     "clones_IGH","clones_IGK","clones_IGL",
                                                     "clones_TRA","clones_TRB","clones_TRD","clones_TRG",
                                                     "entropy_IGH", "entropy_IGK", "entropy_IGL",
                                                     "entropy_TRA", "entropy_TRB", "entropy_TRD", "entropy_TRG")],
                            Validation.repertoire.normal[,c("Ig_Reads","T_Reads","IGH_expression","IGK_expression","IGL_expression",
                                                     "TRA_expression", "TRB_expression", "TRD_expression", "TRG_expression",
                                                     "Alpha_Beta_ratio_expression", "KappaLambda_ratio_expression",
                                                     "clones_IGH","clones_IGK","clones_IGL",
                                                     "clones_TRA","clones_TRB","clones_TRD","clones_TRG",
                                                     "entropy_IGH", "entropy_IGK", "entropy_IGL",
                                                     "entropy_TRA", "entropy_TRB", "entropy_TRD", "entropy_TRG")])

Repertoire.Diversity$outcome<-c(rep("TCGA-PDAC",nrow(PAAD.repertoire.tumor.filter)),rep("GTEX-Normal",nrow(GTEX.repertoire.normal)),
                                rep("Validation-PDAC",nrow(Validation.repertoire.tumor)),rep("Validation-Normal",nrow(Validation.repertoire.normal)))
Repertoire.Diversity$outcome<-factor(Repertoire.Diversity$outcome)


## Heatmap for the BCR and TCR ###
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")
annotation_col = data.frame(Repertoire.Diversity$outcome)
ann_colors = list (outcome = c("GTEX-Normal" = cols[1],
                                "TCGA-PDAC" = cols[2],
                                "Validation-Normal"= cols[3],
                                "Validation-PDAC" = cols[4]))
colnames(annotation_col)<-"outcome"
rownames(annotation_col)<-rownames(Repertoire.Diversity)

#IG
Ig_markers<-c("IGH_expression","IGK_expression", "IGL_expression")
mat<-Repertoire.Diversity[,Ig_markers]
rownames(mat)<-rownames(Repertoire.Diversity)
tiff("Results/ImmuneRep/Comparisons/Ig_expression_heatmap.tiff",h=1500,w=2500,res=300)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(12))
dev.off()
#TR
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")
annotation_col = data.frame(Repertoire.Diversity$outcome)
ann_colors = list (outcome = c("GTEX-Normal" = cols[1],
                               "TCGA-PDAC" = cols[2],
                               "Validation-Normal"= cols[3],
                               "Validation-PDAC" = cols[4]))
colnames(annotation_col)<-"outcome"
rownames(annotation_col)<-rownames(Repertoire.Diversity)

T_markers<-c("TRA_expression","TRB_expression","TRD_expression","TRG_expression")
mat<-Repertoire.Diversity[,T_markers]
rownames(mat)<-rownames(Repertoire.Diversity)
tiff("Results/ImmuneRep/Comparisons//T_expression_heatmap.tiff",h=1500,w=2500,res=300)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(12),fontsize = 12)
dev.off()

#Ig entropy 
Repertoire.Diversity.filter<-Repertoire.Diversity[which(is.na(Repertoire.Diversity$entropy_IGH)==F),]
Ig_markers<-c("entropy_IGH", "entropy_IGK", "entropy_IGL")
mat<-Repertoire.Diversity.filter[,Ig_markers]
rownames(mat)<-rownames(Repertoire.Diversity.filter)
tiff("Results/ImmuneRep/Comparisons/Ig_entropy_heatmap.tiff",h=1500,w=2500,res=300)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))
dev.off()

#T entropy 
Repertoire.Diversity.filter<-Repertoire.Diversity[which(is.na(Repertoire.Diversity$entropy_TRA)==F),]
Ig_markers<-c("entropy_TRA", "entropy_TRB")
mat<-Repertoire.Diversity.filter[,Ig_markers]
rownames(mat)<-rownames(Repertoire.Diversity.filter)
tiff("Results/ImmuneRep/Comparisons/T_entropy_heatmap.tiff",h=1000,w=2100,res=300)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12),fontsize = 12)
dev.off()

######################
##Correlation matrix
######################
All_markers<-c("IGH_expression","IGK_expression", "IGL_expression","entropy_IGH", "entropy_IGK", "entropy_IGL",
               "TRA_expression","TRB_expression","entropy_TRA", "entropy_TRB")
mat_TCGA<-Repertoire.Diversity[which(Repertoire.Diversity$outcome=="TCGA-PDAC"),All_markers]
mat_GTEX<-Repertoire.Diversity[which(Repertoire.Diversity$outcome=="GTEX-Normal"),All_markers]
mat_PDAC_Validation<-Repertoire.Diversity[which(Repertoire.Diversity$outcome=="Validation-PDAC"),All_markers]
mat_normal_Validation<-Repertoire.Diversity[which(Repertoire.Diversity$outcome=="Validation-Normal"),All_markers]

library(corrplot)
tiff("Results/ImmuneRep/Comparisons/CorrelationTCGA.tiff",res=300,w=2200,h=2200)
M <- cor(mat_TCGA)
p <- cor.mtest(mat_TCGA, conf.level = 0.05)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(200),  
         type="upper", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat=p$p,vsig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)

dev.off()

library("PerformanceAnalytics")
tiff("Results/ImmuneRep/Comparisons/Chat_Correlation_GTEX.tiff",res=300,w=2300,h=2300)
chart.Correlation(mat_GTEX, pch=19)
dev.off()



###################
####Summary plots##
####################
Ig_expr<-melt(Repertoire.Diversity[,c("IGH_expression","IGK_expression","IGL_expression","TRA_expression",
                                      "TRB_expression","outcome")])
Ig_expr<-Ig_expr[which(Ig_expr$value!=0),]
Ig_expr$value<-log10(Ig_expr$value)
levels(Ig_expr$variable)<-c("IGH","IGK","IGL","TRA","TRB")
tiff("Results/ImmuneRep/Comparisons/boxplot_Ig_expression.tiff",res=300,h=1200,w=2500)
ggboxplot(Ig_expr, x = "outcome", y = "value",facet.by = "variable",nrow=1,color = "outcome",ggtheme = theme_bw(),xlab = F,
          outlier.size = 0) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  + ylab("Expression(log10)") +
  geom_point(aes(x=outcome, y=value, color=outcome),size=0.5, position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Normal",
                                                  "TCGA-PDAC",
                                                  "Validation-Normal",
                                                  "Validation-PDAC")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Normal","Validation-Normal"),c("GTEX-Normal","TCGA-PDAC"),
                      c("TCGA-PDAC","Validation-PDAC"),c("Validation-Normal","Validation-PDAC"),
                      c("TCGA-PDAC","Validation-Normal")))


dev.off()

Ig_entropy<-melt(Repertoire.Diversity[,c("entropy_IGH","entropy_IGK","entropy_IGL","entropy_TRA","entropy_TRB","outcome")])
Ig_entropy<-Ig_entropy[which(Ig_entropy$value!=0),]
levels(Ig_entropy$variable)<-c("IGH","IGK","IGL","TRA","TRB")
tiff("Results/ImmuneRep/Comparisons/boxplot_Ig_entropy.tiff",res=300,h=1200,w=2500)
ggboxplot(Ig_entropy, x = "outcome", y = "value",facet.by = "variable",nrow=1,color = "outcome",ggtheme = theme_bw(),xlab = F,
          outlier.size = 0)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylab("Shannon Entropy") +
  geom_point(aes(x=outcome, y=value, color=outcome),size=0.5, position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Normal",
                                                  "TCGA-PDAC",
                                                  "Validation-Normal",
                                                  "Validation-PDAC")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Normal","Validation-Normal"),c("GTEX-Normal","TCGA-PDAC"),
                      c("TCGA-PDAC","Validation-PDAC"),c("Validation-Normal","Validation-PDAC"),
                      c("TCGA-PDAC","Validation-Normal")))


dev.off()

kappa_lambda<-melt(Repertoire.Diversity[,c("KappaLambda_ratio_expression","outcome")])
tiff("Results/ImmuneRep/Comparisons/boxplot_kappa_lambda.tiff",res=300,h=2500,w=3000)
ggboxplot(kappa_lambda, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Normal",
                                                  "TCGA-PDAC",
                                                  "Validation-Normal",
                                                  "Validation-PDAC")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Normal","Validation-Normal"),c("GTEX-Normal","TCGA-PDAC"),
                      c("TCGA-PDAC","Validation-PDAC"),c("Validation-Normal","Validation-PDAC"),
                      c("TCGA-PDAC","Validation-Normal")))


dev.off()

alpha_beta_ratio<-melt(Repertoire.Diversity[,c("Alpha_Beta_ratio_expression","outcome")])
tiff("Results/ImmuneRep/Comparisons/boxplot_alpha_beta_ratio.tiff",res=300,h=2500,w=3000)
ggboxplot(alpha_beta_ratio, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Normal",
                                                  "TCGA-PDAC",
                                                  "Validation-Normal",
                                                  "Validation-PDAC")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Normal","Validation-Normal"),c("GTEX-Normal","TCGA-PDAC"),
                      c("TCGA-PDAC","Validation-PDAC"),c("Validation-Normal","Validation-PDAC"),
                      c("TCGA-PDAC","Validation-Normal")))


dev.off()

#############################################################
#### Immune reperotire Analysis TCGA vs. GTEx  ##############
############################################################

#### TCGA - Only PDAC samples
PAAD.repertoire.tumor<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type_4categ=="PDAC"),] 
PAAD.repertoire.tumor$TCGA_sample<-substr(PAAD.repertoire.tumor$TCGA_sample,1,15)
#### GTEX - Pancreas
GTEX.repertoire.normal<-Pancreas.repertoire.diversity

Repertoire.Diversity<-rbind(PAAD.repertoire.tumor[,c("Ig_Reads","T_Reads","IGH_expression","IGK_expression","IGL_expression",
                                                     "TRA_expression", "TRB_expression", "TRD_expression", "TRG_expression",
                                                     "Alpha_Beta_ratio_expression", "KappaLambda_ratio_expression",
                                                     "clones_IGH","clones_IGK","clones_IGL",
                                                     "clones_TRA","clones_TRB","clones_TRD","clones_TRG",
                                                     "entropy_IGH", "entropy_IGK", "entropy_IGL",
                                                     "entropy_TRA", "entropy_TRB", "entropy_TRD", "entropy_TRG")],
                            GTEX.repertoire.normal[,c("Ig_Reads","T_Reads","IGH_expression","IGK_expression","IGL_expression",
                                                      "TRA_expression", "TRB_expression", "TRD_expression", "TRG_expression",
                                                      "Alpha_Beta_ratio_expression", "KappaLambda_ratio_expression",
                                                      "clones_IGH","clones_IGK","clones_IGL",
                                                      "clones_TRA","clones_TRB","clones_TRD","clones_TRG",
                                                      "entropy_IGH", "entropy_IGK", "entropy_IGL",
                                                      "entropy_TRA", "entropy_TRB", "entropy_TRD", "entropy_TRG")])
                          
Repertoire.Diversity$outcome<-c(rep("TCGA-PDAC",nrow(PAAD.repertoire.tumor)),rep("GTEX-Normal",nrow(GTEX.repertoire.normal)))
Repertoire.Diversity$outcome<-factor(Repertoire.Diversity$outcome)


####Summary plots
## Heatmap for the BCR and TCR ###
cols= c("#7FC97F","#BEAED4","#FDC086","#B3CDE3")
annotation_col = data.frame(Repertoire.Diversity$outcome)
ann_colors = list (outcome = c("GTEX-Normal" = cols[1],
                               "TCGA-PDAC" = cols[2]))
colnames(annotation_col)<-"outcome"
rownames(annotation_col)<-rownames(Repertoire.Diversity)

Ig_TR_expr<-melt(Repertoire.Diversity[,c("TRA_expression","TRB_expression","TRD_expression","TRG_expression","IGH_expression","IGK_expression","IGL_expression","outcome")])
Ig_TR_expr<-Ig_TR_expr[which(Ig_TR_expr$value!=0),]
Ig_TR_expr$value<-log10(Ig_TR_expr$value)
tiff("Results/ImmuneRep/Comparisons/boxplot_Ig_TR_expression_GTEX_TCGA.tiff",res=300,h=2500,w=3000)
ggboxplot(Ig_TR_expr, x = "outcome", y = "value",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  facet_wrap(~variable,nrow=2) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Normal",
                                                  "TCGA-PDAC")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Normal","TCGA-PDAC")))


dev.off()

Ig_entropy<-melt(Repertoire.Diversity[,c("entropy_IGH","entropy_IGK","entropy_IGL","entropy_TRA","entropy_TRB","outcome")])
Ig_entropy<-Ig_entropy[which(Ig_entropy$value!=0),]
tiff("Results/ImmuneRep/Comparisons/boxplot_Ig_TR_entropy_GTEX_TCGA.tiff",res=300,h=2500,w=3000)
ggboxplot(Ig_entropy, x = "outcome", y = "value",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  facet_wrap(~variable) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Normal",
                                                  "TCGA-PDAC")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Normal","TCGA-PDAC")))
dev.off()


kappa_lambda<-melt(Repertoire.Diversity[,c("KappaLambda_ratio_expression","outcome")])
tiff("Results/ImmuneRep/Comparisons/boxplot_kappa_lambda_GTEX_TCGA.tiff",res=300,h=2500,w=3000)
ggboxplot(kappa_lambda, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Normal",
                                                  "TCGA-PDAC")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Normal","TCGA-PDAC")))


dev.off()

alpha_beta_ratio<-melt(Repertoire.Diversity[,c("Alpha_Beta_ratio_expression","outcome")])
tiff("Results/ImmuneRep/Comparisons/boxplot_alpha_beta_ratio_GTEX_TCGA.tiff",res=300,h=2500,w=3000)
ggboxplot(alpha_beta_ratio, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Normal",
                                                  "TCGA-PDAC")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Normal","TCGA-PDAC")))


dev.off()

#######Subsampled####
Repertoire.Diversity<-rbind(PAAD.repertoire.tumor[,c( "entropy_20_IGH","entropy_40_IGH","entropy_60_IGH","entropy_80_IGH","entropy_IGH",
                                                      "entropy_20_IGK","entropy_40_IGK","entropy_60_IGK","entropy_80_IGK","entropy_IGK",
                                                      "entropy_20_IGL","entropy_40_IGL","entropy_60_IGL","entropy_80_IGL","entropy_IGL",
                                                      "entropy_20_TRA","entropy_40_TRA","entropy_60_TRA","entropy_80_TRA","entropy_TRA",
                                                      "entropy_20_TRB","entropy_40_TRB","entropy_60_TRB","entropy_80_TRB","entropy_TRB",
                                                      "entropy_20_TRD","entropy_40_TRD","entropy_60_TRD","entropy_80_TRD","entropy_TRD",
                                                      "entropy_20_TRG","entropy_40_TRG","entropy_60_TRG","entropy_80_TRG","entropy_TRG")],
                                                      
                            GTEX.repertoire.normal[,c( "entropy_20_IGH","entropy_40_IGH","entropy_60_IGH","entropy_80_IGH","entropy_IGH",
                                                       "entropy_20_IGK","entropy_40_IGK","entropy_60_IGK","entropy_80_IGK","entropy_IGK",
                                                       "entropy_20_IGL","entropy_40_IGL","entropy_60_IGL","entropy_80_IGL","entropy_IGL",
                                                       "entropy_20_TRA","entropy_40_TRA","entropy_60_TRA","entropy_80_TRA","entropy_TRA",
                                                       "entropy_20_TRB","entropy_40_TRB","entropy_60_TRB","entropy_80_TRB","entropy_TRB",
                                                       "entropy_20_TRD","entropy_40_TRD","entropy_60_TRD","entropy_80_TRD","entropy_TRD",
                                                       "entropy_20_TRG","entropy_40_TRG","entropy_60_TRG","entropy_80_TRG","entropy_TRG")])
                            
###Plot the results
#Entropy IGH
Entropy_PAAD<-melt(PAAD.repertoire.tumor[,c("TCGA_sample","entropy_20_IGH","entropy_40_IGH","entropy_60_IGH","entropy_80_IGH","entropy_IGH")])
Entropy_PAAD$percentage<-c(rep(20,dim(PAAD.repertoire.tumor)[1]),rep(40,dim(PAAD.repertoire.tumor)[1]),rep(60,dim(PAAD.repertoire.tumor)[1]),
                                rep(80,dim(PAAD.repertoire.tumor)[1]),rep(100,dim(PAAD.repertoire.tumor)[1]))
colnames(Entropy_PAAD)[1]<-"sample"

Entropy_GTEX<-melt(GTEX.repertoire.normal[,c("SUBJID","entropy_20_IGH","entropy_40_IGH","entropy_60_IGH","entropy_80_IGH","entropy_IGH")])
Entropy_GTEX$percentage<-c(rep(20,dim(GTEX.repertoire.normal)[1]),rep(40,dim(GTEX.repertoire.normal)[1]),rep(60,dim(GTEX.repertoire.normal)[1]),
                                rep(80,dim(GTEX.repertoire.normal)[1]),rep(100,dim(GTEX.repertoire.normal)[1]))

colnames(Entropy_GTEX)[1]<-"sample"

Entropy<-rbind(Entropy_PAAD,Entropy_GTEX)

tiff("Results/ImmuneRep/Comparisons/Entropy_IGH_subsampled.tiff",res=300,h=2000,w=2000)
ggplot(data = Entropy, aes(x=percentage,y=value,colour=sample)) + geom_line() + geom_point() +
  theme(legend.position = "none") + scale_color_manual(values=c(rep("#7FC97F",dim(GTEX.repertoire.normal)[1]),rep("#BEAED4",dim(PAAD.repertoire.tumor)[1]))) +
  ylab("Entropy IGH") + xlab("proportion of reads samples")
dev.off()

#Entropy IGK
Entropy_PAAD<-melt(PAAD.repertoire.tumor[,c("TCGA_sample","entropy_20_IGK","entropy_40_IGK","entropy_60_IGK","entropy_80_IGK","entropy_IGK")])
Entropy_PAAD$percentage<-c(rep(20,dim(PAAD.repertoire.tumor)[1]),rep(40,dim(PAAD.repertoire.tumor)[1]),rep(60,dim(PAAD.repertoire.tumor)[1]),
                           rep(80,dim(PAAD.repertoire.tumor)[1]),rep(100,dim(PAAD.repertoire.tumor)[1]))
colnames(Entropy_PAAD)[1]<-"sample"

Entropy_GTEX<-melt(GTEX.repertoire.normal[,c("SUBJID","entropy_20_IGK","entropy_40_IGK","entropy_60_IGK","entropy_80_IGK","entropy_IGK")])
Entropy_GTEX$percentage<-c(rep(20,dim(GTEX.repertoire.normal)[1]),rep(40,dim(GTEX.repertoire.normal)[1]),rep(60,dim(GTEX.repertoire.normal)[1]),
                           rep(80,dim(GTEX.repertoire.normal)[1]),rep(100,dim(GTEX.repertoire.normal)[1]))

colnames(Entropy_GTEX)[1]<-"sample"

Entropy<-rbind(Entropy_PAAD,Entropy_GTEX)

tiff("Results/ImmuneRep/Comparisons/Entropy_IGK_subsampled.tiff",res=300,h=2000,w=2000)
ggplot(data = Entropy, aes(x=percentage,y=value,colour=sample)) + geom_line() + geom_point() +
  theme(legend.position = "none") + scale_color_manual(values=c(rep("#7FC97F",dim(GTEX.repertoire.normal)[1]),rep("#BEAED4",dim(PAAD.repertoire.tumor)[1]))) +
  ylab("Entropy IGK") + xlab("proportion of reads samples")
dev.off()

#Entropy IGL
Entropy_PAAD<-melt(PAAD.repertoire.tumor[,c("TCGA_sample","entropy_20_IGL","entropy_40_IGL","entropy_60_IGL","entropy_80_IGL","entropy_IGL")])
Entropy_PAAD$percentage<-c(rep(20,dim(PAAD.repertoire.tumor)[1]),rep(40,dim(PAAD.repertoire.tumor)[1]),rep(60,dim(PAAD.repertoire.tumor)[1]),
                           rep(80,dim(PAAD.repertoire.tumor)[1]),rep(100,dim(PAAD.repertoire.tumor)[1]))
colnames(Entropy_PAAD)[1]<-"sample"

Entropy_GTEX<-melt(GTEX.repertoire.normal[,c("SUBJID","entropy_20_IGL","entropy_40_IGL","entropy_60_IGL","entropy_80_IGL","entropy_IGL")])
Entropy_GTEX$percentage<-c(rep(20,dim(GTEX.repertoire.normal)[1]),rep(40,dim(GTEX.repertoire.normal)[1]),rep(60,dim(GTEX.repertoire.normal)[1]),
                           rep(80,dim(GTEX.repertoire.normal)[1]),rep(100,dim(GTEX.repertoire.normal)[1]))

colnames(Entropy_GTEX)[1]<-"sample"

Entropy<-rbind(Entropy_PAAD,Entropy_GTEX)

tiff("Results/ImmuneRep/Comparisons/Entropy_IGL_subsampled.tiff",res=300,h=2000,w=2000)
ggplot(data = Entropy, aes(x=percentage,y=value,colour=sample)) + geom_line() + geom_point() +
  theme(legend.position = "none") + scale_color_manual(values=c(rep("#7FC97F",dim(GTEX.repertoire.normal)[1]),rep("#BEAED4",dim(PAAD.repertoire.tumor)[1]))) +
  ylab("Entropy IGL") + xlab("proportion of reads samples")
dev.off()

#Entropy TRA
Entropy_PAAD<-melt(PAAD.repertoire.tumor[,c("TCGA_sample","entropy_20_TRA","entropy_40_TRA","entropy_60_TRA","entropy_80_TRA","entropy_TRA")])
Entropy_PAAD$percentage<-c(rep(20,dim(PAAD.repertoire.tumor)[1]),rep(40,dim(PAAD.repertoire.tumor)[1]),rep(60,dim(PAAD.repertoire.tumor)[1]),
                           rep(80,dim(PAAD.repertoire.tumor)[1]),rep(100,dim(PAAD.repertoire.tumor)[1]))
colnames(Entropy_PAAD)[1]<-"sample"

Entropy_GTEX<-melt(GTEX.repertoire.normal[,c("SUBJID","entropy_20_TRA","entropy_40_TRA","entropy_60_TRA","entropy_80_TRA","entropy_TRA")])
Entropy_GTEX$percentage<-c(rep(20,dim(GTEX.repertoire.normal)[1]),rep(40,dim(GTEX.repertoire.normal)[1]),rep(60,dim(GTEX.repertoire.normal)[1]),
                           rep(80,dim(GTEX.repertoire.normal)[1]),rep(100,dim(GTEX.repertoire.normal)[1]))

colnames(Entropy_GTEX)[1]<-"sample"

Entropy<-rbind(Entropy_PAAD,Entropy_GTEX)

tiff("Results/ImmuneRep/Comparisons/Entropy_TRA_subsampled.tiff",res=300,h=2000,w=2000)
ggplot(data = Entropy, aes(x=percentage,y=value,colour=sample)) + geom_line() + geom_point() +
  theme(legend.position = "none") + scale_color_manual(values=c(rep("#7FC97F",dim(GTEX.repertoire.normal)[1]),rep("#BEAED4",dim(PAAD.repertoire.tumor)[1]))) +
  ylab("Entropy TRA") + xlab("proportion of reads samples")
dev.off()

#Entropy TRB
Entropy_PAAD<-melt(PAAD.repertoire.tumor[,c("TCGA_sample","entropy_20_TRB","entropy_40_TRB","entropy_60_TRB","entropy_80_TRB","entropy_TRB")])
Entropy_PAAD$percentage<-c(rep(20,dim(PAAD.repertoire.tumor)[1]),rep(40,dim(PAAD.repertoire.tumor)[1]),rep(60,dim(PAAD.repertoire.tumor)[1]),
                           rep(80,dim(PAAD.repertoire.tumor)[1]),rep(100,dim(PAAD.repertoire.tumor)[1]))
colnames(Entropy_PAAD)[1]<-"sample"

Entropy_GTEX<-melt(GTEX.repertoire.normal[,c("SUBJID","entropy_20_TRB","entropy_40_TRB","entropy_60_TRB","entropy_80_TRB","entropy_TRB")])
Entropy_GTEX$percentage<-c(rep(20,dim(GTEX.repertoire.normal)[1]),rep(40,dim(GTEX.repertoire.normal)[1]),rep(60,dim(GTEX.repertoire.normal)[1]),
                           rep(80,dim(GTEX.repertoire.normal)[1]),rep(100,dim(GTEX.repertoire.normal)[1]))

colnames(Entropy_GTEX)[1]<-"sample"

Entropy<-rbind(Entropy_PAAD,Entropy_GTEX)

tiff("Results/ImmuneRep/Comparisons/Entropy_TRB_subsampled.tiff",res=300,h=2000,w=2000)
ggplot(data = Entropy, aes(x=percentage,y=value,colour=sample)) + geom_line() + geom_point() +
  theme(legend.position = "none") + scale_color_manual(values=c(rep("#7FC97F",dim(GTEX.repertoire.normal)[1]),rep("#BEAED4",dim(PAAD.repertoire.tumor)[1]))) +
  ylab("Entropy TRB") + xlab("proportion of reads samples")
dev.off()

#Entropy TRD
Entropy_PAAD<-melt(PAAD.repertoire.tumor[,c("TCGA_sample","entropy_20_TRD","entropy_40_TRD","entropy_60_TRD","entropy_80_TRD","entropy_TRD")])
Entropy_PAAD$percentage<-c(rep(20,dim(PAAD.repertoire.tumor)[1]),rep(40,dim(PAAD.repertoire.tumor)[1]),rep(60,dim(PAAD.repertoire.tumor)[1]),
                           rep(80,dim(PAAD.repertoire.tumor)[1]),rep(100,dim(PAAD.repertoire.tumor)[1]))
colnames(Entropy_PAAD)[1]<-"sample"

Entropy_GTEX<-melt(GTEX.repertoire.normal[,c("SUBJID","entropy_20_TRD","entropy_40_TRD","entropy_60_TRD","entropy_80_TRD","entropy_TRD")])
Entropy_GTEX$percentage<-c(rep(20,dim(GTEX.repertoire.normal)[1]),rep(40,dim(GTEX.repertoire.normal)[1]),rep(60,dim(GTEX.repertoire.normal)[1]),
                           rep(80,dim(GTEX.repertoire.normal)[1]),rep(100,dim(GTEX.repertoire.normal)[1]))

colnames(Entropy_GTEX)[1]<-"sample"

Entropy<-rbind(Entropy_PAAD,Entropy_GTEX)

tiff("Results/ImmuneRep/Comparisons/Entropy_TRD_subsampled.tiff",res=300,h=2000,w=2000)
ggplot(data = Entropy, aes(x=percentage,y=value,colour=sample)) + geom_line() + geom_point() +
  theme(legend.position = "none") + scale_color_manual(values=c(rep("#7FC97F",dim(GTEX.repertoire.normal)[1]),rep("#BEAED4",dim(PAAD.repertoire.tumor)[1]))) +
  ylab("Entropy TRD") + xlab("proportion of reads samples")
dev.off()

#Entropy TRG
Entropy_PAAD<-melt(PAAD.repertoire.tumor[,c("TCGA_sample","entropy_20_TRG","entropy_40_TRG","entropy_60_TRG","entropy_80_TRG","entropy_TRG")])
Entropy_PAAD$percentage<-c(rep(20,dim(PAAD.repertoire.tumor)[1]),rep(40,dim(PAAD.repertoire.tumor)[1]),rep(60,dim(PAAD.repertoire.tumor)[1]),
                           rep(80,dim(PAAD.repertoire.tumor)[1]),rep(100,dim(PAAD.repertoire.tumor)[1]))
colnames(Entropy_PAAD)[1]<-"sample"

Entropy_GTEX<-melt(GTEX.repertoire.normal[,c("SUBJID","entropy_20_TRG","entropy_40_TRG","entropy_60_TRG","entropy_80_TRG","entropy_TRG")])
Entropy_GTEX$percentage<-c(rep(20,dim(GTEX.repertoire.normal)[1]),rep(40,dim(GTEX.repertoire.normal)[1]),rep(60,dim(GTEX.repertoire.normal)[1]),
                           rep(80,dim(GTEX.repertoire.normal)[1]),rep(100,dim(GTEX.repertoire.normal)[1]))

colnames(Entropy_GTEX)[1]<-"sample"

Entropy<-rbind(Entropy_PAAD,Entropy_GTEX)

tiff("Results/ImmuneRep/Comparisons/Entropy_TRG_subsampled.tiff",res=300,h=2000,w=2000)
ggplot(data = Entropy, aes(x=percentage,y=value,colour=sample)) + geom_line() + geom_point() +
  theme(legend.position = "none") + scale_color_manual(values=c(rep("#7FC97F",dim(GTEX.repertoire.normal)[1]),rep("#BEAED4",dim(PAAD.repertoire.tumor)[1]))) +
  ylab("Entropy TRG") + xlab("proportion of reads samples")
dev.off()

##################################################################
#### Immune reperotire Analysis TCGA GTEx tumor and GTEx Blood ###
#################################################################

#### TCGA - 
PAAD.repertoire<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type_4categ!="PAC-Other"),] #131
PAAD.repertoire$Tumor_type_4categ<-factor(PAAD.repertoire$Tumor_type_4categ)
PAAD.repertoire$TCGA_sample<-substr(PAAD.repertoire$TCGA_sample,1,15)
#### GTEX - Pancreas
GTEX.repertoire.normal<-Pancreas.repertoire.diversity
##GTEx - Blood
load("Data/GTEx/Blood/GTEx_FullData_OnlyDiversity.Rdata")
GTEX.blood.repertoire.diversity<-GTEX.blood.repertoire.diversity[which(is.na(GTEX.blood.repertoire.diversity$totalReads)==F),]

Repertoire.Diversity<-rbind(PAAD.repertoire[,c("Ig_Reads","T_Reads","IGH_expression","IGK_expression","IGL_expression",
                                                     "TRA_expression", "TRB_expression", "TRD_expression", "TRG_expression",
                                                     "Alpha_Beta_ratio_expression", "KappaLambda_ratio_expression",
                                                     "entropy_IGH", "entropy_IGK", "entropy_IGL",
                                                     "entropy_TRA", "entropy_TRB", "entropy_TRD", "entropy_TRG")],
                            GTEX.repertoire.normal[,c("Ig_Reads","T_Reads","IGH_expression","IGK_expression","IGL_expression",
                                                      "TRA_expression", "TRB_expression", "TRD_expression", "TRG_expression",
                                                      "Alpha_Beta_ratio_expression", "KappaLambda_ratio_expression",
                                                      "entropy_IGH", "entropy_IGK", "entropy_IGL",
                                                      "entropy_TRA", "entropy_TRB", "entropy_TRD", "entropy_TRG")],
                            GTEX.blood.repertoire.diversity[,c("Ig_Reads","T_Reads","IGH_expression","IGK_expression","IGL_expression",
                                                           "TRA_expression", "TRB_expression", "TRD_expression", "TRG_expression",
                                                           "Alpha_Beta_ratio_expression", "KappaLambda_ratio_expression",
                                                           "entropy_IGH", "entropy_IGK", "entropy_IGL",
                                                           "entropy_TRA", "entropy_TRB", "entropy_TRD", "entropy_TRG")])
                            
Repertoire.Diversity$outcome<-c(as.character(PAAD.repertoire$Tumor_type_4categ),rep("GTEX-Normal",nrow(GTEX.repertoire.normal)),
                                rep("GTEX-Blood",nrow(GTEX.blood.repertoire.diversity)))
Repertoire.Diversity$outcome<-replace(Repertoire.Diversity$outcome,Repertoire.Diversity$outcome=="normal_pancreas","TCGA-normal-adj-pancreas")
Repertoire.Diversity$outcome<-replace(Repertoire.Diversity$outcome,Repertoire.Diversity$outcome=="PDAC","TCGA-PDAC")
Repertoire.Diversity$outcome<-replace(Repertoire.Diversity$outcome,Repertoire.Diversity$outcome=="pseudonormal_pancreas","TCGA-pseudonormal_pancreas")

Repertoire.Diversity$outcome<-factor(Repertoire.Diversity$outcome)


## Heatmap for the BCR and TCR ##
cols= c("#FBB4AE","#7FC97F","#FDC086","#BEAED4","#B3CDE3")
annotation_col = data.frame(Repertoire.Diversity$outcome)
ann_colors = list (outcome = c("GTEX-Blood" = cols[1],
                               "GTEX-Normal" = cols[2],
                               "TCGA-normal-adj-pancreas" = cols[3],
                               "TCGA-PDAC" = cols[4],
                               "TCGA-pseudonormal_pancreas"= cols[5]))
colnames(annotation_col)<-"outcome"
rownames(annotation_col)<-rownames(Repertoire.Diversity)

#IG
Ig_markers<-c("IGH_expression","IGK_expression", "IGL_expression")
mat<-Repertoire.Diversity[,Ig_markers]
rownames(mat)<-rownames(Repertoire.Diversity)
mat<-mat[which(is.na(mat$IGH_expression)==F),]

tiff("Results/ImmuneRep/Comparisons/Ig_expression_ALL_blood_heatmap.tiff",h=2000,w=4000,res=300)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))
dev.off()

#TR
T_markers<-c("TRA_expression","TRB_expression","TRD_expression","TRG_expression")
mat<-Repertoire.Diversity[,T_markers]
rownames(mat)<-rownames(Repertoire.Diversity)
mat<-mat[which(is.na(mat$TRB_expression)==F),]
tiff("Results/ImmuneRep/Comparisons//T_expression_ALL_blood_heatmap.tiff",h=2000,w=4000,res=300)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))
dev.off()

#Ig entropy 
Repertoire.Diversity.filter<-Repertoire.Diversity[which(is.na(Repertoire.Diversity$entropy_IGH)==F),]
Ig_markers<-c("entropy_IGH", "entropy_IGK", "entropy_IGL")
mat<-Repertoire.Diversity.filter[,Ig_markers]
rownames(mat)<-rownames(Repertoire.Diversity.filter)
tiff("Results/ImmuneRep/Comparisons/Ig_entropy_ALL_blood_heatmap.tiff",h=2000,w=4000,res=300)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))
dev.off()

#T entropy 
Repertoire.Diversity.filter<-Repertoire.Diversity[which(is.na(Repertoire.Diversity$entropy_TRA)==F),]
Ig_markers<-c("entropy_TRA", "entropy_TRB")
mat<-Repertoire.Diversity.filter[,Ig_markers]
rownames(mat)<-rownames(Repertoire.Diversity.filter)
tiff("Results/ImmuneRep/Comparisons/TR_entropy_ALL_blood_heatmap.tiff",h=2000,w=4000,res=300)
pheatmap(t(mat),scale="row",show_colnames = F,border_color=F,annotation_col = annotation_col,
         annotation_colors = ann_colors,color = colorRampPalette(brewer.pal(6,name="PuOr"))(12))
dev.off()


####Summary plots
Ig_expr<-melt(Repertoire.Diversity[,c("IGH_expression","IGK_expression","IGL_expression","outcome")])
Ig_expr$value<-log10(Ig_expr$value)
tiff("Results/ImmuneRep/Comparisons/boxplot_Ig_expression_ALL_blood.tiff",res=300,h=2500,w=3000)
ggboxplot(Ig_expr, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Blood",
                                                  "GTEX-Normal",
                                                  "TCGA-normal-adj-pancreas",
                                                  "TCGA-PDAC",
                                                  "TCGA-pseudonormal_pancreas")) +
  stat_compare_means(
   comparisons =list(c("GTEX-Blood","GTEX-Normal"),c("GTEX-Blood","TCGA-normal-adj-pancreas"),
                     c("GTEX-Blood","TCGA-PDAC"),c("GTEX-Blood","TCGA-pseudonormal_pancreas"),
                    c("TCGA-PDAC","TCGA-normal-adj-pancreas"),c("TCGA-PDAC","TCGA-pseudonormal_pancreas"),
                    c("TCGA-normal-adj-pancreas","TCGA-pseudonormal_pancreas")))


dev.off()

TR_expr<-melt(Repertoire.Diversity[,c("TRA_expression","TRB_expression","TRD_expression","TRG_expression","outcome")])
TR_expr$value<-log10(TR_expr$value)
tiff("Results/ImmuneRep/Comparisons/boxplot_TR_expression_All_blood.tiff",res=300,h=2500,w=3500)
ggboxplot(TR_expr, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Blood",
                                                  "GTEX-Normal",
                                                  "TCGA-normal-adj-pancreas",
                                                  "TCGA-PDAC",
                                                  "TCGA-pseudonormal_pancreas")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Blood","GTEX-Normal"),c("GTEX-Blood","TCGA-normal-adj-pancreas"),
                      c("GTEX-Blood","TCGA-PDAC"),c("GTEX-Blood","TCGA-pseudonormal_pancreas"),
                      c("TCGA-PDAC","TCGA-normal-adj-pancreas"),c("TCGA-PDAC","TCGA-pseudonormal_pancreas"),
                      c("TCGA-normal-adj-pancreas","TCGA-pseudonormal_pancreas")))


dev.off()

Ig_entropy<-melt(Repertoire.Diversity[,c("entropy_IGH","entropy_IGK","entropy_IGL","outcome")])
Ig_entropy<-Ig_entropy[which(Ig_entropy$value!=0),]
tiff("Results/ImmuneRep/Comparisons/boxplot_Ig_entropy_ALL_blood.tiff",res=300,h=2500,w=3000)
ggboxplot(Ig_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Blood",
                                                  "GTEX-Normal",
                                                  "TCGA-normal-adj-pancreas",
                                                  "TCGA-PDAC",
                                                  "TCGA-pseudonormal_pancreas")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Blood","GTEX-Normal"),c("GTEX-Blood","TCGA-normal-adj-pancreas"),
                      c("GTEX-Blood","TCGA-PDAC"),c("GTEX-Blood","TCGA-pseudonormal_pancreas"),
                      c("TCGA-PDAC","TCGA-normal-adj-pancreas"),c("TCGA-PDAC","TCGA-pseudonormal_pancreas"),
                      c("TCGA-normal-adj-pancreas","TCGA-pseudonormal_pancreas")))


dev.off()

TR_entropy<-melt(Repertoire.Diversity[,c("entropy_TRA","entropy_TRB","outcome")])
TR_entropy<-TR_entropy[which(TR_entropy$value!=0),]
tiff("Results/ImmuneRep/Comparisons/boxplot_TR_entropy_blood_ALL.tiff",res=300,h=2500,w=3000)
ggboxplot(TR_entropy, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Blood",
                                                  "GTEX-Normal",
                                                  "TCGA-normal-adj-pancreas",
                                                  "TCGA-PDAC",
                                                  "TCGA-pseudonormal_pancreas")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Blood","GTEX-Normal"),c("GTEX-Blood","TCGA-normal-adj-pancreas"),
                      c("GTEX-Blood","TCGA-PDAC"),c("GTEX-Blood","TCGA-pseudonormal_pancreas"),
                      c("TCGA-PDAC","TCGA-normal-adj-pancreas"),c("TCGA-PDAC","TCGA-pseudonormal_pancreas"),
                      c("TCGA-normal-adj-pancreas","TCGA-pseudonormal_pancreas")))


dev.off()

kappa_lambda<-melt(Repertoire.Diversity[,c("KappaLambda_ratio_expression","outcome")])
kappa_lambda<-kappa_lambda[which(kappa_lambda$value<10),]
tiff("Results/ImmuneRep/Comparisons/boxplot_kappa_lambda_All_blood.tiff",res=300,h=2500,w=3000)
ggboxplot(kappa_lambda, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Blood",
                                                  "GTEX-Normal",
                                                  "TCGA-normal-adj-pancreas",
                                                  "TCGA-PDAC",
                                                  "TCGA-pseudonormal_pancreas")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Blood","GTEX-Normal"),c("GTEX-Blood","TCGA-normal-adj-pancreas"),
                      c("GTEX-Blood","TCGA-PDAC"),c("GTEX-Blood","TCGA-pseudonormal_pancreas"),
                      c("TCGA-PDAC","TCGA-normal-adj-pancreas"),c("TCGA-PDAC","TCGA-pseudonormal_pancreas"),
                      c("TCGA-normal-adj-pancreas","TCGA-pseudonormal_pancreas")))


dev.off()

alpha_beta_ratio<-melt(Repertoire.Diversity[,c("Alpha_Beta_ratio_expression","outcome")])
tiff("Results/ImmuneRep/Comparisons/boxplot_alpha_beta_ratio_ALL_blood.tiff",res=300,h=2500,w=3000)
ggboxplot(alpha_beta_ratio, x = "outcome", y = "value",facet.by = "variable",color = "outcome",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=outcome, y=value, color=outcome), position = position_jitterdodge(dodge.width = 0.8)) +
  scale_color_manual(values = c(cols), labels = c("GTEX-Blood",
                                                  "GTEX-Normal",
                                                  "TCGA-normal-adj-pancreas",
                                                  "TCGA-PDAC",
                                                  "TCGA-pseudonormal_pancreas")) +
  stat_compare_means(
    comparisons =list(c("GTEX-Blood","GTEX-Normal"),c("GTEX-Blood","TCGA-normal-adj-pancreas"),
                      c("GTEX-Blood","TCGA-PDAC"),c("GTEX-Blood","TCGA-pseudonormal_pancreas"),
                      c("TCGA-PDAC","TCGA-normal-adj-pancreas"),c("TCGA-PDAC","TCGA-pseudonormal_pancreas"),
                      c("TCGA-normal-adj-pancreas","TCGA-pseudonormal_pancreas")))
dev.off()



####Summary plots for PDAC
Repertoire.Diversity_PDAC<-Repertoire.Diversity[which(Repertoire.Diversity$outcome=="TCGA-PDAC"),]
Ig_T_expr<-melt(Repertoire.Diversity_PDAC[,c("IGH_expression","IGK_expression","IGL_expression","TRA_expression","TRB_expression","TRD_expression","TRG_expression")])
Ig_T_expr$value<-log10(Ig_T_expr$value)
tiff("Results/ImmuneRep/Comparisons/boxplot_Ig_TR_expression_differences.tiff",res=300,h=2500,w=3000)
ggboxplot(Ig_T_expr, x = "variable", y = "value",color = "variable",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=variable, y=value, color=variable), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(
    comparisons =list(c("IGH_expression","TRA_expression"),c("IGH_expression","TRB_expression"),
                      c("IGH_expression","TRD_expression"),c("IGH_expression","TRG_expression"),
                      c("IGK_expression","TRA_expression"),c("IGK_expression","TRB_expression"),
                      c("IGK_expression","TRD_expression"),c("IGK_expression","TRG_expression"),
                      c("IGL_expression","TRA_expression"),c("IGL_expression","TRB_expression"),
                      c("IGL_expression","TRD_expression"),c("IGL_expression","TRG_expression")))
dev.off()

####Summary plots for BLOOD
Repertoire.Diversity.Blood<-Repertoire.Diversity[which(Repertoire.Diversity$outcome=="GTEX-Blood"),]
Ig_T_expr<-melt(Repertoire.Diversity.Blood[,c("entropy_IGH","entropy_IGK","entropy_IGL","entropy_TRA","entropy_TRB","entropy_TRD","entropy_TRG")])
Ig_T_expr$value<-log10(Ig_T_expr$value)
tiff("Results/ImmuneRep/Comparisons/boxplot_Ig_TR_expression_differences_Blood.tiff",res=300,h=2500,w=3000)
ggboxplot(Ig_T_expr, x = "variable", y = "value",color = "variable",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=variable, y=value, color=variable), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(
    comparisons =list(c("entropy_IGH","entropy_TRA"),c("entropy_IGH","entropy_TRB"),
                      c("entropy_IGH","entropy_TRD"),c("entropy_IGH","entropy_TRG"),
                      c("entropy_IGK","entropy_TRA"),c("entropy_IGK","entropy_TRB"),
                      c("entropy_IGK","entropy_TRD"),c("entropy_IGK","entropy_TRG"),
                      c("entropy_IGL","entropy_TRA"),c("entropy_IGL","entropy_TRB"),
                      c("entropy_IGL","entropy_TRD"),c("entropy_IGL","entropy_TRG")))
dev.off()


#############################################
## Association analysis with mutated genes ##
############################################
PAAD.repertoire.tumor.filter.Ig<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$Ig_Reads>1000),]
PAAD.repertoire.mutated.genes<-PAAD.repertoire.tumor.filter.Ig[,137:670]
PAAD.repertoire.mutated.genes<-PAAD.repertoire.mutated.genes[,which(colSums(PAAD.repertoire.mutated.genes)>0)]

PAAD.repertoire.tumor.filter.T<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$T_Reads>100),]
PAAD.repertoire.mutated.genes<-PAAD.repertoire.tumor.filter.T[,137:670]
PAAD.repertoire.mutated.genes<-PAAD.repertoire.mutated.genes[,which(colSums(PAAD.repertoire.mutated.genes)>0)]

p.value<-NULL
for(i in 1:ncol(PAAD.repertoire.mutated.genes)){
  if(sum(PAAD.repertoire.mutated.genes[,i])>4){
    p.value[i]<-wilcox.test(PAAD.repertoire.tumor.filter.T$entropy_TRB~PAAD.repertoire.mutated.genes[,i])$p.value
  } else {
    p.value[i]<-NA
  }
}

PAAD.repertoire.mutated.genes[,which(p.value<0.05)]
boxplot(PAAD.repertoire.tumor.filter.T$entropy_TRB~PAAD.repertoire.mutated.genes[,which(p.value<0.05)])


p.value[which(p.value<0.05)]
colnames(PAAD.repertoire.mutated.genes)[which(p.value<0.05)]
table(PAAD.repertoire.mutated.genes[which(p.value<0.05)])




#################
## Function for the association analysis
###################
#Function to run the association between clinical outcome and BCR/TCR
association.test.immuneRep<- function (PAAD.repertoire.tumor.clinical.patient,clinical.var){
  PAAD.repertoire.tumor.clinical.patient.Ig<-PAAD.repertoire.tumor.clinical.patient[which(PAAD.repertoire.tumor.clinical.patient$Ig_Reads>1000),]
  Ig_expr<-melt(PAAD.repertoire.tumor.clinical.patient.Ig[,c("TCGA_sample","IGH_expression","IGK_expression","IGL_expression",clinical.var)])
  Ig_expr$value<-log10(Ig_expr$value)
  Ig_expr<-na.omit(Ig_expr)
  tiff(paste0("Results/ImmuneRep/Clinical//boxplot_Ig_expression_TCGA_",clinical.var,".tiff"),res=300,h=1000,w=2000)
  if(length(table(Ig_expr[,clinical.var]))<10){
    print(ggboxplot(Ig_expr, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            geom_point(aes(x=Ig_expr[,clinical.var], y=value,color=Ig_expr[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
            stat_compare_means(label = "p.format"))
  } else {
    print(ggplot(Ig_expr,aes(x = as.numeric(as.character(Ig_expr[,clinical.var])), y = value)) + geom_point() +
            facet_grid(.~Ig_expr$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()

  Ig_entropy<-melt(PAAD.repertoire.tumor.clinical.patient.Ig[,c("TCGA_sample","entropy_IGH","entropy_IGK","entropy_IGL",clinical.var)])
  Ig_entropy<-Ig_entropy[which(Ig_entropy$value!=0),]
  Ig_entropy<-na.omit(Ig_entropy)

  tiff(paste0("Results/ImmuneRep/Clinical//boxplot_Ig_entropy_TCGA_",clinical.var,".tiff"),res=300,h=1000,w=2000)
  if(length(table(Ig_entropy[,clinical.var]))<10){
    print(ggboxplot(Ig_entropy, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            geom_point(aes(x=Ig_entropy[,clinical.var], y=value,color=Ig_entropy[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
            stat_compare_means(label = "p.format"))
  } else {
    print(ggplot(Ig_entropy,aes(x = as.numeric(as.character(Ig_entropy[,clinical.var])), y = value)) + geom_point() +
            facet_grid(.~Ig_entropy$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()
  
  
  Ig_cluster<-melt(PAAD.repertoire.tumor.clinical.patient.Ig[,c("TCGA_sample","cluster_gini_IGH","cluster_gini_IGK","cluster_gini_IGL",clinical.var)])
  Ig_cluster<-na.omit(Ig_cluster)
  tiff(paste0("Results/ImmuneRep/Clinical//boxplot_Ig_cluster_TCGA_",clinical.var,".tiff"),res=300,h=1000,w=2000)
  if(length(table(Ig_cluster[,clinical.var]))<10){
    print(ggboxplot(Ig_cluster, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            geom_point(aes(x=Ig_cluster[,clinical.var], y=value,color=Ig_cluster[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
            stat_compare_means(label = "p.format"))
  } else {
    print(ggplot(Ig_cluster,aes(x = as.numeric(Ig_cluster[,clinical.var]), y = value)) + geom_point() +
            facet_grid(.~Ig_cluster$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()

  Ig_vertex<-melt(PAAD.repertoire.tumor.clinical.patient.Ig[,c("TCGA_sample","vertex_gini_IGH","vertex_gini_IGK","vertex_gini_IGL",clinical.var)])
  Ig_vertex<-na.omit(Ig_vertex)
  tiff(paste0("Results/ImmuneRep/Clinical//boxplot_Ig_vertex_TCGA_",clinical.var,".tiff"),res=300,h=1000,w=2000)
  if(length(table(Ig_vertex[,clinical.var]))<10){
    print(ggboxplot(Ig_vertex, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            geom_point(aes(x=Ig_vertex[,clinical.var], y=value,color=Ig_vertex[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
            stat_compare_means(label = "p.format"))
  } else {
    print(ggplot(Ig_vertex,aes(x = as.numeric(Ig_vertex[,clinical.var]), y = value)) + geom_point() +
            facet_grid(.~Ig_vertex$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()

  
  PAAD.repertoire.tumor.clinical.patient.T<-PAAD.repertoire.tumor.clinical.patient[which(PAAD.repertoire.tumor.clinical.patient$T_Reads>100),]
  T_expr<-melt(PAAD.repertoire.tumor.clinical.patient.T[,c("TCGA_sample","TRA_expression","TRB_expression","TRD_expression","TRG_expression",clinical.var)])
  T_expr$value<-log10(T_expr$value)
  T_expr<-na.omit(T_expr)
  tiff(paste0("Results/ImmuneRep/Clinical//boxplot_T_expr_TCGA_",clinical.var,".tiff"),res=300,h=1000,w=2500)
  if(length(table(T_expr[,clinical.var]))<10){
    print(ggboxplot(T_expr, x = clinical.var, y = "value",color = clinical.var,ggtheme = theme_bw()) + facet_wrap("variable",nrow=1)+
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            geom_point(aes(x=T_expr[,clinical.var], y=value,color=T_expr[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
            stat_compare_means(label = "p.format"))
  } else {
    print(ggplot(T_expr,aes(x = as.numeric(as.character(T_expr[,clinical.var])), y = value)) + geom_point() +
            facet_grid(.~T_expr$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  }
  dev.off()
  # 
  # T_entropy<-melt(PAAD.repertoire.tumor.clinical.patient.T[,c("TCGA_sample","entropy_TRA","entropy_TRB",clinical.var)])
  # T_entropy<-T_entropy[which(T_entropy$value!=0),]
  # T_entropy<-na.omit(T_entropy)
  # tiff(paste0("Results/ImmuneRep/Clinical//boxplot_T_entropy_TCGA_",clinical.var,".tiff"),res=300,h=1000,w=2000)
  # if(length(table(T_entropy[,clinical.var]))<10){  
  #   print(ggboxplot(T_entropy, x = clinical.var, y = "value",color = clinical.var,ggtheme = theme_bw()) +facet_wrap("variable",nrow=1)+
  #           theme(axis.title.x=element_blank(),
  #                 axis.text.x=element_blank(),
  #                 axis.ticks.x=element_blank()) +
  #           geom_point(aes(x=T_entropy[,clinical.var], y=value,color=T_entropy[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
  #           stat_compare_means(label = "p.format"))
  # } else {  
  #   print(ggplot(T_entropy,aes(x = as.numeric(as.character(T_entropy[,clinical.var])), y = value)) + geom_point() + 
  #           facet_grid(.~T_entropy$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  # }
  # dev.off()
  # 
  
  # Alpha_Beta_ratio_expression<-melt(PAAD.repertoire.tumor.clinical.patient.T[,c("TCGA_sample","Alpha_Beta_ratio_expression",clinical.var)])
  # Alpha_Beta_ratio_expression<-na.omit(Alpha_Beta_ratio_expression)
  # tiff(paste0("Results/ImmuneRep/Clinical//boxplot_Alpha_Beta_ratio_expression_TCGA_",clinical.var,".tiff"),res=300,h=1000,w=2000)
  # if(length(table(Alpha_Beta_ratio_expression[,clinical.var]))<10){  
  #   print(ggboxplot(Alpha_Beta_ratio_expression, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
  #           theme(axis.title.x=element_blank(),
  #                 axis.text.x=element_blank(),
  #                 axis.ticks.x=element_blank())+
  #           geom_point(aes(x=Alpha_Beta_ratio_expression[,clinical.var], y=value,color=Alpha_Beta_ratio_expression[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
  #           stat_compare_means(label = "p.format"))
  # } else {  
  #   print(ggplot(Alpha_Beta_ratio_expression,aes(x = as.numeric(Alpha_Beta_ratio_expression[,clinical.var]), y = value)) + geom_point() + 
  #           facet_grid(.~Alpha_Beta_ratio_expression$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson") + xlab(clinical.var))
  # }
  # dev.off()
  # 
  # KappaLambda_ratio_expression<-melt(PAAD.repertoire.tumor.clinical.patient.Ig[,c("TCGA_sample","KappaLambda_ratio_expression",clinical.var)])
  # KappaLambda_ratio_expression<-na.omit(KappaLambda_ratio_expression)
  # tiff(paste0("Results/ImmuneRep/Clinical//boxplot_KappaLambda_ratio_expression_TCGA_",clinical.var,".tiff"),res=300,h=1000,w=2000)
  # if(length(table(KappaLambda_ratio_expression[,clinical.var]))<10){  
  #   print(ggboxplot(KappaLambda_ratio_expression, x = clinical.var, y = "value",facet.by = "variable",color = clinical.var,ggtheme = theme_bw()) +
  #           theme(axis.title.x=element_blank(),
  #                 axis.text.x=element_blank(),
  #                 axis.ticks.x=element_blank()) +
  #           geom_point(aes(x=KappaLambda_ratio_expression[,clinical.var], y=value,color=KappaLambda_ratio_expression[,clinical.var]), position = position_jitterdodge(dodge.width = 0.8)) +
  #           stat_compare_means(label = "p.format"))
  # } else {  
  #   print(ggplot(KappaLambda_ratio_expression,aes(x = as.numeric(KappaLambda_ratio_expression[,clinical.var]), y = value)) + geom_point() + 
  #           facet_grid(.~KappaLambda_ratio_expression$variable) + geom_smooth(method="lm") + stat_cor(method = "pearson")+ xlab(clinical.var))
  # }
  # dev.off()
  # 
  # 
}

############################
### Plots for the paper ###
###########################
PAAD.repertoire.tumor.filter.Ig<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$Ig_Reads>1000),]
df<-PAAD.repertoire.tumor.filter.Ig[,c("entropy_IGH","number_pack_years_smoked")]
tiff("Manuscript/Figures/entropy_IGH_pack_years.tiff",res=300,h=1000,w=1300)
ggplot(df,aes(x = as.numeric(as.character(number_pack_years_smoked)), y = entropy_IGH)) + geom_point() +
 geom_smooth(method="lm") + stat_cor(method = "pearson",size=5)+ xlab("number_pack_years_smoked") +
   theme(axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
         axis.title.x = element_text(size=20),
         axis.title.y = element_text(size=20))+
  xlab("# pack years smoked") + xlim (0,60) +
  ylab("entropy IGH") + ylim (3,11) 
dev.off()  

PAAD.repertoire.tumor.filter.Ig<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$Ig_Reads>1000),]
df<-PAAD.repertoire.tumor.filter.Ig[,c("IGK_expression","gender")]
df$IGK_expression<-log10(df$IGK_expression)
tiff("Manuscript/Figures/IGK_expression_sex.tiff",res=300,h=1000,w=1000)
ggboxplot(df, x = "gender", y = "IGK_expression",color = "gender",ggtheme = theme_bw()) +
  geom_point(aes(x=df$gender, y=df$IGK_expression,color=df$gender), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(label = "p.format",size=8,label.y = -2.3) + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20)) +
  xlab("Sex") +
  ylab("IGK expression") + ylim (-6,-2) +
  theme(legend.position = "none") 

dev.off()  

PAAD.repertoire.tumor.filter.T<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$T_Reads>100),]
df<-PAAD.repertoire.tumor.filter.T[,c("TRB_expression","KRAS.Mutated..1.or.0.")]
df$TRB_expression<-log10(df$TRB_expression)
df$KRAS.Mutated..1.or.0.<-factor(df$KRAS.Mutated..1.or.0.)
tiff("Manuscript/Figures/TRB_expression_KRAS.tiff",res=300,h=1100,w=1100)
ggboxplot(df, x = "KRAS.Mutated..1.or.0.", y = "TRB_expression",color = "KRAS.Mutated..1.or.0.",ggtheme = theme_bw()) +
  geom_point(aes(x=df$KRAS.Mutated..1.or.0., y=df$TRB_expression,color=df$KRAS.Mutated..1.or.0.), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(label = "p.format",size=8,label.y=-4.5) + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20)) +
  xlab("KRAS mutation") +
  ylab("TRB expression") + ylim (-6,-4) +
  theme(legend.position = "none") 
dev.off()  

PAAD.repertoire.tumor.filter$history_of_diabetes<-ifelse(PAAD.repertoire.tumor.filter$history_of_diabetes=="",NA,
                                                         as.character(PAAD.repertoire.tumor.filter$history_of_diabetes))
PAAD.repertoire.tumor.filter.T<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$T_Reads>100),]
df<-PAAD.repertoire.tumor.filter.T[,c("TRB_expression","history_of_diabetes")]
df$TRB_expression<-log10(df$TRB_expression)
df<-na.omit(df)
tiff("Manuscript/Figures/TRB_expression_diabetes.tiff",res=300,h=1100,w=1100)
ggboxplot(df, x = "history_of_diabetes", y = "TRB_expression",color = "history_of_diabetes",ggtheme = theme_bw()) +
  geom_point(aes(x=df$history_of_diabetes, y=df$TRB_expression,color=df$history_of_diabetes), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(label = "p.format",size=8,label.y=-4.5) + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20)) +
  xlab("history of diabetes") +
  ylab("TRB expression") + ylim (-6,-4) +
  theme(legend.position = "none") 
dev.off()  

PAAD.repertoire.tumor.filter$history_of_diabetes<-ifelse(PAAD.repertoire.tumor.filter$history_of_diabetes=="",NA,
                                                         as.character(PAAD.repertoire.tumor.filter$history_of_diabetes))
association.test.immuneRep(PAAD.repertoire.tumor.filter,"history_of_diabetes")


PAAD.repertoire.tumor.filter.T<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$T_Reads>100),]
df<-PAAD.repertoire.tumor.filter.T[,c("TRB_expression","SNV.Neoantigens")]
df$TRB_expression<-log10(df$TRB_expression)
tiff("Manuscript/Figures/TRB_expression_SNV.neoantigens.tiff",res=300,h=1000,w=1400)
ggplot(df,aes(x = as.numeric(as.character(SNV.Neoantigens)), y = TRB_expression)) + geom_point() +
  geom_smooth(method="lm") + stat_cor(method = "pearson",size=5)+ xlab("SNV.Neoantigens") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20)) +
  xlab("SNV.Neoantigens") + 
  ylab("TRB expression") + ylim (-6,-3.5)
dev.off()


PAAD.repertoire.tumor.filter.T<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$T_Reads>100),]
df<-PAAD.repertoire.tumor.filter.T[,c("TRA_expression","ABSOLUTE.Purity")]
df$TRA_expression<-log10(df$TRA_expression)
col<-brewer.pal(9, "Set1")
tiff("Manuscript/Figures/TRA_expression_ABSOLUTE.tiff",res=300,h=1100,w=1500)
ggplot(df,aes(x = as.numeric(as.character(ABSOLUTE.Purity)), y = TRA_expression)) + geom_point(color = col[1]) +
  geom_smooth(method="lm",color = col[1]) + stat_cor(method = "pearson",size=5,color="black",label.y.npc="top", label.x.npc = 0.4)+ xlab("ABSOLUTE.Purity") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))  + ylim(-7,-4) +
  xlab("ABSOLUTE.Purity") + 
  ylab("TRA Expression (log10)") 
dev.off()

PAAD.repertoire.tumor.filter.Ig<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$Ig_Reads>100),]
df<-PAAD.repertoire.tumor.filter.Ig[,c("IGK_clonotypes_cluster","TGF.beta.Response")]
tiff("Manuscript/Figures/IGK_clonotypes_cluster_ABSOLUTE.Purity.tiff",res=300,h=1800,w=1300)
ggboxplot(df, x = "IGK_clonotypes_cluster", y = "TGF.beta.Response",color = "IGK_clonotypes_cluster",ggtheme = theme_bw()) +
  geom_point(aes(x=df$IGK_clonotypes_cluster, y=df$TGF.beta.Response,color=df$IGK_clonotypes_cluster), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(label = "p.format",size=8) + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20)) +
  xlab("IGK_clonotypes_cluster") +
  ylab("TGF.beta.Response") +
  theme(legend.position = "none") 
dev.off()


PAAD.repertoire.tumor.filter.Ig<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$T_Reads>100),]
df<-PAAD.repertoire.tumor.filter.Ig[,c("TRB_expression","Macrophage.Regulation")]
df$TRB_expression<-log10(df$TRB_expression)
col<-brewer.pal(9, "Set1")
tiff("Manuscript/Figures/TRB_expression_Macrophage.Regulation.tiff",res=300,h=1100,w=1400)
ggplot(df,aes(x = as.numeric(as.character(Macrophage.Regulation)), y = TRB_expression)) + geom_point(color = col[4]) +
  geom_smooth(method="lm",color = col[4]) + stat_cor(method = "pearson",size=5,color="black",label.y.npc= "top")+ xlab("Leukocyte_DNAmethylation") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))  + ylim(-7,-4) + xlim(-0.5,2) +
  xlab("Macrophage Regulation") + 
  ylab(NULL) + theme(axis.text.y = element_blank())
dev.off()


PAAD.repertoire.tumor.filter.T<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$T_Reads>100),]
df<-PAAD.repertoire.tumor.filter.T[,c("entropy_TRA","subtypes_Bailey")]
df$TRA_expression<-log10(df$TRA_expression)
df<-df[which(df$entropy_TRA!=0),]
tiff("Manuscript/Figures/entropy_TRA_subtypes_Bailey.tiff",res=300,h=1100,w=1200)
ggboxplot(df, x = "subtypes_Bailey", y = "entropy_TRA",color = "subtypes_Bailey",ggtheme = theme_bw()) +
  geom_point(aes(x=subtypes_Bailey, y=entropy_TRA,color=subtypes_Bailey), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(label = "p.format",size=5) + 
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(face="bold",size=15)) + 
  xlab("subtypes_Bailey") + 
  #theme(legend.position = "top",legend.direction = "horizontal")+
  ylab(NULL) + theme(axis.text.y = element_blank()) + theme(axis.text.x = element_blank()) + xlab(NULL)
dev.off()  


PAAD.repertoire.tumor.filter.Ig<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$Ig_Reads>1000),]
df<-PAAD.repertoire.tumor.filter.Ig[,c("IGH_expression","Th17.Cells")]
df$IGH_expression<-log10(df$IGH_expression)
col<-brewer.pal(9, "Set1")
tiff("Manuscript/Figures/IGH_expression_Th17.Cells.tiff",res=300,h=1100,w=1500)
ggplot(df,aes(x = as.numeric(as.character(Th17.Cells)), y = IGH_expression)) + geom_point(color = col[3]) +
  geom_smooth(method="lm",color = col[3]) + stat_cor(method = "pearson",size=5,color="black",label.y.npc= "top") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))  + ylim(-6,-2) +
  xlab("Th17 Cells") +
  ylab("IGH Expression (log10)")
  #ylab(NULL) + theme(axis.text.y = element_blank())

dev.off()  

PAAD.repertoire.tumor.filter.T<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$T_Reads>100),]
df<-PAAD.repertoire.tumor.filter.T[,c("TRB_expression","IFN.gamma.Response")]
df$TRB_expression<-log10(df$TRB_expression)
col<-brewer.pal(9, "Set1")
tiff("Manuscript/Figures/TRB_expression_IFN.gamma.Response.tiff",res=300,h=1100,w=1500)
ggplot(df,aes(x = as.numeric(as.character(IFN.gamma.Response)), y = TRB_expression)) + geom_point(color = col[4]) +
  geom_smooth(method="lm",color = col[4]) + stat_cor(method = "pearson",size=5,color="black",label.y.npc= "top") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))  + ylim(-7,-4) +
  xlab("IFN gamma Response") +
  ylab("TRB Expression (log10)")
#ylab(NULL) + theme(axis.text.y = element_blank())

dev.off()  



PAAD.repertoire.tumor.filter.T<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$T_Reads>100),]
df<-PAAD.repertoire.tumor.filter.T[,c("entropy_TRB","Nonsilent.Mutation.Rate")]
df<-df[which(df$entropy_TRB!=0),]
#df$TRB_expression<-log10(df$TRB_expression)
col<-brewer.pal(8, "Set3")
tiff("Manuscript/Figures/entropy_TRB_Nonsilent.Mutation.Rates.tiff",res=300,h=1100,w=1500)
ggplot(df,aes(x = as.numeric(as.character(Nonsilent.Mutation.Rate)), y = entropy_TRB)) + geom_point(color = col[4]) +
  geom_smooth(method="lm",color = col[4]) + stat_cor(method = "pearson",size=5,color="black",label.y.npc= "top") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))  + 
  xlab("Non-Silent Mutation Rate") +
  ylab("Shannon Entropy TRB")
#ylab(NULL) + theme(axis.text.y = element_blank())

dev.off()  

PAAD.repertoire.tumor.filter$history_of_diabetes<-factor(PAAD.repertoire.tumor.filter$history_of_diabetes)
PAAD.repertoire.tumor.filter.T<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$T_Reads>100),]
df<-PAAD.repertoire.tumor.filter.T[,c("TRB_expression","history_of_diabetes")]
df<-df[which(is.na(df$history_of_diabetes)==F),]
df$TRB_expression<-log10(df$TRB_expression)
#df<-df[which(df$entropy_TRB!=0),]
tiff("Manuscript/Figures/TRB_expression_history_of_diabetes.tiff",res=300,h=1000,w=1100)
ggboxplot(df, x = "history_of_diabetes", y = "TRB_expression",color = "history_of_diabetes",ggtheme = theme_bw()) +
  geom_point(aes(x=history_of_diabetes, y=TRB_expression,color=history_of_diabetes), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(label = "p.format",size=5) + scale_colour_discrete("Diabetes") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(face="bold",size=15)) + ylim(-7,-4)+
  #theme(legend.position = "top",legend.direction = "horizontal")+
  ylab(NULL) + theme(axis.text.y = element_blank()) + theme(axis.text.x = element_blank()) + xlab(NULL)
dev.off()  



PAAD.repertoire.tumor.filter.Ig<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$Ig_Reads>1000),]
df<-PAAD.repertoire.tumor.filter.Ig[,c("vertex_gini_IGK","Number.of.pack.years.smoked")]
#df<-df[which(df$entropy_IGK!=0),]
#df$IGL_expression<-log10(df$IGL_expression)
col<-brewer.pal(8, "Dark2")
tiff("Manuscript/Figures/vertex_gini_IGK_number_pack_years_smoked.tiff",res=300,h=1100,w=1500)
ggplot(df,aes(x = as.numeric(as.character(Number.of.pack.years.smoked)), y = vertex_gini_IGK)) + geom_point(color = col[1]) +
  geom_smooth(method="lm",color = col[1]) + stat_cor(method = "pearson",size=5,color="black",label.y.npc= "top") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))  + #ylim(-6,-2)+
  xlab("# pack years smoked") +
  ylab("Gini (vertex) IGK")
#ylab(NULL) + theme(axis.text.y = element_blank())

dev.off()  


PAAD.repertoire.tumor.filter.Ig<-PAAD.repertoire.tumor.filter[which(PAAD.repertoire.tumor.filter$Ig_Reads>1000),]
df<-PAAD.repertoire.tumor.filter.Ig[,c("cluster_gini_IGK","Gender")]
#df$IGH_expression<-log10(df$IGH_expression)
#df<-df[which(df$entropy_IGL!=0),]
tiff("Manuscript/Figures/cluster_gini_IGK_gender.tiff",res=300,h=1000,w=1100)
ggboxplot(df, x = "Gender", y = "cluster_gini_IGK",color = "Gender",ggtheme = theme_bw()) +
  geom_point(aes(x=Gender, y=cluster_gini_IGK,color=Gender), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(label = "p.format",size=5) + scale_colour_discrete("Sex") +
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(face="bold",size=15)) + #ylim(-6,-2)+
  #theme(legend.position = "top",legend.direction = "horizontal")+
  ylab(NULL) + theme(axis.text.y = element_blank()) + theme(axis.text.x = element_blank()) + xlab(NULL)
dev.off() 
