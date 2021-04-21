
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
# copy the repository from https://github.com/UVic-omics/CoDA-Penalized-Regression
system('git clone https://github.com/UVic-omics/CoDA-Penalized-Regression')
# cran.packages <- c('knitr', 'glmnet', 'ggplot2', 'gridExtra',
#                    'UpSetR', 'ggforce')
# install.packages(cran.packages)
devtools::install_github(repo = 'UVic-omics/selbal')

library(knitr) # rbookdown, kable
library(glmnet) # glmnet
library(selbal) # selbal
library(ggplot2) # draw selbal
library(gridExtra) # grid.arrange
library(UpSetR) # upset
library(ggforce) # selbal-like plot
library(grid) # grid.draw
library(pheatmap)
library(RColorBrewer)
# source coda-lasso functions
source(file = './CoDA-Penalized-Regression/R/functions_coda_penalized_regression.R')


setwd("~/TCGA-Immune/")
Kraken_counts<-read.csv("Microbiome/kraken_counts_sorted_allreadsBacteria_filter.csv")
rownames(Kraken_counts)<-Kraken_counts$species
Kraken_counts<-Kraken_counts[,-1]

###Annotate the microbiome data with the TCGA samples
annotation_microbiome<-read.csv("Microbiome/gdc_sample_sheet.2020-05-22.csv")
id<-match(substr(rownames(Kraken_counts),1,36),substr(annotation_microbiome$File.Name,1,36))
rownames(Kraken_counts)<-annotation_microbiome$Sample.ID[id]
rownames(Kraken_counts)<-substr(rownames(Kraken_counts),1,15)

load("Data/PAAD/PAAD_FullData.Rdata")
id<-match(rownames(Kraken_counts),substr(PAAD.repertoire.diversity$TCGA_sample,1,15))
Kraken_counts_filter<-Kraken_counts[which(is.na(id)==F),]
PAAD.repertoire.microbiome<-PAAD.repertoire.diversity[na.omit(id),]
rownames(PAAD.repertoire.microbiome)<-substr(PAAD.repertoire.microbiome$TCGA_sample,1,15)

##Build present vs no present
Kraken_counts_filter_present<-apply(Kraken_counts_filter,1,function(x) ifelse(x==0,0,1))
###Filter by clones share at least in 2% of the samples (311*0.05 = 15) or 2 samples
Kraken_counts_filter_2<-Kraken_counts_filter[,which(rowSums(Kraken_counts_filter_present)>dim(Kraken_counts_filter)[1]*0.1)] #
#####################
### Add an offset of 1 to substitute the zeros
Kraken_counts_filter_zerosubs<-Kraken_counts_filter_2+1
################
### CLR-lasso ##
################
##
z_Kraken_counts<-log(Kraken_counts_filter_zerosubs)
clrx_Kraken_counts <- apply(z_Kraken_counts, 2, function(x) x - rowMeans(z_Kraken_counts))

id<-match(rownames(clrx_Kraken_counts),rownames(PAAD.repertoire.microbiome))
set.seed(35)
clrlasso.cv <- cv.glmnet(x = clrx_Kraken_counts, y = PAAD.repertoire.microbiome$Tumor_type_2categ[id], 
                                  family = 'binomial', nfolds = 5, alpha=1)

test_clrlasso <- glmnet(x = clrx_Kraken_counts, y = PAAD.repertoire.microbiome$Tumor_type_2categ[id], 
                            family = 'binomial', alpha=1,lambda=clrlasso.cv$lambda.min)
##1 selected
bacteria_clr<-rownames(test_clrlasso$beta)[which(as.numeric(test_clrlasso$beta)!=0)]
clrx_Kraken_counts[,match(bacteria_clr,colnames(clrx_Kraken_counts))]
clrx_Kraken_counts_sign<-clrx_Kraken_counts[,match(bacteria_clr,colnames(clrx_Kraken_counts))]
boxplot(clrx_Kraken_counts_sign~PAAD.repertoire.microbiome$Tumor_type_2categ[id])

##Heatmap 
brewer.pal(4,name = "Accent")
cols=brewer.pal(4,name = "Accent")

annotation_row = data.frame(PAAD.repertoire.microbiome$Tumor_type_3categ[id])
#ann_colors = list (Sample.Type = c("normal_pancreas (GTEx)" = cols[1],"PDAC (TCGA)" = cols[2]))
colnames(annotation_row)<-"Sample.Type"
rownames(annotation_row)<-rownames(clrx_Kraken_counts)

tiff("Microbiome/Kraken_CLR.tiff",width = 3200, height = 2500, res = 300)
pheatmap(t(clrx_Kraken_counts),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
         color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))
dev.off()

#### Clustering ###
res <- pheatmap(t(clrx_Kraken_counts),scale="row",border_color=F,show_colnames = F, annotation_col = annotation_row,
                color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120))

mat.clust <- as.data.frame(cbind(clrx_Kraken_counts, cluster = cutree(res$tree_col, k = 6)))

annotation_row = data.frame(cluster = factor(mat.clust$cluster))
colnames(annotation_row)<-c("cluster")
rownames(annotation_row)<-rownames(clrx_Kraken_counts)
tiff("Microbiome/Kraken_CLR_clustering.tiff",width = 3200, height = 3000, res = 300)
pheatmap(t(clrx_Kraken_counts),scale="row",border_color=F, show_colnames = F,annotation_col = annotation_row,
         color = colorRampPalette(rev(brewer.pal(6,name="RdGy")))(120),cutree_cols = 6,cutree_rows = 6)
dev.off()

PAAD.repertoire.microbiome$Bacteria_cluster<-mat.clust$cluster

#########################################
## Integration with Immune Repertooire ##
#########################################


#####Filter by only PDAC
id<-match(rownames(Kraken_counts),substr(PAAD.repertoire.tumor.filter$TCGA_sample,1,15))
Kraken_counts_filter<-Kraken_counts[which(is.na(id)==F),]
PAAD.repertoire.microbiome.tumor<-PAAD.repertoire.tumor.filter[na.omit(id),]
rownames(PAAD.repertoire.microbiome.tumor)<-substr(PAAD.repertoire.microbiome.tumor$TCGA_sample,1,15)

##Build present vs no present
Kraken_counts_filter_present<-apply(Kraken_counts_filter,1,function(x) ifelse(x==0,0,1))
###Filter by clones share at least in 2% of the samples (311*0.05 = 15) or 2 samples
Kraken_counts_filter_2<-Kraken_counts_filter[,which(rowSums(Kraken_counts_filter_present)>dim(Kraken_counts_filter)[1]*0.1)] #

####CLR transformation only for PDAC samples
Kraken_counts_filter_zerosubs<-Kraken_counts_filter_2+1
################
### CLR-Transformation ##
################
##
z_Kraken_counts<-log(Kraken_counts_filter_zerosubs)
clrx_Kraken_counts <- apply(z_Kraken_counts, 2, function(x) x - rowMeans(z_Kraken_counts))


p.value=NULL
for(i in 1:ncol(clrx_Kraken_counts)){
  p.value[i]<-kruskal.test(clrx_Kraken_counts[,i]~PAAD.repertoire.microbiome.tumor$IGK_clonotypes_cluster)$p.value
}
clrx_Kraken_counts_sign<-clrx_Kraken_counts[,p.value<0.05]
df<-data.frame(cbind(clrx_Kraken_counts_sign,factor(PAAD.repertoire.microbiome.tumor$IGK_clonotypes_cluster)))
colnames(df)[7]<-"IGK_clonotypes_cluster"
df.m <- melt(df,id.var="IGK_clonotypes_cluster")
df.m$IGK_clonotypes_cluster<-factor(df.m$IGK_clonotypes_cluster)
tiff("Results/Microbiome/clr_IGKcluster.tiff",res=200,h=1000,w=2000)
ggplot(data = df.m, aes(x=variable, y=value)) + geom_boxplot(aes(fill=IGK_clonotypes_cluster))+ coord_flip() + ylab("clr")
dev.off()













PAAD.repertoire.microbiome.tumor$Bacteria_cluster<-mat.clust$cluster

##Subtypes
PAAD.repertoire.microbiome.tumor$Moffitt.clusters..All.150.Samples..1basal..2classical<-factor(PAAD.repertoire.microbiome.tumor$mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical)
PAAD.repertoire.microbiome.tumor$subtypes_Moffit<-factor(ifelse(PAAD.repertoire.microbiome.tumor$Moffitt.clusters..All.150.Samples..1basal..2classical==1,"Basal",
                                                            ifelse(PAAD.repertoire.microbiome.tumor$Moffitt.clusters..All.150.Samples..1basal..2classical==2,"Classical",NA)))

PAAD.repertoire.microbiome.tumor$mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM<-factor(PAAD.repertoire.microbiome.tumor$mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM)
PAAD.repertoire.microbiome.tumor$subtypes_Collisson<-factor(ifelse(PAAD.repertoire.microbiome.tumor$mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM==1,"Classical",
                                                               ifelse(PAAD.repertoire.microbiome.tumor$mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM==2,"Exocrine",
                                                                      ifelse(PAAD.repertoire.microbiome.tumor$mRNA.Collisson.clusters..All.150.Samples..1classical.2exocrine.3QM==3,"QM",NA))))

PAAD.repertoire.microbiome.tumor$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX<-factor(PAAD.repertoire.microbiome.tumor$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX)
PAAD.repertoire.microbiome.tumor$subtypes_Bailey<-factor(ifelse(PAAD.repertoire.microbiome.tumor$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX==1,"Squamous",
                                                            ifelse(PAAD.repertoire.microbiome.tumor$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX==2,"Immunogenic",
                                                                   ifelse(PAAD.repertoire.microbiome.tumor$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX==3,"Progenitor",
                                                                          ifelse(PAAD.repertoire.microbiome.tumor$mRNA.Bailey.Clusters..All.150.Samples..1squamous.2immunogenic.3progenitor.4ADEX==4,"Aberrantly differentiated exocrine",NA)))))

######
##Analysis of Bacteria cluster with IR
######
PAAD.repertoire.microbiome.tumor$Bacteria_cluster<-factor(PAAD.repertoire.microbiome.tumor$Bacteria_cluster)
TR_expr<-melt(PAAD.repertoire.microbiome.tumor[,c("TRA_expression","TRB_expression","Bacteria_cluster")])
TR_expr$value<-log10(TR_expr$value)
tiff("Microbiome/Texpression_bacteriaCluster.tiff",res=300,h=2000,w=2000)
ggboxplot(TR_expr, x = "Bacteria_cluster", y = "value",facet.by = "variable",color = "Bacteria_cluster",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=Bacteria_cluster, y=value, color=Bacteria_cluster), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(
    comparisons =list(c("1","2"),c("1","3"),
                      c("1","4"),c("2","3"),
                      c("2","4"),c("3","4")))

dev.off()

Ig_expr<-melt(PAAD.repertoire.microbiome.tumor[,c("IGH_expression","IGK_expression","IGL_expression","Bacteria_cluster")])
Ig_expr$value<-log10(Ig_expr$value)
tiff("Microbiome/Igexpression_bacteriaCluster.tiff",res=300,h=2000,w=2000)
ggboxplot(Ig_expr, x = "Bacteria_cluster", y = "value",facet.by = "variable",color = "Bacteria_cluster",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=Bacteria_cluster, y=value, color=Bacteria_cluster), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(
    comparisons =list(c("1","2"),c("1","3"),
                      c("1","4"),c("2","3"),
                      c("2","4"),c("3","4")))

dev.off()

TR_entropy<-melt(PAAD.repertoire.microbiome.tumor[,c("entropy_TRA","entropy_TRB","Bacteria_cluster")])
TR_entropy<-TR_entropy[which(TR_entropy$value!=0),]
tiff("Microbiome/Tentropy_bacteriaCluster.tiff",res=300,h=2000,w=2000)
ggboxplot(TR_entropy, x = "Bacteria_cluster", y = "value",facet.by = "variable",color = "Bacteria_cluster",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=Bacteria_cluster, y=value, color=Bacteria_cluster), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(
    comparisons =list(c("1","2"),c("1","3"),
                      c("1","4"),c("2","3"),
                      c("2","4"),c("3","4")))

dev.off()

Ig_entropy<-melt(PAAD.repertoire.microbiome.tumor[,c("entropy_IGH","entropy_IGK","entropy_IGL","Bacteria_cluster")])
Ig_entropy<-Ig_entropy[which(Ig_entropy$value!=0),]
tiff("Microbiome/Igentropy_bacteriaCluster.tiff",res=300,h=2000,w=2000)
ggboxplot(Ig_entropy, x = "Bacteria_cluster", y = "value",facet.by = "variable",color = "Bacteria_cluster",ggtheme = theme_bw(),xlab = F) +
  rotate_x_text() +
  geom_point(aes(x=Bacteria_cluster, y=value, color=Bacteria_cluster), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(
    comparisons =list(c("1","2"),c("1","3"),
                      c("1","4"),c("2","3"),
                      c("2","4"),c("3","4")))

dev.off()

Ig_gini<-melt(PAAD.repertoire.microbiome.tumor[,c("cluster_gini_IGH","vertex_gini_IGH","cluster_gini_IGK","vertex_gini_IGK","cluster_gini_IGL","vertex_gini_IGL","Bacteria_cluster")])
tiff("Microbiome/Ig_gini_bacteriaCluster.tiff",res=300,h=3000,w=2000)
ggboxplot(Ig_gini, x = "Bacteria_cluster", y = "value",color = "Bacteria_cluster",ggtheme = theme_bw(),xlab = F) +
  facet_wrap(~variable,nrow=3) +
  rotate_x_text() +
  geom_point(aes(x=Bacteria_cluster, y=value, color=Bacteria_cluster), position = position_jitterdodge(dodge.width = 0.8)) +
  stat_compare_means(
    comparisons =list(c("1","2"),c("1","3"),
                      c("1","4"),c("2","3"),
                      c("2","4"),c("3","4")))

dev.off()


#########################
### Survival Analysis ###
######################### 
library(survival)
library(survminer)
library(survMisc)
##OS
PAAD.repertoire.microbiome.tumor_fileter<-PAAD.repertoire.microbiome.tumor[which(PAAD.repertoire.microbiome.tumor$Bacteria_cluster=="1" |
                                                                                   PAAD.repertoire.microbiome.tumor$Bacteria_cluster=="2" | 
                                                                                   PAAD.repertoire.microbiome.tumor$Bacteria_cluster=="3"),]
PAAD.repertoire.microbiome.tumor_fileter$Bacteria_cluster<-factor(PAAD.repertoire.microbiome.tumor_fileter$Bacteria_cluster)
surv_object <- Surv(time = PAAD.repertoire.microbiome.tumor_fileter$OS.time, event = PAAD.repertoire.microbiome.tumor_fileter$OS)
res.cox <- coxph(Surv(time = OS.time, event = OS)~Bacteria_cluster,data=PAAD.repertoire.microbiome.tumor_fileter)
summary(res.cox)
ggforest(res.cox)


