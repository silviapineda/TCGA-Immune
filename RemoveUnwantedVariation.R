rm(list = ls(all = TRUE))
x<-date()
print(x)
###########################################################################################
### PROJECT: TCGA Immune project. Compositional Analysis
###
### CITATION: 
###
### PROCESS: 
###           
### DESCRIP: http://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf
### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4404308/
###         
###
### Author: Silvia Pineda
### Date: July, 2019
############################################################################################

case_control_reads<-read.csv("~/TCGA-Immune/Data/case_control_reads.csv")
case_control_reads$Ca_Co<-factor(case_control_reads$Ca_Co)

library(RUVSeq)
matrix<-t(as.matrix(case_control_reads[,c("IGH","IGK","IGL","TRA","TRB","TRD","TRG")]))
colnames(matrix)<-(case_control_reads$sample_SRR)
#filter <- apply(matrix, 1, function(x) length(x[x>5])>=2)
#filtered <- matrix[filter,]
set <- newSeqExpressionSet(matrix,phenoData = data.frame(case_control_reads$Ca_Co, row.names=case_control_reads$sample_SRR))

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, col=colors[case_control_reads$Ca_Co])
plotPCA(set, col=colors[case_control_reads$Ca_Co], cex=1.2)

###upper-quartile normalization
set2 <- betweenLaneNormalization(set, which="upper")
plotRLE(set2, outline=FALSE, col=colors[case_control_reads$Ca_Co])
plotPCA(set2, col=colors[case_control_reads$Ca_Co], cex=1.2)

###RUVr
design <- model.matrix(~case_control_reads$Ca_Co, data=pData(set))
y <- DGEList(counts=counts(set), group=case_control_reads$Ca_Co)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

genes<-c("IGH","IGK","IGL","TRA","TRB","TRD","TRG")
set4 <- RUVr(set, genes, k=1, res)
data<-pData(set4)
plotRLE(set4, outline=FALSE, col=colors[case_control_reads$Ca_Co])
plotPCA(set4, col=colors[case_control_reads$Ca_Co], cex=1.2)

#Edge2 package
design <- model.matrix(~case_control_reads$Ca_Co + W_1, data=pData(set4))
y <- DGEList(counts=counts(set4), group=case_control_reads$Ca_Co)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)

#DSEq2
dds <- DESeqDataSetFromMatrix(countData = counts(set4),
                              colData = pData(set4),
                              design = ~ W_1 + case_control_reads.Ca_Co)
dds <- DESeq(dds)
res <- results(dds)
res
dds <- DESeq(dds, test="LRT", reduced=as.formula("~ W_1"))
res <- results(dds)
res

rlog<-rlog(dds)
data_normalized<-assay(rlog)

cols=brewer.pal(3,name = "Accent")
df<-data.frame(t(data_normalized))
library(ggplot2)
library(ggpubr)
ggplot(df) + 
  geom_boxplot(aes(y=IGH, color=case_control_reads$Ca_Co, x=case_control_reads$Ca_Co), alpha = 0, position = position_dodge(width = .8)) +
  geom_point(aes(y=IGH, color=case_control_reads$Ca_Co,x=case_control_reads$Ca_Co), position = position_jitterdodge(dodge.width = 0.8)) +     
  stat_compare_means(aes(x=case_control_reads$Ca_Co, y=IGH, color=case_control_reads$Ca_Co),label.x.npc = "center",method = "anova") +
  scale_color_manual(name = "otucome",values = c(cols[1], cols[2]), labels = c("Normal_pancreas", "Tumor_pancreas")) +
  scale_x_discrete(name = "outcome", labels=c("Normal_pancreas", "Tumor_pancreas"))


#####################################
###### RemoveUnwantedVariation #####
###################################
setwd("~/TCGA-Immune/")

load("Data/PAAD/PAAD_FullData.Rdata")

##Filter for those that are pancreas
PAAD.repertoire.diversity.tumor<-PAAD.repertoire.diversity[which(PAAD.repertoire.diversity$Tumor_type_2categ=="Tumor_pancreas"),]
data_merge_qc<-data_merge[which(is.na(match(data_merge$sample,rownames(PAAD.repertoire.diversity.tumor)))==F),]
#Prepare vgenes matrix for the TCGA
vgenes<-as.data.frame(unclass(table(data_merge_qc$sample,data_merge_qc$bestVGene)))
vgenes<-vgenes[,-1]
vgenes_filter<-vgenes[apply(vgenes, 1, function(x) length(x[x>5])>2),]

#Prepare vgenes matrix for GTEX
data_merge_gtex<-read.csv("Data/PAAD/GTEX_data_pancreas.csv")
vgenes_gtex<-as.data.frame(unclass(table(data_merge_gtex$sample_SRR,data_merge_gtex$bestVGene)))
vgenes_gtex<-vgenes_gtex[,-1]
vgenes_gtex_filter<-vgenes_gtex[apply(vgenes_gtex, 1, function(x) length(x[x>5])>2),]

##Merge both datasets
id.common<-match(colnames(vgenes),colnames(vgenes_gtex))
vgenes_total<-rbind(vgenes[,which(is.na(id.common)==F)],vgenes_gtex[,na.omit(id.common)])
caco<-c(rep("TCGA",dim(vgenes)[1]),rep("GTEX",dim(vgenes_gtex)[1]))
names(caco)<-rownames(vgenes_total)

vgenes_total_matrix<-t(as.matrix(vgenes_total))
filter <- apply(vgenes_total_matrix, 1, function(x) length(x[x>5])>10)
vgenes_filtered <- vgenes_total_matrix[filter,]


set <- newSeqExpressionSet(vgenes_total_matrix,phenoData = data.frame(caco, row.names=names(caco)))

colors <- brewer.pal(3,name = "Accent")
plotRLE(set, outline=FALSE, col=colors[factor(caco)])
plotPCA(set, col=colors[factor(caco)], cex=1.2)

###upper-quartile normalization
set2 <- betweenLaneNormalization(set, which="upper")
plotRLE(set2, outline=FALSE, col=colors[factor(caco)])
plotPCA(set2, col=colors[factor(caco)], cex=1.2)

###RUVr
design <- model.matrix(~caco, data=pData(set))
y <- DGEList(counts=counts(set), group=caco)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

genes<-colnames(vgenes_total)
set4 <- RUVr(set, genes, k=1, res)
data<-pData(set4)
plotRLE(set4, outline=FALSE, col=colors[case_control_reads$Ca_Co])
plotPCA(set4, col=colors[case_control_reads$Ca_Co], cex=1.2)

#Edge2 package
design <- model.matrix(~case_control_reads$Ca_Co + W_1, data=pData(set4))
y <- DGEList(counts=counts(set4), group=case_control_reads$Ca_Co)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
