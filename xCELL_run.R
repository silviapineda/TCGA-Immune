
PAAD_fpkm<-read.table("~/TCGA-Immune/Data/xCELL/TCGA-PAAD.htseq_fpkm.tsv",header = T)
genecode<-read.table("~/TCGA-Immune/Data/xCELL/gencode.v22.annotation.gene.probeMap",header = T)

id.gene<-match(PAAD_fpkm$Ensembl_ID,genecode$id)
PAAD_fpkm$Gene<-genecode$gene[id.gene]
write.csv(PAAD_fpkm,"~/TCGA-Immune/Data/xCELL/PAAD_fpkm.csv")
