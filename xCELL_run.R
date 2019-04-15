
###Download the RSEM FPKM values from the xenabrowser
#https://xenabrowser.net/datapages/?dataset=tcga_RSEM_gene_fpkm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
xcell<-read.csv("Data/xCELL/tcga_RSEM_gene_fpkm.txt",sep="\t")
xx<-str_replace(colnames(xcell),"\\.","-")
xx<-str_replace(xx,"\\.","-")
xx<-str_replace(xx,"\\.","-")

load("Data/PAAD/PAAD_RepertoireResults_clones.Rdata")
RSEM_PAAD<-xcell[,match(substr(as.character(PAAD_repertoire_diversity$TCGA_sample),1,15),xx)]
rownames(RSEM_PAAD)<-xcell$sample

##Annotation file
genecode<-read.table("~/TCGA-Immune/Data/xCELL/gencode.v23.annotation.gene.probemap.txt",header = T)

id.gene<-match(rownames(RSEM_PAAD),genecode$id)
RSEM_PAAD$Gene<-genecode$gene[id.gene]
write.csv(RSEM_PAAD,"~/TCGA-Immune/Data/xCELL/PAAD_rsem.csv")

###Run xcell tool to obtain the de-convolution
#http://xcell.ucsf.edu/