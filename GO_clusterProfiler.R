if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("Category")
BiocManager::install("AnnotationHub")
install.packages("tidyverse")
# Load library and search for the organism
library(clusterProfiler)
library(pathview)
library(tidyverse)
library(enrichplot)

#Build the Salmon Data Base and import the gene universe (background genes from Salmon genome)
library(Category)
library(AnnotationHub)
hub <- AnnotationHub()
yes
query(hub, c("salmo salar","orgdb"))
salmodb <- hub[["AH72154"]]
DatPkgFactory(salmodb)
columns(salmodb)

ovary_gene_universe <- read.table("/Users/moh034/Desktop/Ovary_gene_universe.txt", header=T, sep="\t")
head(ovary_gene_universe)
#Ovary hypermethylated genes
Ovary_hyperM_genic <- read.table("/Users/moh034/Desktop/ovary_hyper_GBM_clusterprofiler.txt", header=T, sep="\t")
dim(Ovary_hyperM_genic)
head(Ovary_hyperM_genic)
rownames(Ovary_hyperM_genic)<- NULL
Ovary_hyperM_genic_FC <- Ovary_hyperM_genic$meanMethy
names(Ovary_hyperM_genic_FC) <- Ovary_hyperM_genic$geneID
head(Ovary_hyperM_genic_FC)
Ovary_hyperM_genic_FC <- sort(Ovary_hyperM_genic_FC, decreasing = TRUE)
ego_Ovary_hyperM_genic <- enrichGO(gene = names(Ovary_hyperM_genic_FC), universe = as.character(ovary_gene_universe$ENTREZID), keyType = 'ENTREZID', OrgDb = salmodb, ont="ALL", pAdjustMethod = "fdr",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_Ovary_hyperM_genic_Summary <- data.frame(ego_Ovary_hyperM_genic)
write.csv(ego_Ovary_hyperM_genic_Summary, file = "/Users/moh034/Desktop/ego_Ovary_hyperM_genic.csv", row.names=FALSE)
dotplot(ego_Ovary_hyperM_genic, showCategory=20, font.size = 10, orderBy = "x")

#Ovary 148 Upregulated and Hypermethylated genes 
common148 <- read.table("/Users/moh034/Desktop/148OvaryUpHyper_common_genes.txt", header=T, sep="\t")
dim(common148)
head(common148)
rownames(common148)<- NULL
common148_FC <- common148$logFC
names(common148_FC) <- common148$geneID
head(common148_FC)

common148_FC <- sort(common148_FC, decreasing = TRUE)
ego_common148 <- enrichGO(gene = names(common148_FC), universe = as.character(ovary_gene_universe$ENTREZID), keyType = 'ENTREZID', OrgDb = salmodb, ont="ALL", pAdjustMethod = "fdr",pvalueCutoff  = 0.05,qvalueCutoff  = 0.05)
ego_common148_Summary <- data.frame(ego_common148)
write.csv(ego_common148_Summary, file = "/Users/moh034/Desktop/ego_common148.csv", row.names=FALSE)
dotplot(ego_common148, showCategory=20, font.size = 10)
dev.off()
