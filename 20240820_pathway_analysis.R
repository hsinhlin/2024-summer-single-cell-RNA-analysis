#-----------pathway analysis---------------------------------------------------#
#KEGG database analyze cluster 2 -> 對應去define cluster (include gene function)
#------- local KEGG database---------------------------------------------------#
remotes::install_github("YuLab-SMU/createKEGGdb") #import the library and create a KEGG database locally 
createKEGGdb::create_kegg_db(species) #You will get KEGG.db_1.0.tar.gz file in your working directory

install.packages("KEGG.db_1.0.tar.gz", repos=NULL,type="source")
species <-c("ath","hsa","mmu","rno","dre","dme","cel")

#BiocManager::install("org.Mm.eg.db") #human:org.Hs.eg.db
library(createKEGGdb)
library(KEGG.db)
library(clusterProfiler)

library(dplyr)
library(tidyr)
library(tidyverse)

library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_scout()

library(EnsDb.Mmusculus.v79)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)
library(KEGGREST)

###KEGG
#-----input gene list----------------------------------------------------------#
aa <- read.csv("e:/R_analysis/8_1.scRNA/aa_markers_20240820.csv")
aa_a <- aa %>% filter (cluster == 2, p_val_adj < 0.05) #choose cluster 2 only, add condition if needed

#aa_a <- aa  %>% filter (cluster == 2 ,avg_log2FC >= 1)

#-------prepare for pathway analysis-------------------------------------------#
aa_a <- aa_a  %>% select (gene)
aa_a = bitr(aa_a$gene, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Mm.eg.db")

#-----input defult of gene list-------------------------------------------------#
data(geneList, package="DOSE") #Disease Ontology Semantic and Enrichment analysis (run first time)
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = "org.Hs.eg.db")
head(gene.df,10)

#----KEGG- output---------------------------------------------------------------#
kk <- enrichKEGG(gene = aa_a$ENTREZID, organism = 'mmu', pvalueCutoff = 0.05,  use_internal_data=T)
head(kk,10)

dotplot(kk,showCategory=15, title = "cluster_2[p<0.05,3031 genes]")

k1 <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
k1 = as.data.frame(k1) 

write.csv(k1,"k1.csv")

#ShinyGO 0.80
#http://bioinformatics.sdstate.edu/go/