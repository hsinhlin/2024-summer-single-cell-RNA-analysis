library(tidyverse)
library(dplyr)
library(ggplot2)

library(Seurat)
library(patchwork)
#library(SeuratWrappers)

#-------------------------------------------------------------------------------#
subset_aa <- readRDS("20240712_subset_aa.rds")

#Subset on value in the object meta data
subset_aa <- subset(x = subset_aa, 
                    subset = (
                      orig.ident == "E10_forebrain_2"|
                      orig.ident == "E11_forebrain_1"|
                      orig.ident == "E12.5_forebraindorsal_1"|
                      orig.ident == "E13.5_forebraindorsal_1"| 
                      orig.ident == "E15.5_forebraindorsal_1")
                    )                                                         

DimPlot(subset_aa, split.by = "source", reduction = "umap.cca", label = T, label.size = 3) + NoLegend()
#orig
#-------------------------------------------------------------------------------#
p1 <- DimPlot(subset_aa,reduction = "umap.cca",group.by = c("source", "cca_clusters"),
              combine = FALSE, label = T, label.size = 1,)

p2 <- DimPlot(subset_aa,reduction = "umap.harmony",group.by = c("source", "harmony_clusters"),
              combine = FALSE, label.size = 1)

wrap_plots(c(p1,p2), ncol = 2, byrow = F)

#---------JoinLayers--FindAllMarkers to Find DEG--------------------------------------------#
aa <- subset_aa 
aa <- JoinLayers(aa) #MUST join layer!!! #JoinLayers: combine batch effect method e.g. cca, harmony 

# Set the identity of your cells to the desired column
Idents(aa = aa) <- "cca_clusters"  #choose the method of integration #utilize 23 cluster to analyze, not 5 samples 

aa_markers <- FindAllMarkers(aa) #long time
#saveRDS(aa, "aa.rds")
#aa <- readRDS("aa.rds")

write.csv(aa_markers,"aa_markers.csv", row.names = F) #write excel table, already have it, go to visualization

aa_markers_filter <- aa_markers %>%  group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
write.csv(aa_markers_filter,"aa_markers_filterr.csv", row.names = F)


#-------visualizations of marker feature expression-----------------------------
##group.by = cluster or sample E10 etc.  #idnt.
FeaturePlot(aa, features = c("Sox2", "Pax6","Tbr1","Eomes"), reduction = "umap.cca", cols = c("grey","red"))

#Violin plot
VlnPlot(aa, features = c("Sox2", "Pax6","Tbr1","Eomes"), group.by = "source", ncol = 2) 
VlnPlot(aa, features = c("Bcl11b", "Sox5", "Lbd2")) + NoLegend()


#Draw Heat Map (choose top 10 by log2FC per cluster)
aa_markers <- readRDS("aa_markers.rds")

filtered_markers <- aa_markers %>% filter(cluster == "0")

top10  <- filtered_markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup() #need aa_markers.rds
DoHeatmap(aa, features = top10$gene) + NoLegend()
#set top 10 by ourself, 橫cluster直gene

### Feature plot split by source -------- 
#OPC Oligodendrocyte precursor cells 
FeaturePlot(aa, features = c("Olig1", "Olig2", "Pdgfra"), split.by = "source", reduction = "umap.cca", label = T, label.size = 2, cols = c("grey","red"))

#Astrocyte
FeaturePlot(aa, features = c("Apoe", "Aldh1l1", "Slcla13", "Slc1a", "Gfap", "Sparcl1"), split.by = "source", reduction = "umap.cca", cols = c("grey","red"))  #no Slcla13, Slc1a

#Microglia
FeaturePlot(aa, features = c("Aif1", "Tmem119"), split.by = "source", reduction = "umap.cca", cols = c("grey","red"))

#RBC (hemoglobin gene) 
FeaturePlot(aa, features = c("Car2", "Hemgn"), split.by = "source", reduction = "umap.cca", cols = c("grey","red"))

#Endothelial cells
FeaturePlot(aa, features = c("Cldn5", "Mcam"), split.by = "source", reduction = "umap.cca", cols = c("grey","red"))

#Pericytes
FeaturePlot(aa, features = c("Cspg4", "Pdgfrb", "Rgs5", "Notch3"), split.by = "source", reduction = "umap.cca", cols = c("grey","red"), label = T, label.size = 3)

#Vascular& leptomeningeal cells VLMC
FeaturePlot(aa, features = c("Col1a1", "Vtn", "Lgals1"), split.by = "source", reduction = "umap.cca", cols = c("grey","red"))

#---harmony_cluster----------------------------------------------------------------------------#
aa <- subset_aa
aa <- JoinLayers(aa)

Idents(aa = aa) <- "harmony_clusters"

aa_markers_har <- FindAllMarkers(aa)
#aa_har <- readRDS("aa_har.rds")

VlnPlot(aa_har, features = c("Sox2", "Pax6","Tbr1", "Eomes"), ncol = 2) # Tbr2=Eomes
 
FeaturePlot(aa_har, features = c("Sox2", "Pax6","Tbr1","Eomes"),reduction = "umap.harmony", min.cutoff = 3 ,max.cutoff = 5, cols = c("grey","red"))

top10_har  <- aa_markers_har %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup() 
DoHeatmap(aa_har, features = top10_har$gene) + NoLegend()

#--cell type annotation---------------------------------------------------------#
##Mark-based annotation 
##Reference-based annotation[singleR]

#define cell cluster (eg. microglia )

library(tidyverse)
library(dplyr)
library(ggplot2)

library(Seurat)
library(patchwork)
#library(SeuratWrappers)


#-------------------------------------------------------------------------------#
merged_aa <- readRDS("20240509_merged_forebrain_integrated_cca_harmony.rds")

#Subset on value in the object meta data
subset_aa <- subset(x = merged_aa, subset = (orig.ident == "E10_forebrain_2"|orig.ident == "E11_forebrain_1"|orig.ident == "E12.5_forebraindorsal_1"|orig.ident == "E13.5_forebraindorsal_1"| orig.ident == "E15.5_forebraindorsal_1"))

DimPlot(subset_aa, split.by = "source", reduction = "umap.cca", label = T, label.size = 2.5) + NoLegend()

FeaturePlot(subset_aa, split.by = "source", features = "Sox2", reduction = "umap.cca", cols = c("grey","red"))
