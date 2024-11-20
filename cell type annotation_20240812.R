#---------------------------------------------------------------------------------------#
library(tidyverse)
library(dplyr)
library(ggplot2)

library(Seurat)
library(patchwork)

#-------------------------------------------------------------------------------#
merged_aa <- readRDS("20240509_merged_forebrain_integrated_cca_harmony.rds")

#Subset on value in the object meta data
subset_aa <- subset(x = merged_aa, subset = (orig.ident == "E10_forebrain_2"|orig.ident == "E11_forebrain_1"|orig.ident == "E12.5_forebraindorsal_1"|orig.ident == "E13.5_forebraindorsal_1"| orig.ident == "E15.5_forebraindorsal_1"))
rm(merged_aa)

#-------------------------------------------------------------------------------#
p1 <- DimPlot(subset_aa,reduction = "umap.cca",group.by = c("source", "cca_clusters"),
              combine = FALSE, label.size = 1)

p2 <- DimPlot(subset_aa,reduction = "umap.harmony",group.by = c("source", "harmony_clusters"),
              combine = FALSE, label.size = 1)

wrap_plots(c(p1,p2), ncol = 2, byrow = F)

#-------------------------------------------------------------------------------#
DimPlot(subset_aa, split.by = "source", reduction = "umap.cca", label = T, label.size = 2.5) + NoLegend()
DimPlot(subset_aa, split.by = "source", reduction = "umap.harmony", label = T, label.size = 2.5) + NoLegend()

FeaturePlot(subset_aa, split.by = "source", features = "Dlx5", reduction = "umap.cca", cols = c("grey","red"))

DimPlot(subset_aa,reduction = "umap.cca", group.by ="cca_clusters", label = T, label.size = 4) + NoLegend()

#---------JoinLayers--FindAllMarkers--------------------------------------------#
aa <- subset_aa
aa <- JoinLayers(aa)
rm(subset_aa)

# Set the identity of your cells to the desired column
Idents(aa) <- "cca_clusters"  #choose the method of integration

aa_markers <- FindAllMarkers(aa)
#saveRDS(aa, "aa.rds")
#aa <- readRDS("aa.rds")

write.csv(aa_markers,"aa_markers.csv", row.names = F)
#saveRDS(aa_markers, "aa_markers.rds")
#aa_markers <- readRDS("aa_markers.rds")

aa_markers_filter <- aa_markers %>%  group_by(cluster) %>% filter(p_val_adj < 0.05)

aa_markers_filter_up <- aa_markers %>%  group_by(cluster) %>% dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05)
aa_markers_filter_down <- aa_markers %>%  group_by(cluster) %>% dplyr::filter(avg_log2FC < 1, p_val_adj < 0.05)

write.csv(aa_markers_filter,"aa_markers_filter.csv", row.names = F)
write.csv(aa_markers_filter_up,"aa_markers_filter_up.csv", row.names = F)
write.csv(aa_markers_filter_down,"aa_markers_filter_down.csv", row.names = F)

#-------visualizations of marker feature expression-----------------------------#
VlnPlot(aa, features = c("Sox2", "Pax6","Tbr1", "Eomes"), ncol = 2) # Tbr2=Eomes
VlnPlot(aa, features = c("Sox2", "Pax6","Tbr1", "Eomes"), group.by = "source", ncol = 2) 
VlnPlot(aa, features = c("Sox2", "Pax6","Tbr1", "Eomes"), group.by = "harmony_clusters", ncol = 2) 

VlnPlot(aa, features = c("Sox2","Pax6","Tbr1")) # Tbr2=Eomes

FeaturePlot(aa, features = c("Espn", "Ap1m2","Gm16136","Eno3"), reduction = "umap.cca", cols = c("grey","red"))

aa_markers_filter_top20  <- aa_markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1,p_val_adj < 0.05) %>%
  slice_head(n = 20) %>% ungroup() 
write.csv(aa_markers_filter_top20,"aa_markers_filter_top20.csv", row.names = F)

DoHeatmap(aa, features = top10$gene) + NoLegend()

#---harmony_clusters----------------------------------------------------------------------------#
Idents(aa) <- "harmony_clusters"

aa_markers_har <- FindAllMarkers(aa)
#aa_har <- readRDS("aa_har.rds")

top10_har  <- aa_markers_har %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup() 
DoHeatmap(aa_har, features = top10_har$gene) + NoLegend()

#--cell type annotation---------------------------------------------------------#
##Marker-based annotation
##Reference-based annotation[singleR]

#https://biostatsquid.com/singler-tutorial/
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scRNAseq")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scran")
BiocManager::install("scuttle")
BiocManager::install("pheatmap")

library(SingleR)
library(celldex)
library(scRNAseq)
library(SingleCellExperiment)
library(scran)
library(pheatmap)

#----#Reference-based annotation[singleR]--Mouse--------------------------------#
###1. Get normalised or raw counts
norm_counts <- LayerData(aa, assay = "RNA", layer = 'data')

###2.Get the reference dataset
ref <- celldex::MouseRNAseqData()
unique(ref$label.main)
unique(ref$label.fine)

###3.Run SingleR
ct_ann <- SingleR(test = norm_counts, ref = ref, labels = ref$label.main, de.method = 'wilcox')

###4.Inspect quality of the predictions
plotScoreHeatmap(ct_ann) 
plotDeltaDistribution(ct_ann, ncol = 4, dots.on.top = FALSE)

###5.Add SingleR predictions to Seurat object
rownames(ct_ann)[1:5] # make sure you have cell IDs
bb <- AddMetaData(aa, ct_ann$pruned.labels, col.name = 'SingleR_HCA') 

# Visualise them on the UMAP
bb <- SetIdent(bb, value = "SingleR_HCA")
DimPlot(bb, reduction = "umap.cca", label = T , repel = T, label.size = 4) #+ NoLegend()

#demo was finished, then remove bb object
rm(bb, ct_ann, norm_counts, ref)

#--Marker-based annotation------------------------------------------------------#
#Marker-based annotation
FeaturePlot(aa, features = c("Sox2","Pax6","Hes5"),label = T, label.size = 4, reduction = "umap.cca", 
            cols = c("grey","red"),ncol = 3) + (plot_annotation(title = "Markers of Apical Progenitors", 
                                                                theme=theme(plot.title=element_text(hjust=0.5,size = 20))))

FeaturePlot(aa, features = c("Eomes","Neurog2","Btg2"),label = T, label.size = 4, reduction = "umap.cca", 
            cols = c("grey","red"),ncol = 3) + (plot_annotation(title = "Markers of Intermediate Progenitors", 
                                                                theme=theme(plot.title=element_text(hjust=0.5,size = 20))))

FeaturePlot(aa, features = c("Neurod2","Tubb3","Neurod6"),label = T, label.size = 4, reduction = "umap.cca", 
            cols = c("grey","red"),ncol = 3) + (plot_annotation(title = "Markers of Projection neuron", 
                                                                theme=theme(plot.title=element_text(hjust=0.5,size = 20))))

VlnPlot(aa, features =c("Neurod2","Tubb3")) + NoLegend()


#---------Rename cluster-------------------------------------------------------#
cluster_names <- c("Callosal neurons", "Apical progenitors", "Unknow_2", "Unknow_3", "Apical progenitors","Projection neuron",
                   "Apical progenitors","Intermediate progenitors", "Putative near-projecting neurons","MGE-derived interneurons",
                   "Callosal neurons",
                   "Inhibitory interneurons","Non-neuronal","Corticofugal neurons","Unknow_14","Inhibitory interneurons",
                   "Non-neuronal","Non-neuronal","Non-neuronal","Unknow_19", "Non-neuronal",
                   "Non-neuronal","Ependymocytes")

names(cluster_names) <- levels(aa)  

aa <- RenameIdents(object = aa, cluster_names) #assign condition
DimPlot(aa, reduction = "umap.cca", label = TRUE, repel = T, label.size = 4) #+ NoLegend()

#---------Remove cluster--------------------------------------------------------#
bb <- aa

remove <- WhichCells(object = bb, idents = "Non-neuronal",invert = T) # choose not non-neuron, invert = T
bb <- subset(x = bb, cells = remove)
DimPlot(bb, reduction = "umap.cca", label = T, repel = T, label.size = 4) #+ NoLegend()

or

p <- DimPlot(bb, reduction = "umap.cca")+ NoLegend()
LabelClusters(p, id = "ident",  repel = T, fontface = "bold", color = "black") #change text size

#------trajectory analysis---monocle3-------------------------------------------#
#https://cole-trapnell-lab.github.io/monocle3/docs/installation/

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")#3.19=1.30.20

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'lme4', 'S4Vectors', 
                       'SingleCellExperiment','terra', 'ggrastr'), force = TRUE)

install.packages("devtools") 

#https://cran.r-project.org/bin/windows/Rtools/   #for_windown run Rtools
devtools::install_github('cole-trapnell-lab/monocle3')
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)

library(monocle3)
library(SeuratWrappers)

#ref:https://rpubs.com/mahima_bose/Seurat_and_Monocle3_p

#Building trajectories with Monocle 3
cc <- as.cell_data_set(bb) #Seurat object-> CellDataSet object 
head(colData(cc)) #Get cell metadata (show format)

fData(cc) #Get feature metadata;fdata:Converts other data classes into fdata class. 
rownames(fData(cc))[1:10] # gene name, make a gene id 
fData(cc)$gene_id <- rownames(fData(cc))
head(fData(cc))
head(counts(cc)) #Get counts

#Retrieve clustering information from Surat object
#Assign partitions (cluster) to monocle (convert Seurat object to monocle object)
partitions <- c(rep(1, length(cc@colData@rownames)))
names(partitions) <- cc@colData@rownames
partitions <- as.factor(partitions)
partitions
cc@clusters@listData[["UMAP"]][["partitions"]] <- partitions 

#Assign cluster information
cluster <- bb@active.ident
cc@clusters@listData[["UMAP"]][["clusters"]] <- cluster

#Assign UMAP coordinates to monocle
cc@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- bb@reductions$umap.cca@cell.embeddings
cluster.before.traj <-plot_cells(cc, color_cells_by = "cluster", label_groups_by_cluster = F, group_label_size = 5) + NoLegend() 

#Learn Trajectory
cc <- learn_graph(cc, use_partition = F) #connect dots between cluster, depends on cluster number, assign start point, e.g. AP start)
plot_cells(cc, color_cells_by = "cluster", label_groups_by_cluster = F,label_branch_points = T, label_roots = T, label_leaves = F,group_label_size = 5)

#Order cells in pseudotime 
cc <- order_cells(cc, reduction_method = "UMAP", root_cells = NULL)

plot_cells(cc, color_cells_by = "pseudotime", label_groups_by_cluster = T,label_branch_points = T, label_roots = F, label_leaves = F) + scale_color_gradient(low = "mediumpurple1", high = "khaki") + labs(color = "pseudotime")
