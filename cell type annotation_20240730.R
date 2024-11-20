library(tidyverse)
library(dplyr)
library(ggplot2)

library(Seurat)
library(patchwork)
library(SeuratWrappers)

#-------------------------------------------------------------------------------#
merged_aa <- readRDS("20240509_merged_forebrain_integrated_cca_harmony.rds")

#Subset on value in the object meta data
subset_aa <- subset(x = merged_aa, subset = (orig.ident == "E10_forebrain_2"|orig.ident == "E11_forebrain_1"|orig.ident == "E12.5_forebraindorsal_1"|orig.ident == "E13.5_forebraindorsal_1"| orig.ident == "E15.5_forebraindorsal_1"))
rm(merged_aa)

#-------------------------------------------------------------------------------#
p1 <- DimPlot(subset_aa,reduction = "umap.cca",group.by = c("source", "cca_clusters"), combine = FALSE, label.size = 1)

p2 <- DimPlot(subset_aa,reduction = "umap.harmony",group.by = c("source", "harmony_clusters"), combine = FALSE, label.size = 1)

wrap_plots(c(p1,p2), ncol = 2, byrow = F)

#-------------------------------------------------------------------------------#
DimPlot(subset_aa, split.by = "source", reduction = "umap.cca", label = T, label.size = 2.5) + NoLegend()
DimPlot(subset_aa, split.by = "source", reduction = "umap.harmony", label = T, label.size = 2.5) + NoLegend()

FeaturePlot(subset_aa, split.by = "source", features = "Sox2", reduction = "umap.cca", cols = c("grey","red"))

DimPlot(subset_aa,reduction = "umap.cca", group.by ="cca_clusters", label = T, label.size = 4) + NoLegend()

#---------JoinLayers--FindAllMarkers--------------------------------------------#
aa <- subset_aa
aa <- JoinLayers(aa)

# Set the identity of your cells to the desired column
Idents(aa = aa) <- "cca_clusters"  #choose the method of integration

aa_markers <- FindAllMarkers(aa)
#saveRDS(aa, "aa.rds")
#aa <- readRDS("aa.rds")

write.csv(aa_markers,"aa_markers.csv", row.names = F)
#saveRDS(aa_markers, "aa_markers.rds")
#aa_markers <- readRDS("aa_markers.rds")
aa_markers_filter <- aa_markers %>%  group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)
write.csv(aa_markers_filter,"aa_markers_filterr.csv", row.names = F)

#-------visualizations of marker feature expression-----------------------------#
VlnPlot(aa, features = c("Sox2", "Pax6","Tbr1", "Eomes"), ncol = 2) # Tbr2=Eomes
VlnPlot(aa, features = c("Sox2", "Pax6","Tbr1", "Eomes"), group.by = "source", ncol = 2) 
VlnPlot(aa, features = c("Sox2", "Pax6","Tbr1", "Eomes"), group.by = "harmony_clusters", ncol = 2) 

VlnPlot(aa, features = c("Sox2","Pax6","Tbr1")) # Tbr2=Eomes

FeaturePlot(aa, features = c("Sox2", "Pax6","Tbr1","Eomes"), reduction = "umap.cca", split.by = "source", cols = c("grey","red")) # label = T, label.size = 2

top10  <- aa_markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup()  #pos = T
DoHeatmap(aa, features = top10$gene) + NoLegend()

#---harmony_clusters----------------------------------------------------------------------------#
Idents(aa = aa) <- "harmony_clusters"

aa_markers_har <- FindAllMarkers(aa)
#aa_har <- readRDS("aa_har.rds")

top10_har  <- aa_markers_har %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup() 
DoHeatmap(aa_har, features = top10_har$gene) + NoLegend()

#--cell type annotation---------------------------------------------------------#
##Mark-based annotation
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
ref <- celldex::MouseRNAseqData()#celldex library
unique(ref$label.main)  
####shows the annotated cell unique Function: This function extracts the unique cell type labels from the label.main column.
unique(ref$label.fine)
#Fine-Grained Labels: Provide more detailed cell type annotations, useful for distinguishing between similar cell types.

###3.Run SingleR
ct_ann <- SingleR(test = norm_counts, ref = ref, labels = ref$label.main, de.method = 'wilcox') #ref can change
#de.method: The differential expression method used for comparing the test and reference datasets. 'wilcox' refers to the Wilcoxon rank-sum test.


###4.Inspect quality of the predictions
plotScoreHeatmap(ct_ann) #compare with reference
plotDeltaDistribution(ct_ann, ncol = 4, dots.on.top = FALSE)

###5.Add SingleR predictions to Seurat object
rownames(ct_ann)[1:5] # make sure you have cell IDs
bb <- AddMetaData(aa, ct_ann$pruned.labels, col.name = 'SingleR_HCA') #insert cell name column

# Visualise them on the UMAP
bb <- SetIdent(bb, value = "SingleR_HCA")
DimPlot(bb, reduction = "umap.cca", label = T , repel = T, label.size = 3.5) + NoLegend()


#----#Reference-based annotation[singleR]--Human--------------------------------#
###1. Get normalised or raw counts
sce <- scRNAseq::KotliarovPBMCData(mode = c('rna'))
seu <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
#The “Feature names cannot have underscores …” warning message can be ignored. It’s caused by difference between dynverse and Seurat.
rm(sce)

seu <- NormalizeData(object = seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:10)

#saveRDS(seu, "seu_20240729.rds")
#seu <- readRDS("seu_20240729.rds")

raw_counts <- LayerData(seu, assay = "RNA", layer = 'counts') #
raw_counts[c('VIM', 'BCL2', 'TP53', 'CD4'),1:5]
norm_counts <- LayerData(seu, assay = "RNA", layer = 'data') #
norm_counts[c('VIM', 'BCL2', 'TP53', 'CD4'),1:5]

###2.Get the reference dataset
ref <- celldex::HumanPrimaryCellAtlasData()
unique(ref$label.main)
unique(ref$label.fine)

ref <- ref[,grepl('DC|B_cell|Neutrophils|T_cells|Monocyte|Erythroblast|
                 Macrophage|NK_cell|Platelets|Myelocyte', ref$label.main)]
unique(ref$label.main)

###3.Run SingleR
ct_ann <- SingleR(test = norm_counts, # we could also use sce or raw_counts
                  ref = ref, 
                  labels = ref$label.main,
                  de.method = 'wilcox')

###4.Inspect quality of the predictions
plotScoreHeatmap(ct_ann)
plotDeltaDistribution(ct_ann, ncol = 4, dots.on.top = FALSE)

###5.Add SingleR predictions to Seurat object
# Add to seurat object
rownames(ct_ann)[1:5] # make sure you have cell IDs
seu <- AddMetaData(seu, ct_ann$pruned.labels, col.name = 'SingleR_HCA')

# Visualise them on the UMAP
seu <- SetIdent(seu, value = "SingleR_HCA")
DimPlot(seu, reduction = "umap", label = T , repel = T, label.size = 3) + NoLegend()

#--Mark-based annotation---------------------------------------------------------#
#Mark-based annotation
Idents(aa = aa) <- "cca_clusters" 
p1 <- FeaturePlot(aa, features = "Sox2",label = T, label.size = 4,  reduction = "umap.cca", split.by = "source", cols = c("grey","red"))
p2 <- FeaturePlot(aa, features = "Pax6",label = T, label.size = 4, reduction = "umap.cca",  cols = c("grey","red"))
p3 <- FeaturePlot(aa, features = "Hes5",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
wrap_plots(p1,p2,p3)+(plot_annotation(title = "Markers of Apical Progenitors", caption = 'made with patchwork',
                                      theme=theme(plot.title=element_text(hjust=0.5,size = 20))))


p4 <- FeaturePlot(aa, features = "Eomes",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p5 <- FeaturePlot(aa, features = "Neurog2",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p6 <- FeaturePlot(aa, features = "Btg2",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
wrap_plots(p4,p5,p6)+(plot_annotation(title = "Markers of Intermediate Progenitors", caption = 'made with patchwork',
                                      theme=theme(plot.title=element_text(hjust=0.5,size = 20))))

p7 <- FeaturePlot(aa, features = "Tubb3",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p8 <- FeaturePlot(aa, features = "Neurog2",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p9 <- FeaturePlot(aa, features = "Neurod6",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
wrap_plots(p7,p8,p9)+(plot_annotation(title = "Markers of Projection Neurons", caption = 'made with patchwork',
                                      theme=theme(plot.title=element_text(hjust=0.5,size = 20))))

p10 <- FeaturePlot(aa, features = "Satb2",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p11 <- FeaturePlot(aa, features = "Cux2",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
wrap_plots(p10,p11)+(plot_annotation(title = "Markers of Callosal Neurons", caption = 'made with patchwork',
                                     theme=theme(plot.title=element_text(hjust=0.5,size = 20))))

p12 <- FeaturePlot(aa, features = "Fezf2",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p13 <- FeaturePlot(aa, features = "Tle4",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p14 <- FeaturePlot(aa, features = "Pcp4",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p15 <- FeaturePlot(aa, features = "Tcerg1l",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
w <- p12|p13|p14|p15
wrap_plots(w)+(plot_annotation(title = "Markers of Corticofugal Neurons", caption = 'made with patchwork',
                               theme=theme(plot.title=element_text(hjust=0.5,size = 20))))

p16 <- FeaturePlot(aa, features = "Sst",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p17 <- FeaturePlot(aa, features = "Npy",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p18 <- FeaturePlot(aa, features = "Lhx6",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p19 <- FeaturePlot(aa, features = "Nxph1",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p20 <- FeaturePlot(aa, features = "Nxph2",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
x <- p16|p17|p18|p19|p20
wrap_plots(x)+(plot_annotation(title = "Markers of MGE derived interneurons", caption = 'made with patchwork',
                               theme=theme(plot.title=element_text(hjust=0.5,size = 20))))

p21 <- FeaturePlot(aa, features = "Sst",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))

##Update from peiyi
p1 <- FeaturePlot(aa, features = "Sox2",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","red"))
p2 <- FeaturePlot(aa, features = "Pax6",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","mediumorchid2"))
p3 <- FeaturePlot(aa, features = "Hes5",label = T, label.size = 4, reduction = "umap.cca", cols = c("grey","royalblue"))
wrap_plots(p1,p2,p3)+(plot_annotation(title = "Markers of Intermediate Progenitors",
                                      theme=theme(plot.title=element_text(hjust=0.5,size = 20))))


FeaturePlot(aa, features = c("Sox2","Pax6","Hes5"),label = T, label.size = 4, reduction = "umap.cca", 
            cols = c("grey","red"),ncol = 3) + (plot_annotation(title = "Markers of Apical Progenitors", 
                                                                theme=theme(plot.title=element_text(hjust=0.5,size = 20))))

#--------
tab_number <- table(aa$group, aa$cca_clusters)
print(tab_number)

#---------Rename cluster-------------------------------------------------------#
cluster_names <- c("Callosal neurons", "Apical progenitors", "Unknow_2", "Unknow_3", "Apical progenitors","Projection neuron",
                   "Apical progenitors","Intermediate progenitors", "Putative near-projecting neurons","MGE-derived interneurons",
                   "Callosal neurons",
                   "Inhibitory interneurons","Non-neuronal","Corticofugal neurons","Unknow_14","Inhibitory interneurons",
                   "Non-neuronal","Non-neuronal","Non-neuronal","Unknow_19", "Non-neuronal",
                   "Non-neuronal","Ependymocytes")

names(cluster_names) <- levels(aa)

aa <- RenameIdents(object = aa, cluster_names)
DimPlot(aa, reduction = "umap.cca", label = TRUE, repel = T, label.size = 4) #+ NoLegend()

#---------Remove cluster--------------------------------------------------------#
bb <- aa

remove <- WhichCells(object = bb, idents = "Non-neuronal",invert = T)
bb <- subset(x = bb, cells = remove)

DimPlot(bb, reduction = "umap.cca", label = T, repel = T, label.size = 4) #+ NoLegend()

##-----write differebt clusters into seperate excel page----------------#
library(writexl)
library(openxlsx)

filtered_markers <- aa_markers %>% filter(cluster == "22")

write_xlsx(filtered_markers, path = "/Users/rusher/Desktop/R 0723/C0.xlsx")
workbook <- loadWorkbook("/Users/rusher/Desktop/R 0723/C0.xlsx")
addWorksheet(workbook, sheetName = "C22")
writeData(workbook, sheet = "C22", x = filtered_markers)
saveWorkbook(workbook, "/Users/rusher/Desktop/R 0723/C0.xlsx", overwrite = TRUE)
