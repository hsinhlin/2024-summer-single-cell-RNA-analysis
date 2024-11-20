#------trajectory analysis---monocle3-------------------------------------------#
#https://cole-trapnell-lab.github.io/monocle3/docs/installation/

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")# 3.19=1.30.20

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'lme4', 'S4Vectors', 'SingleCellExperiment','terra', 'ggrastr'), force = TRUE)

install.packages("devtools") 

#https://cran.r-project.org/bin/windows/Rtools/   #for_windown run Rtools
devtools::install_github('cole-trapnell-lab/monocle3')
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)

library(monocle3)
library(SeuratWrappers)

#----------Convert  Seurat into CellDataSet object----------
cds_obj <- SeuratWrappers::as.cell_data_set(aa)
set.seed(42)
# Perform dimensionality reduction using UMAP
cds_obj <- preprocess_cds(cds_obj, num_dim = 30)  # Adjust num_dim as necessary
cds_obj <- reduce_dimension(cds_obj, reduction_method = "UMAP", preprocess_method = "PCA")
cds_obj <- cluster_cells(cds_obj)
cds_obj <- learn_graph(cds_obj, use_partition=FALSE, close_loop=FALSE)
p1 <- plot_cells(cds_obj, color_cells_by="seurat_clusters",  group_label_size=4, graph_label_size=3, label_cell_groups=FALSE, label_principal_points=TRUE, label_groups_by_cluster=FALSE)
scale_color_manual(values = my_colors_15)

#pathway enrichment test 
