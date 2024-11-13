#!/usr/bin/env Rscript 

# Set the library path to the system library path within the Singularity container
#lib_path <- "/home/p717a/R/x86_64-pc-linux-gnu-library/4.3"
#if (!dir.exists(lib_path)) {
#  stop("Library path does not exist:", lib_path)
#}
#.libPaths(lib_path)

#pacman::p_load(
#Seurat, ggplot2, patchwork, tidyverse, hdf5r, ggsci, 
#celldex, RColorBrewer, SingleCellExperiment, glmGamPoi, 
#reticulate, cowplot, viridis, pheatmap, scran, SingleR, 
#BiocParallel, DoubletFinder,presto,argparse, ggrepel) 

Packages <- c("VGAM", "R.utils", "Rfast2", "ape", "mixtools", "spatstat.explore", 
              "spatstat.geom", "hdf5r", "cowplot", "pacman", "patchwork", "pheatmap", "ggplot2", 
              "ggrepel", "Matrix", "remotes", "Seurat", "ggsci", 
              "viridis", "dplyr", "BiocManager", "multtest", "S4Vectors", "SummarizedExperiment", 
              "SingleCellExperiment", "DESeq2", "BiocGenerics", "GenomicRanges", "presto", 
              "celldex", "Biobase", "limma", "glmGamPoi", "SingleR", "SeuratObject", "scran")


lapply(Packages, library, character.only = TRUE)


args <- commandArgs(trailingOnly = TRUE)
h5_file <- args[1]
output_prefix <- gsub('.h5', '', h5_file)

#seuratobject

adj_matrix <- Read10X_h5(h5_file,use.names = T) 
seurat_obj <- CreateSeuratObject(counts = adj_matrix,project = 'TISCH2', min.cells = 3, min.features = 200)
adj.matrix <- NULL
seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT')
seurat_obj[['percent.rb']] <- PercentageFeatureSet(seurat_obj, pattern = '^RP[SL]')
plot1<- VlnPlot(object = seurat_obj, features = c('nFeature_RNA','nCount_RNA','percent.mt','percent.rb'), ncol = 4)
plot2<-FeatureScatter(object = seurat_obj, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
plot3<-FeatureScatter(object = seurat_obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
plot4<-FeatureScatter(object = seurat_obj, feature1 = 'nCount_RNA', feature2 = 'percent.rb')
plot5<-FeatureScatter(object = seurat_obj, feature1 = 'percent.rb', feature2 = 'percent.mt')
output_prefix <- gsub('.h5', '', h5_file)
ggsave(filename = paste0(output_prefix, '_vlnplot.png'), plot = plot1)
ggsave(filename = paste0(output_prefix, '_featurescatter1.png'), plot = plot2)
ggsave(filename = paste0(output_prefix, '_featurescatter2.png'), plot = plot3)
ggsave(filename = paste0(output_prefix, '_featurescatter3.png'), plot = plot4)
ggsave(filename = paste0(output_prefix, '_featurescatter4.png'), plot = plot5)
saveRDS(seurat_obj, file = paste0(output_prefix, '_seurat_results.rds'))

#normalization

seurat_obj <- SCTransform(object = seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
top10 <- VariableFeatures(object = seurat_obj, 10)

#PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
plot8 <- DimPlot(object = seurat_obj, reduction = "pca", label = TRUE)
plot9<- ElbowPlot(object = seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, verbose = F)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:20, verbose = F, check_duplicates = FALSE)
plot10 <- DimPlot(object = seurat_obj, label.size = 4,repel = T,label = T)

output_prefix <- gsub('.h5', '', h5_file)
ggsave(filename = paste0(output_prefix, "_PCA_dimplot.png"), plot = plot8)
ggsave(filename = paste0(output_prefix, "_PCA_elbowplot.png"), plot = plot9)
ggsave(filename = paste0(output_prefix, "_DimPlot_UMAP.png"), plot = plot10)

#cellcycle
#cc.genes.updated.2019
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
#assign Cell Score
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes)
plot11<- FeaturePlot(object = seurat_obj,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
seurat_obj <- SCTransform(seurat_obj, method = "glmGamPoi", vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)
seurat_obj <- RunPCA(seurat_obj, verbose = F)
seurat_obj <- RunPCA(seurat_obj, verbose = F)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = F)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = F)
seurat_obj <- FindClusters(seurat_obj, verbose = F)
table(seurat_obj[[]]$seurat_clusters)
plot12<- DimPlot(object = seurat_obj, label.size = 4,repel = T,label = T)


output_prefix <- gsub('.h5', '', h5_file)
ggsave(filename = paste0(output_prefix, "_feature_cellcycle.png"), plot = plot11)
ggsave(filename = paste0(output_prefix, "_DimPlot_cellcycle+percentmt.png"), plot = plot12)

#markers
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
all.markers <- FindAllMarkers(seurat_obj, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)

top3_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers

all.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 3) %>%
    ungroup() -> top3

plot13 <- DoHeatmap(object = seurat_obj, features = top3$gene) + NoLegend()


#write csv
output_prefix <- gsub('.h5', '', h5_file)
ggsave(filename = paste0(output_prefix, '_heatmap.png'), plot = plot13)
write.csv(all.markers, file = paste0(output_prefix, '_output.csv'), row.names = FALSE)
saveRDS(seurat_obj, file = paste0(output_prefix, '_seurat_results_final.rds'))
