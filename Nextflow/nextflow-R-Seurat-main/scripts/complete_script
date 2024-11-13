#!/usr/bin/env Rscript
pacman::p_load(
Seurat, ggplot2, patchwork, tidyverse, hdf5r, ggsci, 
celldex, RColorBrewer, SingleCellExperiment, glmGamPoi, 
reticulate, cowplot, viridis, pheatmap, scran, SingleR, 
BiocParallel, DoubletFinder,presto,argparse, ggrepel) 

args <- commandArgs(trailingOnly = TRUE)
h5_file <- args[1]

adj_matrix <- Read10X_h5(h5_file,use.names = T) 
seurat_obj <- CreateSeuratObject(counts = adj_matrix,project = 'TISCH2', min.cells = 3, min.features = 200)
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

seurat_obj <- SCTransform(object = seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
#seurat_obj <- NormalizeData(seurat_obj, normalization.method = 'LogNormalize', scale.factor = 10000)
#seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = 'vst', nfeatures = 2000)


seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, verbose=FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)


plot6 <- DimPlot(object = seurat_obj, label = TRUE)
ggsave(filename = paste0(output_prefix, "_dimplot.png"), plot = plot6)

top10 <- VariableFeatures(object = seurat_obj, 10)
plot7 <- VariableFeaturePlot(object = seurat_obj)
LabelPoints(plot = plot7, points = top10, repel = TRUE, xnudge = 0.2, ynudge = 0.2) 
ggsave(filename = paste0(output_prefix, '_top10_variablefeatures.png'), plot = plot7, width = 8, height = 6)


all.markers <- FindAllMarkers(seurat_obj, only.pos = T)
all.markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1)

all.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

plot8 <- DoHeatmap(object = seurat_obj, features = top10$gene) + NoLegend()
ggsave(filename = paste0(output_prefix, '_heatmap.png'), plot = plot8)

#write csv
write.csv(all.markers, file = paste0(output_prefix, '_output.csv'), row.names = FALSE)


cc.genes.updated.2019
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
#assign Cell Score
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE,use.synonyms = TRUE)
seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))

#Dimplot
plot9 <-DimPlot(object = seurat_obj)
#FeaturePlot
plot10 <-FeaturePlot(object = seurat_obj,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size=10))
#ViolinPlot
plot11 <-VlnPlot(seurat_obj,features = c("S.Score","G2M.Score")) & theme(plot.title = element_text(size=10))
ggsave(paste0(output_prefix, '_dimplot_cell_cycle.png'), plot = plot9, width = 11, height = 8.5, dpi = 300)
ggsave(paste0(output_prefix, '_feature_cell_cycle.png'), plot = plot10, width = 11, height = 8.5, dpi = 300)
ggsave(paste0(output_prefix, '_violin_cell_cycle.png'), plot = plot11, width = 11, height = 8.5, dpi = 300)

# CELL Type Annotation
#Get Data From reference database
monaco.ref <- celldex::MonacoImmuneData()
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
dice.ref <- celldex::DatabaseImmuneCellExpressionData()

sce <- as.SingleCellExperiment(DietSeurat(seurat_obj))

monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)

#Add annotations to the metadata
seurat_obj@meta.data$monaco.main <- monaco.main$pruned.labels
seurat_obj@meta.data$monaco.fine <- monaco.fine$pruned.labels
seurat_obj@meta.data$hpca.main   <- hpca.main$pruned.labels
seurat_obj@meta.data$hpca.fine   <- hpca.fine$pruned.labels
seurat_obj@meta.data$dice.main   <- dice.main$pruned.labels
seurat_obj@meta.data$dice.fine   <- dice.fine$pruned.labels

# Original seurat_obj
original_seurat_obj <- seurat_obj

# Plot for "monaco.fine"
seurat_obj_monaco <- SetIdent(original_seurat_obj, value = "monaco.fine")
plot12 <- DimPlot(object = seurat_obj_monaco, label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
ggsave(filename = paste0(output_prefix, '_monaco_dimplot.png'), plot = plot12, width = 11, height = 8.5, dpi = 300)

# Plot for "hpca.fine"
seurat_obj_hpca <- SetIdent(original_seurat_obj, value = "hpca.fine")
plot13 <- DimPlot(object = seurat_obj_hpca, label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
ggsave(filename = paste0(output_prefix, '_hpca_dimplot.png'), plot = plot13, width = 11, height = 8.5, dpi = 300)

# Plot for "dice.fine"
seurat_obj_dice <- SetIdent(original_seurat_obj, value = "dice.fine")
plot14 <- DimPlot(object = seurat_obj_dice, label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
ggsave(filename = paste0(output_prefix, '_dice_dimplot.png'), plot = plot14, width = 11, height = 8.5, dpi = 300)





saveRDS(seurat_obj, file = paste0(output_prefix, '_seurat_results_final.rds'))

