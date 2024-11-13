#!/usr/bin/env Rscript

# Load required libraries
library(Seurat)
library(future.apply)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(ggsci)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(glmGamPoi)
library(cowplot)
library(pheatmap)
library(scran)
library(SingleR)
library(presto)
library(ggrepel)
library(argparse)
library(dplyr)

# Set up argument parsing
args <- commandArgs(trailingOnly = TRUE)
h5_file <- args[1]
output_prefix <- gsub('_expression\\.h5$', '', h5_file)
output_dir <- paste0(output_prefix, "_results")

# Create output directory
dir.create(output_dir, showWarnings = FALSE)

# Define a function to handle errors
handle_error <- function(e) {
  message("An error occurred: ", e$message)
  traceback()
  quit(status = 1)
}

# Main processing wrapped in a tryCatch block
tryCatch({

  message("Reading H5 file...")
  adj_matrix <- Read10X_h5(h5_file, use.names = TRUE)
  
  message("Creating Seurat object...")
  seurat_obj <- CreateSeuratObject(counts = adj_matrix, project = 'TISCH2', min.cells = 3, min.features = 200)
  
  message("Calculating mitochondrial and ribosomal percentages...")
  seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = '^MT')
  seurat_obj[['percent.rb']] <- PercentageFeatureSet(seurat_obj, pattern = '^RP[SL]')
  
  message("Generating quality control plots...")
  plot1 <- VlnPlot(object = seurat_obj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'percent.rb'), ncol = 4)
  plot2 <- FeatureScatter(object = seurat_obj, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
  plot3 <- FeatureScatter(object = seurat_obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
  plot4 <- FeatureScatter(object = seurat_obj, feature1 = 'nCount_RNA', feature2 = 'percent.rb')
  plot5 <- FeatureScatter(object = seurat_obj, feature1 = 'percent.rb', feature2 = 'percent.mt')
  
  ggsave(filename = file.path(output_dir, paste0(output_prefix, '_vlnplot.png')), plot = plot1)
  ggsave(filename = file.path(output_dir, paste0(output_prefix, '_featurescatter1.png')), plot = plot2)
  ggsave(filename = file.path(output_dir, paste0(output_prefix, '_featurescatter2.png')), plot = plot3)
  ggsave(filename = file.path(output_dir, paste0(output_prefix, '_featurescatter3.png')), plot = plot4)
  ggsave(filename = file.path(output_dir, paste0(output_prefix, '_featurescatter4.png')), plot = plot5)
  
  message("Performing SCTransform normalization...")
  seurat_obj <- SCTransform(object = seurat_obj, assay = "RNA", vars.to.regress = "percent.mt", verbose = FALSE)
  
  message("Running PCA...")
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  
  message("Running UMAP...")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  
  message("Finding neighbors...")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
  
  message("Finding clusters...")
  seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
  
  message("Generating dimension plot...")
  plot6 <- DimPlot(object = seurat_obj, label = TRUE)
  ggsave(filename = file.path(output_dir, paste0(output_prefix, "_dimplot.png")), plot = plot6)
  
  message("Identifying top variable features...")
  top10 <- VariableFeatures(object = seurat_obj)[1:10]
  plot7 <- VariableFeaturePlot(object = seurat_obj)
  LabelPoints(plot = plot7, points = top10, repel = TRUE, xnudge = 0.2, ynudge = 0.2)
  ggsave(filename = file.path(output_dir, paste0(output_prefix, '_top10_variablefeatures.png')), plot = plot7, width = 8, height = 6)
  
  message("Finding all markers...")
  all.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
  write.csv(all.markers, file = file.path(output_dir, paste0(output_prefix, '_output.csv')), row.names = FALSE)
  
  message("Filtering top markers for heatmap...")
  top10 <- all.markers %>%
    group_by(cluster) %>%
    filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  message("Generating heatmap...")
  plot8 <- DoHeatmap(object = seurat_obj, features = top10$gene) + NoLegend()
  ggsave(filename = file.path(output_dir, paste0(output_prefix, '_heatmap.png')), plot = plot8)
  
  message("Cell cycle scoring...")
  cc.genes.updated.2019
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, use.synonyms = TRUE)
  seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))
  
  plot9 <- DimPlot(object = seurat_obj)
  plot10 <- FeaturePlot(object = seurat_obj, features = c("S.Score", "G2M.Score"), label.size = 4, repel = TRUE, label = TRUE) + theme(plot.title = element_text(size = 10))
  plot11 <- VlnPlot(seurat_obj, features = c("S.Score", "G2M.Score")) + theme(plot.title = element_text(size = 10))
  
  ggsave(file.path(output_dir, paste0(output_prefix, '_dimplot_cell_cycle.png')), plot = plot9, width = 11, height = 8.5, dpi = 300)
  ggsave(file.path(output_dir, paste0(output_prefix, '_feature_cell_cycle.png')), plot = plot10, width = 11, height = 8.5, dpi = 300)
  ggsave(file.path(output_dir, paste0(output_prefix, '_violin_cell_cycle.png')), plot = plot11, width = 11, height = 8.5, dpi = 300)
  
#  cache_file <- function(filename, compute_function, ...) {
#  if (file.exists(filename)) {
#    message(paste("Loading from cache:", filename))
#    readRDS(filename)
#  } else {
#    result <- compute_function(...)
#    saveRDS(result, filename)
#    result
#  }
#}
#
#message("Cell type annotation...")
#
## Define paths for cached results
#monaco_ref_cache <- file.path(output_dir, "monaco_ref.rds")
#hpca_ref_cache <- file.path(output_dir, "hpca_ref.rds")
#dice_ref_cache <- file.path(output_dir, "dice_ref.rds")
#sce_cache <- file.path(output_dir, "sce.rds")
#
## Load or compute references
#monaco.ref <- cache_file(monaco_ref_cache, celldex::MonacoImmuneData)
#hpca.ref <- cache_file(hpca_ref_cache, celldex::HumanPrimaryCellAtlasData)
#dice.ref <- cache_file(dice_ref_cache, celldex::DatabaseImmuneCellExpressionData)
#
## Convert Seurat object to SingleCellExperiment and cache
#sce <- cache_file(sce_cache, as.SingleCellExperiment, DietSeurat(seurat_obj))
#
## Annotation cache file paths
#monaco_main_cache <- file.path(output_dir, "monaco_main.rds")
#monaco_fine_cache <- file.path(output_dir, "monaco_fine.rds")
#hpca_main_cache <- file.path(output_dir, "hpca_main.rds")
#hpca_fine_cache <- file.path(output_dir, "hpca_fine.rds")
#dice_main_cache <- file.path(output_dir, "dice_main.rds")
#dice_fine_cache <- file.path(output_dir, "dice_fine.rds")
#
## Run SingleR and cache results
#monaco.main <- cache_file(monaco_main_cache, SingleR, test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.main)
#monaco.fine <- cache_file(monaco_fine_cache, SingleR, test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
#hpca.main <- cache_file(hpca_main_cache, SingleR, test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.main)
#hpca.fine <- cache_file(hpca_fine_cache, SingleR, test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.fine)
#dice.main <- cache_file(dice_main_cache, SingleR, test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.main)
#dice.fine <- cache_file(dice_fine_cache, SingleR, test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.fine)
#
## Store the results in the Seurat object
#seurat_obj@meta.data$monaco.main <- monaco.main$pruned.labels
#seurat_obj@meta.data$monaco.fine <- monaco.fine$pruned.labels
#seurat_obj@meta.data$hpca.main <- hpca.main$pruned.labels
#seurat_obj@meta.data$hpca.fine <- hpca.fine$pruned.labels
#seurat_obj@meta.data$dice.main <- dice.main$pruned.labels
#seurat_obj@meta.data$dice.fine <- dice.fine$pruned.labels
#
## Generate cell type annotation plots
#message("Generating cell type annotation plots...")
#seurat_obj_monaco <- SetIdent(seurat_obj, value = "monaco.fine")
#plot12 <- DimPlot(object = seurat_obj_monaco, label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
#ggsave(filename = file.path(output_dir, paste0(output_prefix, '_monaco_dimplot.png')), plot = plot12, width = 11, height = 8.5, dpi = 300)
#
#seurat_obj_hpca <- SetIdent(seurat_obj, value = "hpca.fine")
#plot13 <- DimPlot(object = seurat_obj_hpca, label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
#ggsave(filename = file.path(output_dir, paste0(output_prefix, '_hpca_dimplot.png')), plot = plot13, width = 11, height = 8.5, dpi = 300)
#
#seurat_obj_dice <- SetIdent(seurat_obj, value = "dice.fine")
#plot14 <- DimPlot(object = seurat_obj_dice, label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
#ggsave(filename = file.path(output_dir, paste0(output_prefix, '_dice_dimplot.png')), plot = plot14, width = 11, height = 8.5, dpi = 300)
#
message("Saving final Seurat object...")
saveRDS(seurat_obj, file = file.path(output_dir, paste0(output_prefix, '_seurat_results_final.rds')))
#
#message("Script completed successfully.")
#  #message("Cell type annotation...")
#  #monaco.ref <- celldex::MonacoImmuneData()
#  #hpca.ref <- celldex::HumanPrimaryCellAtlasData()
#  #dice.ref <- celldex::DatabaseImmuneCellExpressionData()
#  #
#  #sce <- as.SingleCellExperiment(DietSeurat(seurat_obj))
#  #
#  #monaco.main <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.main)
#  #monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
#  #hpca.main <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.main)
#  #hpca.fine <- SingleR(test = sce, assay.type.test = 1, ref = hpca.ref, labels = hpca.ref$label.fine)
#  #dice.main <- SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.main)
#  #dice.fine <- SingleR(test = sce, assay.type.test = 1, ref = dice.ref, labels = dice.ref$label.fine)
#  #
#  #seurat_obj@meta.data$monaco.main <- monaco.main$pruned.labels
#  #seurat_obj@meta.data$monaco.fine <- monaco.fine$pruned.labels
#  #seurat_obj@meta.data$hpca.main <- hpca.main$pruned.labels
#  #seurat_obj@meta.data$hpca.fine <- hpca.fine$pruned.labels
#  #seurat_obj@meta.data$dice.main <- dice.main$pruned.labels
#  #seurat_obj@meta.data$dice.fine <- dice.fine$pruned.labels
#  #
#  ##original_seurat_obj <- seurat_obj
#  #
#  #message("Generating cell type annotation plots...")
#  #seurat_obj_monaco <- SetIdent(seurat_obj, value = "monaco.fine")
#  #plot12 <- DimPlot(object = seurat_obj_monaco, label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
#  #ggsave(filename = file.path(output_dir, paste0(output_prefix, '_monaco_dimplot.png')), plot = plot12, width = 11, height = 8.5, dpi = 300)
#  #
#  #seurat_obj_hpca <- SetIdent(seurat_obj, value = "hpca.fine")
#  #plot13 <- DimPlot(object = seurat_obj_hpca, label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
#  #ggsave(filename = file.path(output_dir, paste0(output_prefix, '_hpca_dimplot.png')), plot = plot13, width = 11, height = 8.5, dpi = 300)
#  #
#  #seurat_obj_dice <- SetIdent(seurat_obj, value = "dice.fine")
#  #plot14 <- DimPlot(object = seurat_obj_dice, label = TRUE, repel = TRUE, label.size = 3) + NoLegend()
#  #ggsave(filename = file.path(output_dir, paste0(output_prefix, '_dice_dimplot.png')), plot = plot14, width = 11, height = 8.5, dpi = 300)
#  #
#  #message("Saving final Seurat object...")
#  #saveRDS(seurat_obj, file = file.path(output_dir, paste0(output_prefix, '_seurat_results_final.rds')))
##
#  #
#  #message("Script completed successfully.")

}, error = handle_error)