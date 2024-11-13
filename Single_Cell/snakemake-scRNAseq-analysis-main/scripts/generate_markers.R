pacman::p_load(
    Seurat, ggplot2, patchwork, tidyverse, hdf5r, ggsci, 
    celldex, RColorBrewer, SingleCellExperiment, glmGamPoi, 
    reticulate, cowplot, viridis, pheatmap, scran, SingleR, 
    BiocParallel, DoubletFinder
)



# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]  # Seurat object
output_csv <- args[2]  # Output path for markers
seurat_obj <- readRDS(args[1])

# Ensure the output directory exists
output_dir <- dirname(output_csv)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(seurat_obj), 10)

all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15) #Depending on Dimensionality,clustering will bet set accordingly
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)


seurat_obj <- RunPCA(seurat_obj, verbose = F)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:15)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

all.markers <- FindAllMarkers(seurat_obj, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(all.markers, file = output_csv, row.names = FALSE)
