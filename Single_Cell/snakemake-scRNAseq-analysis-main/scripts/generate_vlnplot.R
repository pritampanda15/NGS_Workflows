pacman::p_load(
    Seurat, ggplot2, patchwork, tidyverse, hdf5r, ggsci, 
    celldex, RColorBrewer, SingleCellExperiment, glmGamPoi, 
    reticulate, cowplot, viridis, pheatmap, scran, SingleR, 
    BiocParallel, DoubletFinder
)


tryCatch({
    args <- commandArgs(trailingOnly = TRUE)
    h5_file <- args[1]
    output_plot <- args[2]  # This will now include the dataset-specific path
    output_rds <- args[3]   # Adjusted to include the dataset name in the path

    # Ensure the output directory exists
    output_dir <- dirname(output_plot)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    # Load data
    adj_matrix <- Seurat::Read10X_h5(filename = h5_file, use.names = TRUE)
    seurat_obj <- CreateSeuratObject(counts = adj_matrix, project = 'TISCH2', min.cells = 3, min.features = 200)

    # Additional processing and plot generation
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
    seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
    
    plot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.1) + 
            theme(plot.title = element_text(size=10))

    # Save plot
    ggsave(filename = output_plot, plot = plot, width = 11, height = 8.5, dpi = 300)

    # Save the Seurat object as an RDS file
    saveRDS(seurat_obj, file = output_rds)

}, error = function(e) {
    cat("Error in processing Seurat object: ", e$message, "\n")
    quit(status = 1)
})
