# r-seurat
Seurat R pipeline to process 10x .h5 files

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install CRAN packages
install.packages(c("hdf5r", "pacman", "ggsci", "RColorBrewer", "ggrepel", 
                   "pheatmap", "viridis", "argparse", "cowplot", "dplyr", 
                   "ggplot2", "patchwork", "matrix.utils", "reticulate"))

# Install package from a specific repository
install.packages("RPresto", repos = "https://conda-forge.org/r")

# Install Bioconductor packages
BiocManager::install(c("celldex", "SingleR", "SingleCellExperiment", 
                       "glmGamPoi", "scran"))
```
