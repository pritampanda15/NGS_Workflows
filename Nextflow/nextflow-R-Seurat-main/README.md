# Post Processing of 10xscRNAseq *.h5 files in Seurat
This pipeline comprised of reading 10x h5 files to annotation using SingleR including Cell cycle sorting.

## Make sure you have .h5 files inside a folder called data/
```
nextflow -log logs/seurat.nf run nf_run.nf -c nextflow.config --with-report scrnaseq_seurat.html
```
## You can specify continer image (docker) in the nextflow.config file
```
apptainer pull  --name pritampkp15-scrna-seurat-analysis.img  docker://pritampkp15/scrna-seurat-analysis
```
