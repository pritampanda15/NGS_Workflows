DATASETS = ["AEL_GSE142213",
"BLCA_GSE149652",
"BRCA_GSE110686",
"ESCA_GSE154763",
"Glioma_GSE102130",
"KIRC_GSE111360",
"LIHC_GSE166635",
"MCC_GSE117988_aPD1aCTLA4",
"MM_GSE117156",
"NET_GSE140312",
"NHL_GSE147944",
"NSCLC_GSE127471",
"OV_GSE115007",
"PAAD_GSE141017",
"PPB_GSE163678",
"PRAD_GSE150692",
"STAD_GSE167297",
"UVM_GSE169609"]

rule all:
    input:
        expand("../results/{dataset}/vlnplot.png", dataset=DATASETS),
        expand("../results/{dataset}/seurat_markers.csv", dataset=DATASETS)

rule generate_vlnplot:
    input:
        h5_file = "../data/{dataset}_expression.h5"
    output:
        plot = "../results/{dataset}/vlnplot.png",
        rds = "../results/{dataset}/seurat_object.rds"
    log:
        "logs/generate_vlnplot.{dataset}.log"
        "logs/seurat_obj.{dataset}.log"
    conda:
        "seurat_analysis"
    shell:
        "Rscript ../scripts/generate_vlnplot.R {input.h5_file} {output.plot} {output.rds} {wildcards.dataset} --log {log}"

rule generate_markers:
    input:
        rds = "../results/{dataset}/seurat_object.rds"
    output:
        markers_csv = "../results/{dataset}/seurat_markers.csv"
    log:
        "logs/generate_markers.{dataset}log"
    conda:
        "seurat_analysis"
    shell:
        "Rscript ../scripts/generate_markers.R  {input.rds} {output.markers_csv} --log {log}"

