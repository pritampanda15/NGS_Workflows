import os
import shutil

DATASETS = {
"blood": ["AEL_GSE142213"],
"bladder":["BLCA_GSE149652"],
"breast":["BRCA_GSE110686"],
"esophagus":["ESCA_GSE154763"],
"brain":["Glioma_GSE102130"],
"kidney":["KIRC_GSE111360"],
"liver":["LIHC_GSE166635"],
"skin":["MCC_GSE117988_aPD1aCTLA4"],
"bone":["MM_GSE117156"],
"nervous":["NET_GSE140312"],
"lymph":["NHL_GSE147944"],
"lung":["NSCLC_GSE127471"],
"pelvic":["OV_GSE115007"],
"pancreas":["PAAD_GSE141017"],
"soft_tissue":["PPB_GSE163678"],
"prostate":["PRAD_GSE150692"],
"stomach":["STAD_GSE167297"],
"eye":["UVM_GSE169609"]
}

output_base_path = "../results"  # Adjust this path as necessary

for category, datasets in DATASETS.items():
    for dataset in datasets:
        # Define current and target paths
        current_path = os.path.join(output_base_path, dataset)
        target_path = os.path.join(output_base_path, category, dataset)
        
        # Check if the current dataset output exists
        if os.path.exists(current_path):
            # Move the dataset folder to the correct category folder
            shutil.move(current_path, target_path)
            print(f"Moved {current_path} to {target_path}")
