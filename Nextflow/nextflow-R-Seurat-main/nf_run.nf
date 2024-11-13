#!/usr/bin/env nextflow
// Nextflow DSL2 pipeline for processing scRNA data using Seurat
//shebang line to define the environment

//enable DSL2 using NXF_VER=23.10.1
nextflow.enable.dsl=2


/**
 * Define the parameters for the script using the params object 
    * In Nextflow, `params` is a special object used to handle parameters that can be passed to the script at runtime. 
    * These parameters can be used to customize the behavior of the script without having to modify the script itself.
    * nextflow run myscript.nf --data_folder /path/to/data --output_folder /path/to/results (cn also run the script like this)
    * In this command, the --data_folder and --output_folder options are used to provide new values for the params.data_folder and params.output_folder parameters, 
    * which will override the default values set in the script.
 */

// Define parameters
params.data_folder = './data/'
params.output_folder = './results'

/**
    * "Processing data from folder: ${params.data_folder}"
    * "Saving results to folder: ${params.output_folder}"
    * The log message is defined using a multi-line string, which is enclosed in triple quotes (""").
    * This allows the log message to span multiple lines, which makes it easier to read and understand.
*/

log.info """\
    S E U R A T - R  - P I P E L I N E
    ===================================
    data: ${params.data_folder}
    outdir: ${params.output_folder}
    """
    .stripIndent(true) // Strip leading whitespace from the log message

/**
 * Define a process named SeuratAnalysis to process the scRNA data using Seurat starts with process
    * In Nextflow, `process` is a special keyword used to define a process that will be executed as part of the workflow.
    * Each process is defined using a block of code that contains the process name, input and output declarations, and a script to be executed.
    * The input and output declarations are used to define the input and output channels for the process, which are used to connect the process to other processes in the workflow.
    * The script block contains the code that will be executed when the process is run.
    * The process name is used to identify the process in the workflow, and can be used to tag the process for better traceability.
 */
    process SeuratAnalysis { 
    tag "$h5_id"    // Tag the process with the sample ID (in this case h5_id) for better traceability otherwise sample_id can be used
    publishDir "${params.output_folder}/${h5_file.baseName}", mode: 'copy'

    input:
    path h5_file // Input channel for the process, which accepts a list of files matching the pattern /path/to/data/*.h5 you can also name it as path input_file

    output:
    file("${h5_file.baseName}_seurat_results.rds")
    file("${h5_file.baseName}_vlnplot.png")
    file("${h5_file.baseName}_featurescatter1.png")
    file("${h5_file.baseName}_featurescatter2.png")
    file("${h5_file.baseName}_featurescatter3.png")
    file("${h5_file.baseName}_featurescatter4.png")
    file("${h5_file.baseName}_PCA_dimplot.png")
    file("${h5_file.baseName}_PCA_elbowplot.png")
    file("${h5_file.baseName}_DimPlot_UMAP.png")
    file("${h5_file.baseName}_DimPlot_cellcycle+percentmt.png")
    file("${h5_file.baseName}_feature_cellcycle.png")
    file("${h5_file.baseName}_heatmap.png")
    file("${h5_file.baseName}_output.csv")
    file("${h5_file.baseName}_seurat_results_final.rds")

    script:
    h5_id = h5_file.baseName
    """
    mkdir -p ${h5_id}
    seurat.R ${h5_file}
    """
    }
    workflow {
    Channel.fromPath("${params.data_folder}/*.h5") | SeuratAnalysis
}