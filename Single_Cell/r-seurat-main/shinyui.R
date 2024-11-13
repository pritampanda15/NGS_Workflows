library(shiny)
library(Seurat)
library(ggplot2)
library(patchwork)

# Set max upload size to 2GB
options(shiny.maxRequestSize = 2000 * 1024^2)

# Define UI
ui <- fluidPage(
  titlePanel("Seurat Processing App"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload h5 file"),
      actionButton("process_sctransform", "Perform SCTransform"),
      actionButton("process_pca", "Run PCA"),
      actionButton("process_umap", "Run UMAP"),
      actionButton("process_clustering", "Find Clusters"),
      downloadButton("downloadRDS", "Download Seurat Object")
    ),
    mainPanel(
      plotOutput("vlnplot"),
      plotOutput("dimplot")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  seurat_obj <- reactiveVal(NULL)
  
  observeEvent(input$file, {
    req(input$file)
    h5_file <- input$file$datapath
    
    # Read uploaded h5 file and create Seurat object
    adj_matrix <- Read10X_h5(h5_file, use.names = TRUE)
    seurat_obj(CreateSeuratObject(counts = adj_matrix, project = 'TISCH2', min.cells = 3, min.features = 200))
    
    # Calculate mitochondrial and ribosomal percentages
    obj <- seurat_obj()
    obj[['percent.mt']] <- PercentageFeatureSet(obj, pattern = '^MT')
    obj[['percent.rb']] <- PercentageFeatureSet(obj, pattern = '^RP[SL]')
    seurat_obj(obj)
    
    output$vlnplot <- renderPlot({
      VlnPlot(object = seurat_obj(), features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'percent.rb'), ncol = 4)
    })
  })
  
  observeEvent(input$process_sctransform, {
    req(seurat_obj())
    
    # Perform SCTransform normalization
    obj <- SCTransform(object = seurat_obj(), assay = "RNA", vars.to.regress = "percent.mt", verbose = FALSE)
    seurat_obj(obj)
  })
  
  observeEvent(input$process_pca, {
    req(seurat_obj())
    
    # Run PCA
    obj <- RunPCA(seurat_obj(), features = VariableFeatures(object = seurat_obj()))
    seurat_obj(obj)
    
    output$dimplot <- renderPlot({
      DimPlot(object = seurat_obj(), reduction = "pca", label = TRUE)
    })
  })
  
  observeEvent(input$process_umap, {
    req(seurat_obj())
    
    # Run UMAP
    obj <- RunUMAP(seurat_obj(), dims = 1:30)
    seurat_obj(obj)
    
    output$dimplot <- renderPlot({
      DimPlot(object = seurat_obj(), reduction = "umap", label = TRUE)
    })
  })
  
  observeEvent(input$process_clustering, {
    req(seurat_obj())
    
    # Find neighbors and clusters
    obj <- FindNeighbors(seurat_obj(), dims = 1:30)
    obj <- FindClusters(obj, verbose = FALSE)
    seurat_obj(obj)
    
    output$dimplot <- renderPlot({
      DimPlot(object = seurat_obj(), reduction = "umap", label = TRUE)
    })
  })
  
  output$downloadRDS <- downloadHandler(
    filename = function() {
      paste("seurat_object_", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(seurat_obj(), file)
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
