source("global.R")
library("shinyWidgets")
library("shinyFiles")


ui <- dashboardPage(
  dashboardHeader(title = 'DR-scRNAseq'),
  dashboardSidebar(
    sidebarMenu(id="tab",
                useShinyjs(),
                menuItem("Homepage", tabName = "home", icon = icon("house")),
                menuItem("Upload and Convert", tabName = "upload", icon = icon("upload")),
                menuItem("Analyze and Process", tabName = "analyze", icon = icon("edit")),
                menuItem("Plot and Explore", tabName = "plot", icon = icon("chart-line")),
                conditionalPanel(condition = "input.tab == 'plot'",
                div(
                    fileInput("fileplot", "Upload File", multiple = FALSE, accept = c(".rds")),
                    actionButton("resetplot", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                    actionButton("runplot", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%")
                    ),
                div(
                  id = "dropdownDiv",
                  selectInput(
                    inputId = "dimplot",
                    label = "Dimensionality Reduction Method:",
                    choices = c("PCA", "UMAP", "t-SNE")
                             )
                  )          ),
                
                
                conditionalPanel(condition = "input.tab == 'analyze'",
                div(
                    fileInput("fileanalyze", "Upload File", multiple = FALSE, accept = c(".rds")),
                    actionButton("resetanalyze", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                    actionButton("runanalyze", "Process", icon = icon("gears"), style = "color: #fff; background-color: #28a745; width: 87.25%")
                    ),          
                div(
                  textInput("ncounts", "nCounts threshold:", value = "200"),
                  div(
                    id = "numericInputDiv",
                    uiOutput("numeric_input"))            
                   ),
                div(
                  textInput("nfeatures", "nFeatures threshold:", value = "2500"),
                  div(
                    id = "numericInputDiv",
                    uiOutput("numeric_input"))            
                   ),
                div(
                  textInput("mtcounts", "mt.counts threshold:", value = "5"),
                  div(
                    id = "numericInputDiv",
                    uiOutput("numeric_input"))            
                   ),
                div(
                  id = "dropdownDiv",
                  selectInput(
                    inputId = "norm",
                    label = "Normalization Method:",
                    choices = c("LogNormalize", "CLR", "RC")
                  ),
                  #downloadButton("download_drobj", "Download Seurat Object")
                  div(
                    style = "text-align: center; margin-top: 20px;",
                    downloadButton("download_drobj", "Download Seurat Object", style = "color: #ffffff; background-color: #007bff; border-color: #007bff;")
                  )
                  )          ),
                
                
                
                
                conditionalPanel(condition = "input.tab == 'upload'",
                div(
                    fileInput("fileupload", "Upload h5 File", multiple = FALSE, accept = c(".h5", ".h5ad")),
                    actionButton("resetupload", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                    actionButton("convertupload", "Convert", icon = icon("arrows-rotate"), style = "color: #fff; background-color: #28a745; width: 87.25%"),
                    #downloadButton("downloaddata", "Download Seurat Object")
                    ),
                div(
                  style = "text-align: center; margin-top: 20px;",
                  downloadButton("downloaddata", "Download Seurat Object", style = "color: #ffffff; background-color: #007bff; border-color: #007bff;")
                ),
                div(
                  style = "text-align: left; margin-top: 50px; margin-left: 10px;",
                  h5(strong("Upload a folder"))),
                div(
                  shinyDirButton("dir1", "Select data directory", "Upload"),
                  div(
                    style = "margin: 10px auto; width: 85%; background-color: #f8f9fa; border: 2px solid #dee2e6; padding: 5px; border-radius: 2px;",
                    verbatimTextOutput("dirpath", placeholder = TRUE)
                  )
                ),
                div(
                  style = "text-align: left; margin-top: 40px;",
                  actionButton("resetdir", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                  actionButton("convertdir", "Convert", icon = icon("arrows-rotate"), style = "color: #fff; background-color: #28a745; width: 87.25%"),
                  #downloadButton("downloaddir", "Download Seurat Object")
                  ),
                div(
                  style = "text-align: center; margin-top: 20px;",
                  downloadButton("downloaddir", "Download Seurat Object", style = "color: #ffffff; background-color: #007bff; border-color: #007bff;")
                )
                          )
                )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "home",
              tabsetPanel(id = "main_tabs",
                          tabPanel("Instructions", 
                                   includeMarkdown("instructions.Rmd"),
                                   conditionalPanel(condition = "input.tab == 'upload'",
                                      tags$a(
                                        href = "https://drive.google.com/uc?export=download&id=1l1iFIx05iu_5pHw1ciWrhogPdpk0L9K5",  # The URL you want to link to
                                        class = "btn btn-primary",        # Bootstrap button classes
                                        "Download sample h5 file"                    # Button text
                                      ),
                                      tags$a(
                                        href = "https://drive.google.com/uc?export=download&id=17IRibxal4LdWzfJKbJpCNk4RqHyLx9x9",  # The URL you want to link to
                                        class = "btn btn-primary",        # Bootstrap button classes
                                        "Download sample folder"                    # Button text
                                      )
                                   )
                                  
                                   
                        ))
              )
            )
              )
)

server <- function(input, output, session){
  options(shiny.maxRequestSize = 1024*1024^2) #1GB cap
  
  shinyjs::disable("runplot")
  
#if a file is uploaded, run-plot button will be available
  observe({
    if (is.null(input$fileplot) != TRUE){
      shinyjs::enable("runplot")
    } else {
      shinyjs::disable("runplot")
    }
  })
  
#if a file is uploaded, run-analyze button will be available
  observe({
    if (is.null(input$fileanalyze) != TRUE){
      shinyjs::enable("runanalyze")
    } else {
      shinyjs::disable("runanalyze")
    }
  })
  
#if a file is uploaded, upload button will be available
  observe({
    if (is.null(input$fileupload) != TRUE){
      shinyjs::enable("convertupload")
    } else {
      shinyjs::disable("convertupload")
    }
  })
  
#if reset for run is clicked, run-plot button will not be available
  observeEvent(input$resetplot, {
    shinyjs::reset("fileplot")
    shinyjs::disable("runplot")
  })
  
#if reset for run is clicked, run-analyze button will not be available
  observeEvent(input$resetanalyze, {
    shinyjs::reset("fileanalyze")
    shinyjs::disable("runanalyze")
  })
  
#if reset for upload is clicked, upload button will not be available
  observeEvent(input$resetupload, {
    shinyjs::reset("fileupload")
    shinyjs::disable("convertupload")
  })

#for analyzing
  observeEvent(input$runanalyze,{
    shinyjs::disable("runanalyze")
    show_modal_spinner(text = "Processing the data...")
    obj <- load_seurat_obj(input$fileanalyze$datapath)
    
    if (is.vector(obj)){
      showModal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
             paste(unlist(obj), collapse = "<br><br>"))
      ))
      shinyjs::enable("runanalyze")
  } else {
    
    #data processing
    
    ncounts <- as.numeric(input$ncounts)
    nfeatures <- as.numeric(input$nfeatures)
    mtcounts <- as.numeric(input$mtcounts)
    dimred_obj <- seurat_processing(obj, ncounts, nfeatures, mtcounts, input$norm)
    
    #download
    output$download_drobj <- downloadHandler(
      filename = function(){
        paste0('SeuratObj', '.rds')
          },
      content = function(file){
        saveRDS(dimred_obj, file = file)
          }
        )
    
    remove_modal_spinner()
    shinyjs::enable("runanalyze")
  }})
  
#for plotting 
  observeEvent(input$runplot, {
    shinyjs::disable("runplot")
    
    show_modal_spinner(text = "Preparing plots...")
    obj <- load_seurat_obj(input$fileplot$datapath)
    if (is.vector(obj)){
      showModal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
             paste(unlist(obj), collapse = "<br><br>"))
      ))
      shinyjs::enable("runplot")
      
    } else {
      
      #for UMAP
      if (input$dimplot == "UMAP"){
      output$umap <- renderPlotly({
        if (!is.null(input$metadata_col)) {
          #create_metadata_umap_hover(obj, input$metadata_col)
          create_metadata_umap_hover(obj, input$metadata_col)$plotly
        }
      })
      
      output$umap3d <- renderPlotly({
        if (!is.null(input$metadata_col)) {
          create_metadata_umap_3d(obj, input$metadata_col)$plotly
        }
      })
      
      output$featurePlotUMAP <- renderPlot({
        if (!is.null(input$gene)) {
          create_feature_plot_umap(obj, input$gene)
        }
      })
      
      insertTab(
        inputId = "main_tabs",
        tabPanel(
          "2D UMAP",
          fluidRow(
            column(
              width = 8,
              plotlyOutput(outputId = 'umap'),
              downloadButton("download_umaps", "Download UMAP")
            ),
            column(
              width = 4,
              selectizeInput("metadata_col", 
                             "Metadata Column", 
                             colnames(obj@meta.data),
                             selected = "cell_type"
              )
            )
          ),
          style = "height: 90%; width: 95%; padding-top: 5%;"
        ),
        select = TRUE
      )
      
      insertTab(
        inputId = "main_tabs",
        tabPanel(
          "3D UMAP",
          fluidRow(
            column(
              width = 8,
              plotlyOutput(outputId = 'umap3d')
            ),
            column(
              width = 4,
              selectizeInput("metadata_col", 
                             "Metadata Column", 
                             colnames(obj@meta.data),
                             selected = "cell_type"
              )
            )
          ),
          style = "height: 90%; width: 95%; padding-top: 5%;"
        ),
        select = TRUE
      )
      
      insertTab(
        inputId = "main_tabs",
        tabPanel(
          "Gene Expression",
          fluidRow(
            column(
              width = 8,
              plotOutput(outputId = 'featurePlotUMAP'),
              downloadButton("downloadFeaturePlot", "Download Feature Plot")
            ),
            column(
              width = 4,
              selectizeInput("gene", 
                             "Genes", 
                             rownames(obj),
                             selected = "ISG15"
              )
            )
          ),
          style = "height: 90%; width: 95%; padding-top: 5%;"
        )
      )
      
      remove_modal_spinner()
      shinyjs::enable("runplot")
      
      ###############
      output$download_umaps <- downloadHandler(
        filename = function(){
          "myplot.png"
        },
        content = function(file){
          #p <- create_metadata_umap_hover(obj, input$metadata_col)$x$attrs[[1]]$plot
          #ggsave(filename = file, plot = p, device = "png")
          plot_data <- create_metadata_umap_hover(obj, input$metadata_col)
          ggsave(filename = file, plot = plot_data$ggplot, device = "png")
        }
      )
      ##############
      
      }
      
      #for PCA
      else if (input$dimplot == "PCA"){
        output$pca <- renderPlotly({
          if (!is.null(input$metadata_col)) {
            create_metadata_pca_hover(obj, input$metadata_col)$plotly
          }
        })
        
        output$pca3d <- renderPlotly({
          if (!is.null(input$metadata_col)) {
            create_metadata_pca_3d(obj, input$metadata_col)
          }
        })
        
        output$featurePlotPCA <- renderPlot({
          if (!is.null(input$gene)) {
            create_feature_plot_pca(obj, input$gene)
          }
        })
        
        insertTab(
          inputId = "main_tabs",
          tabPanel(
            "2D PCA",
            fluidRow(
              column(
                width = 8,
                plotlyOutput(outputId = 'pca'),
                downloadButton("download_pca", "Download PCA")
              ),
              column(
                width = 4,
                selectizeInput("metadata_col", 
                               "Metadata Column", 
                               colnames(obj@meta.data),
                               selected = "cell_type"
                )
              )
            ),
            style = "height: 90%; width: 95%; padding-top: 5%;"
          ),
          select = TRUE
        )
        
        insertTab(
          inputId = "main_tabs",
          tabPanel(
            "3D PCA",
            fluidRow(
              column(
                width = 8,
                plotlyOutput(outputId = 'pca3d')
              ),
              column(
                width = 4,
                selectizeInput("metadata_col", 
                               "Metadata Column", 
                               colnames(obj@meta.data),
                               selected = "cell_type"
                )
              )
            ),
            style = "height: 90%; width: 95%; padding-top: 5%;"
          ),
          select = TRUE
        )
        
        insertTab(
          inputId = "main_tabs",
          tabPanel(
            "Gene Expression",
            fluidRow(
              column(
                width = 8,
                plotOutput(outputId = 'featurePlotPCA'),
                downloadButton("downloadFeaturePlot", "Download Feature Plot")
              ),
              column(
                width = 4,
                selectizeInput("gene", 
                               "Genes", 
                               rownames(obj),
                               selected = "ISG15"
                )
              )
            ),
            style = "height: 90%; width: 95%; padding-top: 5%;"
          )
        )
        
        remove_modal_spinner()
        shinyjs::enable("runplot")
        
        ###############
        output$download_pca <- downloadHandler(
          filename = function(){
            "my_pca_plot.png"
          },
          content = function(file){
            #p <- create_metadata_umap_hover(obj, input$metadata_col)$x$attrs[[1]]$plot
            #ggsave(filename = file, plot = p, device = "png")
            plot_data <- create_metadata_pca_hover(obj, input$metadata_col)
            ggsave(filename = file, plot = plot_data$ggplot, device = "png")
          }
        )
        ##############
        
      }
      
      #for t-SNE
      else if (input$dimplot == "t-SNE"){
        output$tsne <- renderPlotly({
          if (!is.null(input$metadata_col)) {
            create_metadata_tsne_hover(obj, input$metadata_col)$plotly
          }
        })
        
        output$tsne3d <- renderPlotly({
          if (!is.null(input$metadata_col)) {
            create_metadata_tsne_3d(obj, input$metadata_col)
          }
        })
        
        output$featurePlotTSNE <- renderPlot({
          if (!is.null(input$gene)) {
            create_feature_plot_tsne(obj, input$gene)
          }
        })
        
        insertTab(
          inputId = "main_tabs",
          tabPanel(
            "2D t-SNE",
            fluidRow(
              column(
                width = 8,
                plotlyOutput(outputId = 'tsne'),
                downloadButton("download_tsne", "Download t-SNE")
              ),
              column(
                width = 4,
                selectizeInput("metadata_col", 
                               "Metadata Column", 
                               colnames(obj@meta.data),
                               selected = "cell_type"
                )
              )
            ),
            style = "height: 90%; width: 95%; padding-top: 5%;"
          ),
          select = TRUE
        )
        
        insertTab(
          inputId = "main_tabs",
          tabPanel(
            "3D t-SNE",
            fluidRow(
              column(
                width = 8,
                plotlyOutput(outputId = 'tsne3d')
              ),
              column(
                width = 4,
                selectizeInput("metadata_col", 
                               "Metadata Column", 
                               colnames(obj@meta.data),
                               selected = "cell_type"
                )
              )
            ),
            style = "height: 90%; width: 95%; padding-top: 5%;"
          ),
          select = TRUE
        )
        
        insertTab(
          inputId = "main_tabs",
          tabPanel(
            "Gene Expression",
            fluidRow(
              column(
                width = 8,
                plotOutput(outputId = 'featurePlotTSNE'),
                downloadButton("downloadFeaturePlot", "Download Feature Plot")
              ),
              column(
                width = 4,
                selectizeInput("gene", 
                               "Genes", 
                               rownames(obj),
                               selected = "ISG15"
                )
              )
            ),
            style = "height: 90%; width: 95%; padding-top: 5%;"
          )
        )
        
        remove_modal_spinner()
        shinyjs::enable("runplot")
        
        ###############
        output$download_tsne <- downloadHandler(
          filename = function(){
            "my_tsne_plot.png"
          },
          content = function(file){
            #p <- create_metadata_umap_hover(obj, input$metadata_col)$x$attrs[[1]]$plot
            #ggsave(filename = file, plot = p, device = "png")
            plot_data <- create_metadata_tsne_hover(obj, input$metadata_col)
            ggsave(filename = file, plot = plot_data$ggplot, device = "png")
          }
        )
        ##############
        
      }
      
      
      }
  }
  
  )
  
  
#for uploading h5
  observeEvent(input$convertupload, {
    shinyjs::disable("convertupload")
    show_modal_spinner(text = "Converting...")
    converted_obj <- load_h5(input$fileupload$datapath)
      
      
    #download
    output$downloaddata <- downloadHandler(
      filename = function(){
        paste0('SeuratObjh5', '.rds')
      },
      content = function(file){
        saveRDS(converted_obj, file = file)
      }
      )
      
      remove_modal_spinner()
      shinyjs::enable("convertupload")
      
  })

   
#clear all sidebar inputs when 'Reset' button is clicked for run
  observeEvent(input$resetplot, {
    shinyjs::reset("file")
    removeTab("main_tabs", "UMAP")
    removeTab("main_tabs", "PCA")
    removeTab("main_tabs", "t-SNE")
    removeTab("main_tabs", "Gene Expression")
    shinyjs::disable("run")
  })
  

  roots <- c(home = normalizePath("~"))
  shinyDirChoose(
    input,
    'dir1',
    roots = roots,
    filetypes = c('', 'mtx', "tsv", "csv")
  )
  
  
  output$dirpath <- renderPrint({
    req(input$dir1)
    dirpath <- parseDirPath(roots, input$dir1)
    dirpath
  })
  
  
  #for uploading gz
  observeEvent(input$convertdir, {
    shinyjs::disable("convertdir")
    show_modal_spinner(text = "Converting...")
    selected_dir <- parseDirPath(roots, input$dir1)
    obj <- load_gz(selected_dir)
    
    #download
    output$downloaddir <- downloadHandler(
      filename = function(){
        paste0('SeuratObjh5', '.rds')
      },
      content = function(file){
        saveRDS(obj, file = file)
      }
    )
    
    remove_modal_spinner()
    shinyjs::enable("convertupload")
  })
  
  
  
  
}

shinyApp(ui, server)
