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
                    fileInput("file", "Upload File", multiple = FALSE, accept = c(".rds")),
                    actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                    actionButton("run-analyze", "Process", icon = icon("gears"), style = "color: #fff; background-color: #28a745; width: 87.25%")
                    ),          
                div(
                  textInput("label_input", "nCounts threshold:", value = " "),
                  div(
                    id = "numericInputDiv",
                    uiOutput("numeric_input"))            
                   ),
                div(
                  textInput("label_input", "nFeatures threshold:", value = " "),
                  div(
                    id = "numericInputDiv",
                    uiOutput("numeric_input"))            
                   ),
                div(
                  textInput("label_input", "mt.counts threshold:", value = " "),
                  div(
                    id = "numericInputDiv",
                    uiOutput("numeric_input"))            
                   ),
                div(
                  id = "dropdownDiv",
                  selectInput(
                    inputId = "myDropdown",
                    label = "Normalization Method:",
                    choices = c("LogNormalize", "CLR", "RC")
                  )
                  )          ),
                
                
                conditionalPanel(condition = "input.tab == 'upload'",
                div(
                    fileInput("fileupload", "Upload h5 File", multiple = FALSE, accept = c(".h5")),
                    actionButton("resetupload", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                    actionButton("convertupload", "Convert", icon = icon("arrows-rotate"), style = "color: #fff; background-color: #28a745; width: 87.25%")
                    ),
                div(
                  h5(strong("Upload a folder")),
                  shinyDirButton("dir1", "Select data directory", "Upload"),
                  verbatimTextOutput("dir1", placeholder = TRUE),
                  actionButton("resetdir", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                  actionButton("convertdir", "Convert", icon = icon("arrows-rotate"), style = "color: #fff; background-color: #28a745; width: 87.25%"),
                  downloadButton("downloaddata", "Download Seurat Object")
                  )
                          )
                )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "home",
              tabsetPanel(id = "main_tabs",
                          tabPanel("Instructions", 
                                   includeMarkdown("instructions.Rmd")
                        ))
              )
            )
              )
)

server <- function(input, output, session){
  options(shiny.maxRequestSize = 300*1024^2)
  
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
    if (is.null(input$file) != TRUE){
      shinyjs::enable("run-analyze")
    } else {
      shinyjs::disable("run-analyze")
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
  observeEvent(input$reset, {
    shinyjs::reset("file")
    shinyjs::disable("run-analyze")
  })
  
#if reset for upload is clicked, upload button will not be available
  observeEvent(input$resetupload, {
    shinyjs::reset("fileupload")
    shinyjs::disable("convertupload")
  })
  
#input
 # output$numeric_input <- renderUI({
 #    numericInput(
 #      inputId = "numeric_input",
 #      label = input$label_input,
 #      value = 0
 #    )
 #  })

  
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
      output$umap <- renderPlot({
        if (!is.null(input$metadata_col)) {
          create_metadata_umap(obj, input$metadata_col)
        }
      })
      
      output$featurePlotUMAP <- renderPlot({
        if (!is.null(input$gene)) {
          create_feature_plot(obj, input$gene)
        }
      })
      
      insertTab(
        inputId = "main_tabs",
        tabPanel(
          "UMAP",
          fluidRow(
            column(
              width = 8,
              plotOutput(outputId = 'umap'),
              downloadButton("download_umap", "Download UMAP")
            ),
            column(
              width = 4,
              selectizeInput("metadata_col", 
                             "Metadata Column", 
                             colnames(obj@meta.data)
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
                             rownames(obj)
              )
            )
          ),
          style = "height: 90%; width: 95%; padding-top: 5%;"
        )
      )
      
      remove_modal_spinner()
      shinyjs::enable("runplot")
      
      }
      
      #for PCA
      else if (input$dimplot == "PCA"){
        output$pca <- renderPlot({
          if (!is.null(input$metadata_col)) {
            create_metadata_pca(obj, input$metadata_col)
          }
        })
        
        output$featurePlotPCA <- renderPlot({
          if (!is.null(input$gene)) {
            create_feature_plot(obj, input$gene)
          }
        })
        
        insertTab(
          inputId = "main_tabs",
          tabPanel(
            "PCA",
            fluidRow(
              column(
                width = 8,
                plotOutput(outputId = 'pca'),
                downloadButton("download_pca", "Download PCA")
              ),
              column(
                width = 4,
                selectizeInput("metadata_col", 
                               "Metadata Column", 
                               colnames(obj@meta.data)
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
                               rownames(obj)
                )
              )
            ),
            style = "height: 90%; width: 95%; padding-top: 5%;"
          )
        )
        
        remove_modal_spinner()
        shinyjs::enable("runplot")
        
      }
      
      #for t-SNE
      else if (input$dimplot == "t-SNE"){
        output$tsne <- renderPlot({
          if (!is.null(input$metadata_col)) {
            create_metadata_tsne(obj, input$metadata_col)
          }
        })
        
        output$featurePlotTSNE <- renderPlot({
          if (!is.null(input$gene)) {
            create_feature_plot(obj, input$gene)
          }
        })
        
        insertTab(
          inputId = "main_tabs",
          tabPanel(
            "t-SNE",
            fluidRow(
              column(
                width = 8,
                plotOutput(outputId = 'tsne'),
                downloadButton("download_tsne", "Download t-SNE")
              ),
              column(
                width = 4,
                selectizeInput("metadata_col", 
                               "Metadata Column", 
                               colnames(obj@meta.data)
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
                               rownames(obj)
                )
              )
            ),
            style = "height: 90%; width: 95%; padding-top: 5%;"
          )
        )
        
        remove_modal_spinner()
        shinyjs::enable("runplot")
        
      }
      }
  })
  
#for converting to h5
  # observeEvent(input$convertupload, {
  #   shinyjs::disable("upload")
  #   
  #   show_modal_spinner(text = "Converting...")
  #   h5 <- load_h5(input$upload$datapath)
  #   if (is.vector(obj)){
  #     showModal(modalDialog(
  #       title = "Error with file",
  #       HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
  #            paste(unlist(obj), collapse = "<br><br>"))
  #     ))
  #     shinyjs::enable("run")
  #     
  #   } else {
  #     
  #     #h5 to seurat object
  #     output$seuratobj <- h5
  # 
  #     
  #     #download the rds file
  #     output$download_seuratobj <- downloadHandler(
  #       filename = function(){
  #         paste0(input$metadata_col, '_SeuratObj', '.rds')
  #       },
  #       content = function(file){
  #         saveRDS(seuratobj, file = file)
  #       }
  #     )
  #     
  #     
  #     
  #     remove_modal_spinner()
  #     shinyjs::enable("convertupload") #stophere
  #     
  #   }
  # })
  
  observeEvent(input$convertupload, {
    req(input$fileupload)
    shinyjs::disable("convertupload")
    show_modal_spinner(text = "Converting...")
    
    coundata <- Seurat::Read10X_h5(input$fileupload$datapath)
    seurat_obj$data <- Seurat::CreateSeuratObject(counts = count_data, project = "MySeuratObject")
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("seurat_object_", Sys.Date(), ".rds", sep = "")
      },
      content = function(file) {
        saveRDS(seurat_obj$data, file = file)
      }
    )
  }
              )

  
#clear all sidebar inputs when 'Reset' button is clicked for run
  observeEvent(input$resetplot, {
    shinyjs::reset("file")
    removeTab("main_tabs", "UMAP")
    removeTab("main_tabs", "PCA")
    removeTab("main_tabs", "t-SNE")
    removeTab("main_tabs", "Gene Expression")
    shinyjs::disable("run")
  })
  
  
  
  shinyDirChoose(
    input,
    'dir1',
    roots = c(home = '~'),
    filetypes = c('', 'mtx', "tsv", "csv", "bw")
  )
  
  global <- reactiveValues(datapath = getwd())
  
  dir1 <- reactive(input$dir1)
  
  output$dir1 <- renderText({
    global$datapath
  })
  
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir1
               },
               handlerExpr = {
                 if (!"path" %in% names(dir1())) return()
                 home <- normalizePath("~")
                 global$datapath <-
                   file.path(home, paste(unlist(dir1()$path[-1]), collapse = .Platform$file.sep))
               })
  
  observeEvent(input$convertdir, {
    req(dir1())
    shinyjs::disable("convertdir")
    show_modal_spinner(text = "Processing...")
    
    count_data <- Read10X(dir1())
    seurat_obj$data <- CreateSeuratObject(counts = count_data, project = "MySeuratObject")
    remove_modal_spinner()
    shinyjs::enable("convertdir")
  })
  
  output$downloaddata <- downloadHandler(
    filename = function() {
      paste("seurat_object_", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(seurat_obj$data, file = file)
    })
}

shinyApp(ui, server)
