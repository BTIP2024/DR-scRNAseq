source("global.R")


ui <- dashboardPage(
  dashboardHeader(title = 'DR-scRNAseq'),
  dashboardSidebar(
    sidebarMenu(id="tab",
                useShinyjs(),
                menuItem("Homepage", tabName = "home", icon = icon("house")),
                menuItem("Upload and Convert", tabName = "upload", icon = icon("upload")),
                menuItem("Analyze and Process", tabName = "analyze", icon = icon("edit")),
                menuItem("Plot and Explore", tabName = "explore", icon = icon("chart-line")),
                conditionalPanel(condition = "input.tab == 'explore'",
                div(
                    fileInput("file", "Upload File", multiple = FALSE, accept = c(".rds")),
                    actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                    actionButton("run-plot", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%")
                    )            ),
                
                
                conditionalPanel(condition = "input.tab == 'analyze'",
                div(
                    fileInput("file", "Upload File", multiple = FALSE, accept = c(".rds")),
                    actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                    actionButton("run-analyze", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%")
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
                   )         ),
                
                
                conditionalPanel(condition = "input.tab == 'upload'",
                div(
                    fileInput("file", "Upload File", multiple = FALSE, accept = c(".h5")),
                    actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                    actionButton("upload", "Convert", icon = icon("arrows-rotate"), style = "color: #fff; background-color: #28a745; width: 87.25%")
                    )            )
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
  
  shinyjs::disable("run")
  
#if a file is uploaded, run-plot button will be available
  observe({
    if (is.null(input$file) != TRUE){
      shinyjs::enable("run-plot")
    } else {
      shinyjs::disable("run-plot")
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
    if (is.null(input$file) != TRUE){
      shinyjs::enable("upload")
    } else {
      shinyjs::disable("upload")
    }
  })
  
#if reset for run is clicked, run-plot button will not be available
  observeEvent(input$reset, {
    shinyjs::reset("file")
    shinyjs::disable("run-plot")
  })
  
#if reset for run is clicked, run-plot button will not be available
  observeEvent(input$reset, {
    shinyjs::reset("file")
    shinyjs::disable("run-analyze")
  })
  
#if reset for upload is clicked, upload button will not be available
  observeEvent(input$reset, {
    shinyjs::reset("file")
    shinyjs::disable("upload")
  })
  
#input
 # output$numeric_input <- renderUI({
 #    numericInput(
 #      inputId = "numeric_input",
 #      label = input$label_input,
 #      value = 0
 #    )
 #  })

  
#for plotting tSNE
  observeEvent(input$run, {
    shinyjs::disable("run")
    
    show_modal_spinner(text = "Preparing plots...")
    obj <- load_seurat_obj(input$file$datapath)
    if (is.vector(obj)){
      showModal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
             paste(unlist(obj), collapse = "<br><br>"))
      ))
      shinyjs::enable("run")
      
    } else {
      
      output$umap <- renderPlot({
        if (!is.null(input$metadata_col)) {
          create_metadata_umap(obj, input$metadata_col)
        }
      })
      
      output$featurePlot <- renderPlot({
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
              plotOutput(outputId = 'featurePlot'),
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
      shinyjs::enable("run")
      
    }
  })
  
#for converting to h5
  observeEvent(input$upload, {
    shinyjs::disable("upload")
    
    show_modal_spinner(text = "Converting...")
    obj <- load_seurat_obj(input$upload$datapath)
    if (is.vector(obj)){
      showModal(modalDialog(
        title = "Error with file",
        HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
             paste(unlist(obj), collapse = "<br><br>"))
      ))
      shinyjs::enable("run")
      
    } else {
      
      #h5 to seurat object
      output$seuratobj <- load_h5(path)
      
      
      #download the rds file
      output$download_seuratobj <- downloadHandler(
        filename = function(){
          paste0(input$metadata_col, '_SeuratObj', '.rds')
        },
        content = function(file){
          saveRDS(seuratobj, file = file)
        }
      )
      
      
      
      remove_modal_spinner()
      shinyjs::enable("run")
      
    }
  })
  
#clear all sidebar inputs when 'Reset' button is clicked for run
  observeEvent(input$reset, {
    shinyjs::reset("file")
    removeTab("main_tabs", "UMAP")
    removeTab("main_tabs", "Gene Expression")
    shinyjs::disable("run")
  })
  
}

shinyApp(ui, server)
