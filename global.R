library(shiny)
library(shinydashboard)
library(markdown)
library(shinyjs)
library(shinybusy)
library(Seurat)
library(ggplot2)
library(plotly)
library(tools)
library(dplyr)
library(DT)
library(shinydashboardPlus)
library(glue)
library(markdown)
library(ggthemes)

load_seurat_obj <- function(path){
  errors <- c()
  #check file ext
  if(!tolower(tools::file_ext(path)) == 'rds'){
    errors <- c(errors, "Invalid rds file")
    return(errors)
  }
  
  #try to read in file
  tryCatch(
    {
      obj <- readRDS(path)
    },
    error = function(e){
      errors <- c(errors, "Invalid rds file")
      return(errors)
    }
  )
  #validate obj is a seurat obj
  if (!inherits(obj, "Seurat")){
    errors <- c(errors, "File is not a seurat object")
    return(errors)
  }
  return(obj)
}

create_feature_plot_pca <- function(obj, gene) {
  if (gene %in% rownames(obj)) {
    FP <- Seurat::FeaturePlot(obj, features = gene, pt.size = 0.1, combine = FALSE, reduction = "pca")
  } else {
    FP <- ggplot() +
      theme_void() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(FP)
}

create_feature_plot_tsne <- function(obj, gene) {
  if (gene %in% rownames(obj)) {
    FP <- Seurat::FeaturePlot(obj, features = gene, pt.size = 0.001, combine = FALSE, reduction = "tsne")
  } else {
    FP <- ggplot() +
      theme_void() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(FP)
}

create_feature_plot_umap <- function(obj, gene) {
  if (gene %in% rownames(obj)) {
    FP <- Seurat::FeaturePlot(obj, features = gene, pt.size = 0.001, combine = FALSE, reduction = "umap")
  } else {
    FP <- ggplot() +
      theme_void() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  return(FP)
}

seurat_processing <- function(obj, qc1, qc2, qc3, norm){
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nFeature_RNA > qc1 & nFeature_RNA < qc2 & percent.mt < qc3)
  obj <- Seurat::NormalizeData(obj, normalization.method = norm)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(obj)
  obj <- Seurat::ScaleData(obj, features = all.genes)
  obj <- Seurat::RunPCA(obj, features = VariableFeatures(object = obj))
  obj <- Seurat::FindNeighbors(obj, dims = 1:10)
  obj <- Seurat::FindClusters(obj, resolution = 0.5)
  new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
  names(new.cluster.ids) <- levels(obj)
  obj <- RenameIdents(obj, new.cluster.ids)
  obj$cell_type = Idents(obj)
  obj <- Seurat::RunTSNE(obj, dims = 1:10)
  objfinal <- Seurat::RunUMAP(obj, dims = 1:10)
  return(objfinal)
}

load_h5 <- function(path){
  obj <- Seurat::Read10X_h5(path)
  obj <- Seurat::CreateSeuratObject(obj)
  return(obj)
}

load_gz <- function(path){
  obj <- Seurat::Read10X(path)
  obj <- Seurat::CreateSeuratObject(obj)
  return(obj)
}

create_metadata_pca_hover <- function(obj, col){
  if (col %in% colnames(obj@meta.data)) {
    pca <- DimPlot(obj, reduction = "pca", label = TRUE, group.by = col) + xlab("PCA 1") + ylab("PCA 2") + 
      theme(axis.title = element_text(size = 18)) +  
      guides(colour = guide_legend(override.aes = list(size = 10)))
  } else {
    pca <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "col doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  #ggplotly(pca)
  list(ggplot = pca, plotly = ggplotly(pca))
}

create_metadata_umap_hover <- function(obj, col){
  if (col %in% colnames(obj@meta.data)) {
    umap <- DimPlot(obj, reduction = "umap", label = TRUE, group.by = col) + xlab("UMAP 1") + ylab("UMAP 2") + 
      theme(axis.title = element_text(size = 18)) +  
      guides(colour = guide_legend(override.aes = list(size = 10)))
  } else {
    umap <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "col doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  #ggplotly(umap)
  list(ggplot = umap, plotly = ggplotly(umap))
}

create_metadata_tsne_hover <- function(obj, col){
  if (col %in% colnames(obj@meta.data)) {
    tsne <- DimPlot(obj, reduction = "tsne", label = TRUE, group.by = col) + xlab("t-SNE 1") + ylab("t-SNE 2") + 
      theme(axis.title = element_text(size = 18)) +  
      guides(colour = guide_legend(override.aes = list(size = 10)))
  } else {
    tsne <- ggplot() +
      theme_void() +
      geom_text(aes(x = 0.5, y = 0.5, label = "col doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  }
  #ggplotly(tsne)
  list(ggplot = tsne, plotly = ggplotly(tsne))
}

create_metadata_pca_3d <- function(obj, col){
  plotting.data <- FetchData(object = obj, vars = c("PC_1", "PC_2", "PC_3", "seurat_clusters"))
  plot_ly(data = plotting.data, x = ~PC_1, y = ~PC_2, z = ~PC_3, color = ~seurat_clusters,
          type = "scatter3d", mode = "markers")
}

create_metadata_tsne_3d <- function(obj, col){
  obj <- RunTSNE(obj, dims = 1:10, dim.embed = 3)
  
  # Extract tSNE information from Seurat Object
  tsne_1 <- obj[["tsne"]]@cell.embeddings[,1]
  tsne_2 <- obj[["tsne"]]@cell.embeddings[,2]
  tsne_3 <- obj[["tsne"]]@cell.embeddings[,3]
  
  # Prepare a dataframe for cell plotting
  plot.data <- FetchData(object = obj, vars = c("tSNE_1", "tSNE_2", "tSNE_3", "seurat_clusters"))
  
  # Make a column of row name identities (these will be your cell/barcode names)
  plot.data$label <- paste(rownames(plot.data))
  
  # Plot your data, in this example my Seurat object had 21 clusters (0-20)
  plot_ly(data = plot.data, 
          x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3, 
          color = ~seurat_clusters, 
          type = "scatter3d", 
          mode = "markers", 
          marker = list(size = 5, width=2), # controls size of points
          text=~label, #This is that extra column we made earlier for which we will use
          hoverinfo="text")
}

create_metadata_umap_3d <- function(obj, col){
  obj <- RunUMAP(obj, dims = 1:10, n.components = 3L)
  
  #Embeddings(object = pbmc, reduction = "umap")
  
  # Prepare a dataframe for cell plotting
  plot.data <- FetchData(object = obj, vars = c("umap_1", "umap_2", "umap_3", "seurat_clusters"))
  
  # Make a column of row name identities (these will be your cell/barcode names)
  plot.data$label <- paste(rownames(plot.data))
  
  # Plot your data, in this example my Seurat object had 21 clusters (0-20)
  umap_3d <- plot_ly(data = plot.data, 
                 x = ~umap_1, y = ~umap_2, z = ~umap_3, 
                 color = ~seurat_clusters, 
                 type = "scatter3d", 
                 mode = "markers", 
                 marker = list(size = 5, width=2), # controls size of points
                 text=~label, #This is that extra column we made earlier for which we will use for cell ID
                 hoverinfo="text")
  list(ggplot = umap_3d, plotly = ggplotly(umap_3d))
}

# create_feature_plot_tsne_hover <- function(obj, gene) {
#   if (gene %in% rownames(obj)) {
#     FP <- Seurat::FeaturePlot(obj, features = gene, pt.size = 0.001, combine = FALSE)
#   } else {
#     FP <- ggplot() +
#       theme_void() + 
#       geom_text(aes(x = 0.5, y = 0.5, label = "Gene doesn't exist"), size = 20, color = "gray73", fontface = "bold") +
#       theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#   }
#   return(FP)
# }

