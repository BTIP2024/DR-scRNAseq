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
    FP <- Seurat::FeaturePlot(obj, features = gene, pt.size = 10, combine = FALSE, reduction = "pca")
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
  ggplotly(pca)
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
  ggplotly(umap)
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
  ggplotly(tsne)
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

