#' Randomly subsample single cell data set
#' 
#' Randomly subsample a proportion of data, 
#' cells must be in the same order in counts 
#' and metadata. 
#'
#' @param counts a _m_ x _n_ matrix of raw count data where _m_ = no. genes, and _n_ = no. cells 
#' @param metadata matrix _n_ x _j_ of metadata where _n_ = no. cells, and _j_ = no. features
#' @param p proportion of data to subsample (p = [0, 1])
#' @return a list with objects named "c" and "m" corresponding to your subsampled counts and metadata 
#' @export

subsample <- function(counts, metadata, p=0.5){
  n_cells <- dim(counts)[2]
  cell_idx <- sample(1:n_cells, size = p*n_cells)
  subsample_counts <- counts[, cell_idx]
  subsample_meta <- metadata[cell_idx,]
  out <- list()
  out[["c"]] <- subsample_counts
  out[["m"]] <- subsample_meta
  return(out)
}

#' Seuratify your Single Cell Data 
#' 
#' This function takes in raw counts data and 
#' associated metadata and creates a SeuratObject
#'
#' @param counts a _m_ x _n_ matrix of raw count data where _m_ = no. genes, and _n_ = no. cells 
#' @param metadata matrix of metadata with at least a column cell identifier named "cellid" 
#' @param prj_name project name string 
#' @return a Seurat object of your count matrix with metadata 
#' @export
seuratify <- function(counts, prj_name = ""){
  start_time <- Sys.time()
  
  so <- Seurat::CreateSeuratObject(counts, project = prj_name)
  so@meta.data$cellid <- row.names(so@meta.data)
  
  end_time <- Sys.time()
  run_time <- end_time - start_time
  print(run_time)
  
  return(so)
}

#' Standard single cell RNA preprocessing 
#' 
#' normalize, transform, and scale
#'
#' @param so a SeuratObject of your data 
#' @return a preprocessed SeuratObject
#' @export
#' 
std_preprocess <- function(so){
  start_time <- Sys.time()

  so <- Seurat::NormalizeData(so)
  so <- Seurat::FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
  so <- Seurat::ScaleData(so)
  
  end_time <- Sys.time()
  run_time <- end_time - start_time
  print(run_time)
  return(so)
}

#' Embedding methods: PCA, tSNE, UMAP
#' 
#' functionality for returning PCA, tSNE, 
#' and UMAP projections 
#'
#' @param so a SeuratObject of your data 
#' @param approx_pca.flag approximate PCA calculation? 
#' @param tsne.flag boolean, run tSNE? 
#' @param umap.flag boolean, run umap?
#' @param n_pcs number of PCs to run umap and tsne on 
#' @param resolution clustering resolution, default = 0.4
#' @param plot plot initial clustering?
#' @return a list of all parameters calculated 
#' @export
#' 
embed <- function(so,
                  approx_pca.flag = T,
                  tsne.flag = F,
                  umap.flag = F,
                  n_pcs = 30,
                  resolution = 0.4,
                  plot = T, 
                  plot.out = ""){
  out <- list()
  
  # run pca
  start_time <- Sys.time()
  print("Running PCA....")
  so <- Seurat::RunPCA(so, npcs= n_pcs, approx = approx_pca.flag)
  
  # calculate variance explained
  mat <- Seurat::GetAssayData(so, assay = "RNA", slot = "scale.data")
  total_variance <- sum(matrixStats::rowVars(mat))
  eigValues = (so[["pca"]]@stdev)^2  ## EigenValues
  varExplained = eigValues / total_variance
  
  print(paste0("PC1 Variance Explained: ", round(varExplained[1]*100, 2), "%"))
  print(paste0("PC2 Vairance Explained: ", round(varExplained[2]*100, 2), "%"))
  print(paste0("Variance Explained in ", n_pcs, " PCs : ", round(sum(varExplained*100), 2), "%"))
  
  end_time <- Sys.time()
  run_time <- end_time - start_time
  print(run_time)
  
  out[["pca"]] <- so[["pca"]]@cell.embeddings
  
  # run tsne
  if (tsne.flag){
    start_time <- Sys.time()
    
    print("Running TSNE...")
    so <- Seurat::RunTSNE(object = so, reduction.use = "pca",dims = 1:n_pcs, do.fast = TRUE)
     out[["tsne"]] <- so[["tsne"]]@cell.embeddings
    
    end_time <- Sys.time()
    run_time <- end_time - start_time
    print(run_time)
  }
  
  # run umap 
  if(umap.flag){
    start_time <- Sys.time()
    
    print("Running UMAP...")
    so <- Seurat::RunUMAP(object = so, dims = 1:n_pcs, )
    out[["umap"]] <- so[["umap"]]@cell.embeddings
    
    end_time <- Sys.time()
    run_time <- end_time - start_time
    print(run_time)
  }
  
  # clustering
  print("Running Clustering...")
  start_time <- Sys.time()
  so <- Seurat::FindNeighbors(so, reduction = "pca",dims = 1:n_pcs)
  so <- Seurat::FindClusters(so, resolution = resolution)
  end_time <- Sys.time()
  run_time <- end_time - start_time
  print(run_time)
  
  out[["cluster"]] <- so@meta.data$seurat_clusters
  
  if (plot) {
    if (tsne.flag) {
      p <- Seurat::DimPlot(so, reduction = "tsne", group.by = "RNA_snn_res.0.4")
      png(paste0(plot.out, "tsne.png"))
      print(p)
      dev.off()
    }
    if (umap.flag) {
      p <- Seurat::DimPlot(so, reduction = "umap", group.by = "RNA_snn_res.0.4")
      png(paste0(plot.out, "umap.png"), width = 6, height = 5, res = 200, units = "in")
      print(p)
      dev.off()
    }
  }
  
  out[["so"]] <- so
  
  return(out)
}
