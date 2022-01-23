'%!in%' <- function(x,y)!('%in%'(x,y))

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

#' Wrapper around PCA 
#' 
#' functionality for returning PCA results in Seurat Object
#'
#' @param so a SeuratObject of your data 
#' @param pca.assay = which assay to use, SCT or RNA
#' @param approx_pca.flag approximate PCA calculation? 
#' @export
#' 
run_pca <- function(so,
                    approx_pca.flag = F,
                    pca.assay = "RNA",
                    n_pcs = 30) {
  
  out <- list()
  
  # run pca
  start_time <- Sys.time()
  print("Running PCA....")
  so <- Seurat::RunPCA(so, 
                       npcs= length(row.names(so)), 
                       assay = pca.assay, 
                       approx = approx_pca.flag) # always run on all PCs
  
  # calculate variance explained
  mat <- Seurat::GetAssayData(so, assay = pca.assay, slot = "scale.data")
  total_variance <- sum(matrixStats::rowVars(mat))
  eigValues = (so[["pca"]]@stdev)^2  ## EigenValues
  varExplained = eigValues / total_variance
  
  print(paste0("PC1 Variance Explained: ", round(varExplained[1]*100, 2), "%"))
  print(paste0("PC2 Vairance Explained: ", round(varExplained[2]*100, 2), "%"))
  print(paste0("Variance Explained in ", n_pcs, " PCs : ", round(sum(varExplained[1:n_pcs]*100), 2), "%"))
  
  end_time <- Sys.time()
  run_time <- end_time - start_time
  print(run_time)
  
  
  out[["pca"]] <- so[["pca"]]@cell.embeddings
  out[["so"]] <- so
  
  return(out)
  
}


#' Wrapper around UMAP and/or TSNE
#' 
#' functionality for returning UMAP an TSNE results 
#'
#' @param so a SeuratObject of your data
#' @param tsne.flag boolean, run tSNE? 
#' @param umap.flag boolean, run umap?
#' @param n_pcs number of PCs to run umap and tsne on 
#' @return a list of all parameters calculated 

embed <- function(so, 
                  tsne.flag = F,
                  umap.flag = F,
                  n_pcs = 30) {
  out <- list()
  
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
  
  
  out[["so"]] <- so
  
  return(out)
}


#' Wrapper around sNN clustering in Seurat
#' 
#' functionality for returning clusters 
#'
#' @param so a SeuratObject of your data
#' @param resolution clustering resolution, default = 0.4
#' @return a list of all parameters calculated 

cluster <- function(so, 
                    resolution = 0.4){
  out <- list()
  
  # clustering
  print("Running Clustering...")
  start_time <- Sys.time()
  so <- Seurat::FindNeighbors(so, reduction = "pca",dims = 1:n_pcs_umap)
  so <- Seurat::FindClusters(so, resolution = resolution)
  end_time <- Sys.time()
  run_time <- end_time - start_time
  print(run_time)
  
  out[["cluster"]] <- so@meta.data$seurat_clusters
  out[["so"]] <- so
  
  return(out)
}












                  

                  
                  
                  
        



