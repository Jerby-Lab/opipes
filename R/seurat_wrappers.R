'%!in%' <- function(x,y)!('%in%'(x,y))

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

  so <- Seurat::CreateSeuratObject(counts, project = prj_name)
  so@meta.data$cellid <- row.names(so@meta.data)

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
std_preprocess <- function(so, nfeatures =100){

  so <- Seurat::NormalizeData(so)
  so <- Seurat::FindVariableFeatures(so, selection.method = "vst", nfeatures = 100)
  so <- Seurat::ScaleData(so)

  return(so)
}

#' Wrapper around PCA
#'
#' functionality for returning PCA results in Seurat Object
#'
#' @param so a SeuratObject of your data
#' @param pca.assay which assay to use, SCT or RNA
#' @param approx_pca.flag approximate PCA calculation?
#' @param verbose whether to be verbose with output
#' @export
#'
run_pca <- function(so,
                    approx_pca.flag = F,
                    pca.assay = "RNA",
                    n_pcs = 30,
                    verbose = F) {

  # run pca
  start_time <- Sys.time()
  print("Running PCA....")
  so <- Seurat::RunPCA(so,
                       npcs= n_pcs,
                       assay = pca.assay,
                       approx = approx_pca.flag,
                       verbose = verbose)
  return(so)

}

#' Variance Explained
#'
#' extract variance explained from Seurat Data
#'
#' @param so a SeuratObject of your data (with pca run)
#' @param slot which data slot to use
#' @return a list of all parameters calculated
#' @export

var_explained_seurat <- function(so, slot = "scale.data"){
  # calculate variance explained
  mat <- Seurat::GetAssayData(so, assay = pca.assay, slot = slot)
  total_variance <- sum(matrixStats::rowVars(mat))
  eigValues = (so[["pca"]]@stdev)^2  ## eigenValues
  varExplained = eigValues / total_variance

  return(varExplained)
}


#' Wrapper around UMAP and/or TSNE
#'
#' functionality for returning UMAP an TSNE results
#'
#' @param so a SeuratObject of your data
#' @param tsne.flag boolean, run tSNE?
#' @param umap.flag boolean, run umap?
#' @param n_dims number of dimensions you want to use from your embedding space
#' @param reduction string that matches to the right reduction slot
#' @param verbose whether to be verbose in the calculation
#' @return SeuratObject with computed visualizations
#' @export

umap_tsne <- function(so,
                  tsne.flag = F,
                  umap.flag = F,
                  n_dims = 30,
                  reduction = "pca",
                  verbose = F) {
  # run tsne
  if (tsne.flag){
    cat("\n\tRunning TSNE...\n")
    so <- Seurat::RunTSNE(object = so,
                          reduction.use = reduction,
                          dims = 1:n_dims,
                          do.fast = TRUE,
                          verbose=verbose)
  }

  # run umap
  if(umap.flag){
    cat("\n\tRunning UMAP...\n")
    so <- Seurat::RunUMAP(object = so,
                          dims = 1:n_dims,
                          reduction=reduction,
                          verbose = verbose)
  }

  return(so)
}


#' Wrapper around sNN clustering in Seurat
#'
#' functionality for returning clusters
#'
#' @param so a SeuratObject of your data
#' @param resolution clustering resolution, default = 0.4
#' @param n_pcs how many pcs to cluster on
#' @return a list of all parameters calculated
#' @export
#'
cluster <- function(so,
                    resolution = 0.4,
                    n_dims = 30,
                    reduction = "pca"){
  # clustering
  print("Running Clustering...")
  so <- Seurat::FindNeighbors(so, reduction = reduction ,dims = 1:n_dims)
  so <- Seurat::FindClusters(so, resolution = resolution)

  return(so)
}



















