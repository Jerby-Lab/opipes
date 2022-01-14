subsample <- function(counts, metadata){
  counts
}

#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param counts a _m_ x _n_ matrix of raw count data where _m_ = no. genes, and _n_ = no. cells 
#' @param prj_name project name string 
#' @param metadata matrix of metadata with at least a column cell identifier named "cellid" 
#' @return a Seurat object of your count matrix with metadata 
#' @export
seuratify <- function(counts, prj_name, metadata){
  so <- Seurat::CreateSeuratObject(counts, project = prj_name)
  so@meta.data$cellid <- row.names(so@meta.data)
  so@meta.data <- dplyr::merge(so@meta.data, metadata, by.x = "cellid", all.x = T, sort = F)
  return(so)
}

# std_preprocess

# embed




