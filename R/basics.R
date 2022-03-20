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

subsample <- function(counts, metadata, p=0.5, replace = T){
  n_cells <- dim(counts)[2]
  cell_idx <- sample(1:n_cells, size = p*n_cells, replace = replace)
  subsample_counts <- counts[, cell_idx]
  subsample_meta <- metadata[cell_idx,]
  out <- list()
  out[["c"]] <- subsample_counts
  out[["m"]] <- subsample_meta
  return(out)
}


#' Reconstruction Loss
#'
#' functionality for returning reconstruction loss
#' for any dimension reduction algorithm. Please make sure
#' that the the order of your rows and columns are the same!
#'
#' @param X original matrix _n_ x _m_ matrix
#' @param X_r low-dimensional reconstructed _n_ x _m_ matrix
#' @return a list of reconstruction loss, labeled with cellnames
#' @export
#'
recon_loss <- function(X, X_r){
  recon_loss <- unlist(lapply(1:dim(X)[1], function(i){sqrt(sum((X[i,] - X_r[i,])^2))}))
  names(recon_loss) <- row.names(X_r)
  return(recon_loss)
}
