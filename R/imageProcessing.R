#' RGB to Hexcode
#'
#' Simple function that converts a length 3 list corresponding
#' to RGB code to a color hexcode
#'
#' @param x length three string corresponding to RGB color code
#' @returns  a hexcode (string)
#' @export
rgb2hex <- function(x){
  hex = rgb(x[1], x[2], x[3], maxColorValue=255)
  return(hex)
}

#' Spatial Cell Mask Visualization
#'
#' Plots for each sample, an image of each segmented
#' cell colored by user specified features
#'
#' @param seg_path string of full path to cell segmentation map
#' @param celltypes a list of cell labels with names of the list corresponding to cellid
#' @param cell2rgb a dictionary that matches discrete cell label to a rgb code
#' @param samplename string name of SMI sample
#' @param background (default "black") background color ("black" or "white")
#' @param outpath (default = "~/") string path to output directory
#' @param file_prefix (default = "") string prefix of file
#' @param cont_field the celltype for which the coloring is continuous, not discrete. names must correspond to cellid
#' @param contvals list of the hexcolors corresponding to each cell to be plotted on continuous scale
#' @export
spatial_sample_visualization <- function(
  seg_path,
  celltypes,
  cell2rgb,
  samplename,
  background = "black",
  outpath = "~/",
  file_prefix = "",
  cont_field = "",
  contvals = NULL) {

  # read in a process cellsegmentation mask
  cellseg = read.csv(seg_path)
  colnames(cellseg) <- unlist(lapply(colnames(cellseg), function(x) {strsplit(x , split = "X0.")[[1]][2]}))
  colnames(cellseg)[1] <- "0"
  cellmask <- cellseg
  cellmask[cellmask != 0] <- 0.8
  if (background == "white") {
    cellmask[cellmask == 0] <- 1
  } else {
    cellmask[cellmask == 0] <- 0
  }

  # convert cellsegmentation mask to EBImage object
  cellmask <- EBImage::Image(as.matrix(cellmask))
  cellmask <- EBImage::channel(cellmask, 'rgb')

  # set discrete colors
  levels <- setdiff(unique(celltypes), cont_field)
  for (type in levels) {
    print(type)
    cellids = names(celltypes[celltypes == type])
    print(length(cellids))
    cellidx = unlist(lapply(cellids, function(x) as.integer(strsplit(x, split ="c")[[1]][2])))
    celltype_mask <- t(apply(cellseg, 1, function(x){x %in% cellidx}))

    # change the channels
    celltype_rgb_1 <- cellmask@.Data[,,1]
    celltype_rgb_1[celltype_mask] <- as.numeric(cell2rgb[[type]][1])/255
    cellmask[,,1] <- celltype_rgb_1
    celltype_rgb_2 <- cellmask@.Data[,,2]
    celltype_rgb_2[celltype_mask] <- as.numeric(cell2rgb[[type]][2])/255
    cellmask[,,2] <- celltype_rgb_2
    celltype_rgb_3 <- cellmask@.Data[,,3]
    celltype_rgb_3[celltype_mask] <- as.numeric(cell2rgb[[type]][3])/255
    cellmask[,,3] <- celltype_rgb_3
  }

  # set continuous colors
  if (cont_field != "") {
    cellids = names(celltypes[celltypes == cont_field])
    print(cont_field)
    print(length(cellids))
    cellidx = unlist(lapply(cellids, function(x) as.integer(strsplit(x, split ="c")[[1]][2])))
    celltype_mask <- t(apply(cellseg, 1, function(x){x %in% cellidx}))
    cellids <- paste0(samplename, "_c", cellseg[celltype_mask])

    colormap <- t(sapply(unique(contvals), simplify = T, grDevices::col2rgb))/255

    celltype_rgb_1 <- cellmask@.Data[,,1]
    celltype_rgb_1[celltype_mask] <- colormap[contvals[cellids],1]
    cellmask[,,1] <- celltype_rgb_1
    celltype_rgb_2 <- cellmask@.Data[,,2]
    celltype_rgb_2[celltype_mask] <- colormap[contvals[cellids],2]
    cellmask[,,2] <- celltype_rgb_2
    celltype_rgb_3 <- cellmask@.Data[,,3]
    celltype_rgb_3[celltype_mask] <- colormap[contvals[cellids],3]
    cellmask[,,3] <- celltype_rgb_3
  }

  # write to file
  outfile = paste0(outpath, "/", file_prefix, "_", samplename, ".jpg")
  EBImage::writeImage(cellmask, files = outfile)
}
