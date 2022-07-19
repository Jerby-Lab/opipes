#
#
# rds <- readRDS("/Volumes/ljerby/HGSC_Profiling/Data/R_data/SMI_HGSC_malignant_hotScores.rds")
# rds <- readRDS("/Volumes/ljerby/HGSC_Profiling/Data/R_data/SMI_HGSC_hotInSituPlotColors.rds")
#
# rgb2hex <- function(x){
#   hex = rgb(x[1], x[2], x[3], maxColorValue=255)
#   return(hex)
# }
#
# spatial_sample_visualization <- function(
#   seg_path,
#   celltypes,
#   cell2rgb,
#   cont_df,
#   continuous_field = c("Malignant"),
#   out_path) {
#
#   cellseg = read.csv(seg_path)
#   colnames(cellseg) <- unlist(lapply(colnames(cellseg), function(x) {strsplit(x , split = "X0.")[[1]][2]}))
#   colnames(cellseg)[1] <- "0"
#   cellmask <- cellseg
#   cellmask[cellmask != 0] <- 0.01
#   cellmask[cellmask == 0] <- 1
#   cellmask <- EBImage::Image(as.matrix(cellmask))
#   EBImage::display(cellmask)
#
#   # unique(cell_types)
#   #
#   # for (type in levels) {
#   #   print(type)
#   #   cellids = rds$cells[rds$samples == sample_name & rds$cell.types2 == type]
#   #   print(length(cellids))
#   #   cellids = unlist(lapply(cellids, function(x) as.integer(strsplit(x, split ="c")[[1]][2])))
#   #   celltype_mask <- t(apply(cellseg, 1, function(x){x %in% cellids}))
#   #
#   #   # change the channels
#   #   celltype_rgb_1 <- celltypes@.Data[,,1]
#   #   celltype_rgb_1[celltype_mask] <- as.numeric(cell_2_rgb[[type]][1])/255
#   #   celltypes[,,1] <- celltype_rgb_1
#   #   celltype_rgb_2 <- celltypes@.Data[,,2]
#   #   celltype_rgb_2[celltype_mask] <- as.numeric(cell_2_rgb[[type]][2])/255
#   #   celltypes[,,2] <- celltype_rgb_2
#   #   celltype_rgb_3 <- celltypes@.Data[,,3]
#   #   celltype_rgb_3[celltype_mask] <- as.numeric(cell_2_rgb[[type]][3])/255
#   #   celltypes[,,3] <- celltype_rgb_3
#   # }
#   # writeImage(celltypes, paste0(OUTPATH, DATE, "_", name, "_celltypes.jpg"), quality = 85)
# }
#
#
# seg_path = paste0("/Volumes/ljerby/HGSC_Profiling/Results/Segmentation/SMI/Seg/",
#                           "SMI_T10_F001",
#                           "_whole-cell_03.csv")
#
# cell2rgb <- list("Other" = c(200, 200, 200),
#               "Malignant" = c(7, 224, 0),
#               "T.cell" = c(0, 0, 0))
#
# cell2rgb <- list("Other" = c(200, 200, 200),
#                  "Malignant" = c(7, 224, 0),
#                  "T.cell" = c(255, 247, 0))
#
# r <- readRDS("/Volumes/ljerby/HGSC_Profiling/Data/R_data/SMI_HGSC_seg03_annotated_20220712.rds")
#
# q <- r$cells[r$samples == "SMI_T10_F001"]
#
# spatial_sample_visualization(seg_path, NULL, cell2rgb, NULL, NULL, NULL)
