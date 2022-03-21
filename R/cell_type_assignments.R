#' scRNA Marker Based Cell Type Assignments
#'
#' Using gene expression markers, unsupervised clusters,
#' and gene expression, assign the likely celltype to those
#' clusters
#'
#' @param r data list including the gene expression data
#' @param cell.clusters a _n_-length list to where _n_ is the number of cells
#' @param cell.sig cell type signatures (list of lists of markers associated with each cell-type, default exists)
#' @param immune.markers immune markers (default exists)
#' @param non.immune.cell.types which cell types are not immune cell types?
#' @param EM.flag run EM? (default = F)
#' @param OE.type which version of overall-expression to calculate (default = "V1")
#' @param test.type which kind of test to test marker significance between clusters
#' @param minZ minimum significance cutoff to assign a cell-type (default = 10)
#' @return a list with objects named "c" and "m" corresponding to your subsampled counts and metadata
#' @export

scRNA_markers_assign <- function(r,
                                 cell.clusters,
                                 cell.sig,
                                 immune.markers,
                                 non.immune.cell.types = c('Endothelial','CAF','Malignant'),
                                 EM.flag = F,
                                 OE.type = "V1",
                                 test.type = "ttest",
                                 minZ = 10){

  # -------------------------------------
  # reconcile cell-type signature markers
  # -------------------------------------
  if(missing(cell.sig)){
    cell.sig<-list(T.cell = c("CD2","CD3D","CD3E","CD3G","CD4","CD8A","CD8B"),
                   B.cell = c("CD19","CD79A","CD79B","BLK","CD20","CD22"),
                   Macrophage = c("CD163","CD14","CSF1R","CD4"),
                   Mast = c("ENPP3","KIT"),
                   NK = c("KLRA1","NKG2","KLRB1","KLRD1"),
                   Endothelial =  c("PECAM1", "VWF", "CDH5"),
                   CAF = c("FAP","THY1","DCN", "COL1A1","COL1A2","COL6A1","COL6A2","COL6A3"),
                   Malignant = c("EPCAM","MUC16","CD24"))
  }
  cell.sig<-intersect.list1(cell.sig,r$genes,n1 = 1)

  # ----------------------------------
  # reconcile immune cell type markers
  # ----------------------------------
  if(missing(immune.markers)){
    immune.markers<-c("CD3D","CD3E","CD3G","PTPRC","CD8A","CD8B","CD4","CD79A","CD19")
  }

  # ---------------------------------------------------
  # double check cell markers for non-immune cell types
  # ---------------------------------------------------
  non.immune.cell.types<-intersect(names(cell.sig),non.immune.cell.types)
  if(any(!is.element(non.immune.cell.types,names(cell.sig)))){
    err.msg<-"Error: There are some non-immune cell types without marker genes."
    print(err.msg)
    return(r)
  }

  # --------------------------------------
  # double check clusters are in the input
  # --------------------------------------
  if(missing(cell.clusters)){
    err.msg<-"Error: Cluster information is missing."
    print(err.msg)
    return(r)
  }

  # ----------------------------------------------------
  # identify cells that express canonical immune markers
  # ----------------------------------------------------
  r$b.immune<-colSums(r$tpm[is.element(r$genes,immune.markers),]>0)>0

  # ----------------------------
  # calculate overall expression
  # ----------------------------
  if(is.null(r$binZ)){r<-prep4OE(r)}

  if(OE.type == "V1"){
    markers.scores<-get.OE(r, cell.sig) #c(cell.sig,list(cc = cc.genes)))
    r$markers.scores<-markers.scores[,names(cell.sig)]
    # r$cc<-markers.scores[,"cc"]
  }else{
    r$zscores1<-t(apply(r$tpm,1,function(x) x/max(x)))
    r$markers.scores<-t(laply(cell.sig,function(x) colMeans(r$zscores1[intersect(x,r$genes),])))
    # r$markers.scores<-t(laply(cell.sig,function(x) colMeans(r$tpm[intersect(x,r$genes),])))
    colnames(r$markers.scores)<-names(cell.sig)
    r$markers.scores<-10*center.matrix(r$markers.scores)
  }

  # -------------------------
  # determine CD8/CD4 T-cells
  # -------------------------
  r$cd8<-rep(F,length(r$cells))
  r$cd4<-rep(F,length(r$cells))
  if(is.element("CD8A",rownames(r$tpm))){r$cd8<-r$cd8|r$tpm["CD8A",]>0}
  if(is.element("CD8B",rownames(r$tpm))){r$cd8<-r$cd8|r$tpm["CD8B",]>0}
  if(is.element("CD4",rownames(r$tpm))){r$cd4<-r$tpm["CD4",]>0}

  if(is.element("Cd8a",rownames(r$tpm))){r$cd8<-r$cd8|r$tpm["Cd8a",]>0}
  if(is.element("Cd8b1",rownames(r$tpm))){r$cd8<-r$cd8|r$tpm["Cd8b1",]>0}
  if(is.element("Cd4",rownames(r$tpm))){r$cd4<-r$tpm["Cd4",]>0}

  # ----------------------------
  # compute binary marker scores
  # ----------------------------
  r$markersB<-apply(r$markers.scores,2,function(x1){get.binary.scores(x1, EM.flag)})
  b.non.immune<-rowSums(as.matrix(r$markersB[,non.immune.cell.types]))>0
  r$markersB<-cbind.data.frame(r$markersB,
                               Doublet = b.non.immune&r$b.immune,
                               T.CD8 = r$markersB[,"T.cell"]&r$cd8&!r$cd4,
                               T.CD4 = r$markersB[,"T.cell"]&r$cd4&!r$cd8)
  r$markersB[r$b.immune,c('Endothelial','CAF','Malignant')]<-F
  r$markersB[is.na(r$markersB)]<-0
  b<-r$markersB$malignant
  q1<-0.5

  # --------------------------------------------
  # account for CNAs to identify malignant cells.
  # --------------------------------------------
  if(!is.null(r$cnv.scores)){
    r$cnv.scores$top<-r$cnv.scores$mean>quantile(r$cnv.scores$mean,q1)&r$cnv.scores$R>quantile(r$cnv.scores$R,q1,na.rm = T)
    r$markersB$malignant<-plot.bimodal.distribution(r$cnv.scores$R)$labels
    r$markersB[r$markersB$malignant,setdiff(colnames(r$markersB),"malignant")]<-F
  }

  # --------------------------------------------
  # statistical tests for determining cell-type
  # --------------------------------------------
  if(test.type == "ttest"){
    # Use t-test to identify if the cells in a particular cluster
    # have a higher expression of a specific cell type signature
    # compared to all the other clusters together
    r<-scRNA_cluster.annotation.ttest(r,cell.clusters = cell.clusters,minZ = minZ)
  }else{
    # Use enrichment of a specific cell type annotation in a cluster to annotate clusters.
    r<-scRNA_cluster.annotation.HG(r,cell.clusters = cell.clusters,EM.flag = EM.flag)
  }

  # ----------------------------------------------
  # subtype annotations based on canonical markers
  # ----------------------------------------------
  r$cell.types[r$cell.types=="T.CD8"&(!r$cd8&r$cd4)]<-"T.cell.UD"
  r$cell.types[r$cell.types=="T.CD4"&(r$cd8&!r$cd4)]<-"T.cell.UD"
  r$cell.types[r$cell.types=="NK"&(r$cd8|r$cd4)]<-"T.cell.UD"
  r$cell.types[r$cell.types==""]<-"UD"

  # -------------------------
  # print table of cell-types
  # -------------------------
  print(table(r$cell.types))

  return(r)
}

#' Hypergeometric Test
#'
#' function for computing hypergeometric test
#' on single cell clustered data
#'
#' @param r data list including the gene expression data
#' @param cell.clusters a _n_-length list of cluster assignments, where _n_ = number of cells
#' @param EM.flag whether to compute EM
#' @return data list, but now including a `cell.types` slot
#'
scRNA_cluster.annotation.HG <- function(r,
                                        cell.clusters,
                                        EM.flag){

  B.clusters<-labels.mat.2.logical.mat(cell.clusters)
  # b<-rowSums(r$markersB[,names(cell.sig)])==1
  # P<-get.hyper.p.value.mat(B.clusters[b,],r$markersB[b,],full.flag = T)
  P<-get.hyper.p.value.mat(B.clusters,r$markersB,full.flag = T)
  P$p[is.na(P$p)]<-1
  View(t(P$frc))
  cluster.cell.types<-laply(rownames(P$p),function(x){
    p<-P$p[x,];sort(p)
    ob<-P$ob[x,]
    frc<-P$frc[x,]
    b1<-p==min(p)&p<1e-3&frc>ifelse(EM.flag,0.6,0.08)
    if(!any(b1)|sum(P$frc[x,names(cell.sig)]>0.4)>1){return("UD")}
    b2<-ob==max(ob[b1])
    v1<-names(ob)[b1&b2]
    v<-v1
    # v2<-colnames(r$markersB)[colMeans(r$markersB[B.clusters[,x],])>0.3]
    # v<-paste(intersect(v1,v2),collapse = "_")
    return(v)
  })
  names(cluster.cell.types)<-rownames(P$p)
  r$cell.types<-cluster.cell.types[cell.clusters]
  barplot(sort(table(cell.types)/length(cell.types),decreasing = T),las=2,main = "HG-based")
  return(r)
}

#' Compute _t_-tests
#'
#' function for computing two sample _t_-test
#' on single cell clustered data
#'
#' @param r data list including the gene expression data
#' @param cell.clusters a _n_-length list of cluster assignments, where _n_ = number of cells
#' @param minZ cutoff for significance in the -log10(p) space (default = 10)
#' @return data list, but now including `cell.types`, and `assignments.ttests` slots
#'
scRNA_cluster.annotation.ttest <- function(r,
                                         cell.clusters,
                                         minZ = 10) {
  b<-sample.per.label(r$clusters,500) # boolean for subsampling
  z<-plyr::laply(unique(r$clusters),function(x) t.test.mat(t(r$markers.scores[b,]),r$clusters[b]==x)[,3])
  rownames(z)<-unique(r$clusters)
  z<-cbind.data.frame(z,cell.type1 = colnames(z)[apply(z,1,function(x) which(x==max(x))[1])],
                      n = rowSums(z>minZ))
  b2<-z$n>1
  z$diff<-apply(z[,1:(ncol(z)-2)],1,function(x) sort(x,decreasing = T)[1])-apply(z[,1:(ncol(z)-2)],1,function(x) sort(x,decreasing = T)[2])
  z$unique<-z$n==1|z$diff>10
  z$cell.type<-z$cell.type1
  z$cell.type[!z$unique]<-"UD"
  r$cell.types<-z[r$clusters,"cell.type"]
  b.mal<-r$cell.types=="Malignant"
  b<-b&!b.mal
  z.mal<-plyr::laply(unique(r$clusters),function(x) t.test.mat(t(r$markers.scores[b,]),r$clusters[b]==x)[,3])
  # if(!is.null(r$cc)){
  #   z.cc<-laply(unique(r$clusters),function(x){
  #     p<-t.test.labels(r$cc,r$clusters==x,alternative = "greater")
  #     p<-max(p,1e-30)
  #     z<-ifelse(p<0.5,-log10(p),log10(1-p))
  #     return(z)
  #   })
  # }
  z$Malignant2<-z.mal[,"Malignant"]
  # z$cycling<-z.cc
  z$cell.type[z$cell.type=="UD"&z$n==0&z$Malignant2>minZ]<-"Malignant"
  # z$cell.type[z$cell.type=="UD"&z$n==0&z$cycling>minZ]<-"Cycling"
  r$cell.types<-z[r$clusters,"cell.type"]
  r$assignments.ttests<-z
  return(r)
}

#' Get Binary Scores
#'
#' function for getting binary scores (tends to be more relevant)
#' for the EM algorithm
#'
#' @param x1 list of the markers score for a particular cell-type
#' @return a list of binary scores
#'
get.binary.scores<-function(x1, EM.flag){
  if(!EM.flag){
    return(x1>quantile(x1,0.95))
  }
  b<-x1>quantile(x1,0.2);x<-x1[b]
  f1<-function(x,seed){
    v<-tryCatch({plot.bimodal.distribution(x,seed = seed)},
                error = function(err){
                  return(list(loglik = 1000,labels = rep(NA,length(x))))})
  }
  v1<-f1(x,111);v2<-f1(x,123);v3<-f1(x,222);v<-v1
  v4<-f1(x,493369196)
  if(v2$loglik>v$loglik){v<-v2}
  if(v3$loglik>v$loglik){v<-v3}
  if(v4$loglik>v$loglik){v<-v4}
  b[b]<-v$labels
  return(b)
}

#' Average Matrix Rows
#'
#' function for computing matrix averages
#'
#' @param m matrix
#' @param ids ids of the samples you care about
#' @param f aggregation function
#' @return a list of binary scores
#'
average.mat.rows<-function(m,ids,f = colMeans){
  ids.u<-sort(unique(ids))
  m1<-get.mat(ids.u,colnames(m))
  for(x in ids.u){
    b<-is.element(ids,x)
    if(sum(b)==1){
      m1[x,]<-m[b,]
    }else{
      m1[x,]<-f(m[b,])
    }
  }
  return(m1)
}

#' Get Matrix
#'
#' function for casting data into a matrix.
#'
#' @param m.rows matrix rows
#' @param m.cols matrix columns
#' @param data your data
#' @return your data in matrix format
#'
get.mat<-function(m.rows,m.cols,data = NA){
  m<-matrix(data = data, nrow = length(m.rows),ncol = length(m.cols),
            dimnames = list(m.rows,m.cols))
  return(m)
}
