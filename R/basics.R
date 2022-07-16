#' Randomly subsample single cell data set
#'
#' Randomly subsample a proportion of data,
#' cells must be in the same order in counts
#' and metadata.
#'
#' @param counts a _m_ x _n_ matrix of raw count data where _m_ = no. genes, and _n_ = no. cells
#' @param metadata matrix _n_ x _j_ of metadata where _n_ = no. cells, and _j_ = no. features
#' @param p proportion of data to subsample (p = [0, 1])
#' @param replace boolean to sample with (T) or without (F) replacement
#' @return a list with objects named "c" and "m" corresponding to your subsampled counts and metadata
#' @export
subsample <- function(counts,
                      metadata,
                      p=0.5,
                      replace = F){
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


#' Sample cells per label
#'
#' samples _size_ number of cells from each
#' labelled population
#'
#' @param labels list of cells based on labels
#' @param size number of cells to sample per label
#' @param boolean.flag not sure
#' @param v not sure what this is
#' @param remove.flag not sure
#' @return a list of cell names that correspond to the subsample
#'
sample.per.label<-function(labels,
                           size,
                           boolean.flag = T,
                           v,
                           remove.flag = F){
  if(missing(v)){
    v<-1:length(labels)
  }
  ul<-unique(labels)
  vr<-NULL
  for(i in 1:length(ul)){
    b<-is.element(labels,ul[i])
    if(size<1){
      vr<-c(vr,sample(v[b],size = sum(b)*size))
    }else{
      vr<-c(vr,sample(v[b],size = min(size,sum(b))))
    }
  }
  b<-is.element(v,vr)
  if(remove.flag){
    b.small<-!is.element(labels,get.abundant(labels[b],size-1))
    b[b.small]<-F
    vr<-v[b]
  }
  if(boolean.flag){return(b)}
  return(vr)
}

#' Convert p-values to z-scores
#'
#' samples _size_ number of cells from each
#' labelled population
#'
#' @param p a table of p values where p[,1] is the greater than hypothesis, and the p[,2] is the less than hypothesis
#' @return _z_-scores
#'
get.p.zscores<-function(p){
  b<-p[,1]>0.5
  b[is.na(b)]<-F
  zscores<-(-log10(p[,1]))
  zscores[b]<-log10(p[b,2])
  # signficiant in p[,1] will be positive
  # signficiant in p[,2] will be negative
  return(zscores)
}

#' Intersect Lists of Lists
#'
#' functionality for finding the intersection
#' of lists of lists.
#'
#' @param l first list of lists
#' @param g second lists of lists
#' @param n1 cutoff for length of list
#' @param HG.universe genes for computing GO enrichment
#' @param prf some sort of string parsing
#' @return _z_-scores
#' @export
#'
intersect.list1<-function(l,
                          g,
                          n1=0,
                          HG.universe = NULL,
                          prf = ""){
  l1<-lapply(l, function(x) intersect(x,g))
  l1<-l1[plyr::laply(l1,length)>n1]
  if(prf!=""){
    names(l1)<-paste(prf,names(l1),sep = ".")
  }
  if(!is.null(HG.universe)){
    p<-GO.enrichment.lapply(l[names(l1)],genes = HG.universe,list(g))
    names(l1)<-paste(names(l1),format(p,scientific = T,digits= 3),sep = "P = ")
  }
  return(l1)
}

#' Helper function for Overall Expression
#'
#' prepares data to compute overall expression
#' signature for each cell-type
#'
#' @param r data list
#' @param n.cat (default = 50)
#' @return r, your data list
#'
prep4OE<-function(r,
                  n.cat = 50){

  r$zscores<-center.matrix(r$tpm,dim = 1,sd.flag = T)
  X<-10*((2^r$tpm)-1)
  r$genes.dist<-log2(rowMeans(X,na.rm = T)+1)
  r$genes.dist.q<-discretize.prvt(r$genes.dist,n.cat = n.cat)
  b<-rowSums(is.na(r$zscores))==0
  if(any(!b)){r<-set.list(r,b)}
  r$binZ<-average.mat.rows(r$zscores,r$genes.dist.q,f = colMeans)
  return(r)
}

#' Overall Expression Wrapper
#'
#' wrapper around code that computes overall expression
#' signature for each cell-type
#'
#' @param r data list
#' @param sig (default = 50)
#' @return r, your data list
#' @export
get.OE<-function(r,sig){
  scores<-get.OE1(r,sig)
  names(sig)<-gsub(" ",".",names(sig))
  two.sided<-unique(gsub(".up","",gsub(".down","",names(sig))))
  b<-is.element(paste0(two.sided,".up"),names(sig))&
    is.element(paste0(two.sided,".down"),names(sig))
  if(any(b)){
    two.sided<-two.sided[b]
    scores2<-as.matrix(scores[,paste0(two.sided,".up")]-scores[,paste0(two.sided,".down")])
    colnames(scores2)<-two.sided
    scores<-cbind(scores2,scores)
  }

  if(!is.null(r$cells)){
    rownames(scores)<-r$cells
  }else{
    if(!is.null(r$samples)){
      rownames(scores)<-r$samples
    }
  }
  return(scores)
}

#' Overall Expression Subroutine
#'
#' code that computes overall expression
#' signature for each cell-type
#'
#' @param r data list
#' @param sig (default = 50)
#' @return r, your data list
<<<<<<< HEAD
#' @export
=======
#'
>>>>>>> 8f781af8a6374298c6baaf805667401e5c9a92fb
get.OE1 <- function(r,sig){
  if(is.list(sig)){
    scores<-t(plyr::laply(sig, function(g) get.OE1(r,g)))
    rownames(scores)<-r$cells
    colnames(scores)<-names(sig)
    return(scores)
  }
  g<-sig
  b<-is.element(r$genes,g)
  assertthat::is.string(rownames(r$binZ)[1])
  n1<-plyr::laply(rownames(r$binZ),function(x) sum(b[r$genes.dist.q==x]))
  rand.scores<-t(r$binZ)%*%n1
  if(sum(b)==1){
    raw.scores<-r$zscores[b,]
  }else{
    raw.scores<-colSums(r$zscores[b,])
  }
  scores<-(raw.scores-rand.scores)/sum(b)
  return(scores)
}

#' Matrix centering
#'
#' Code to center your matrix
#'
#' @param m matrix object
#' @param dim which dimension, rows = 1, columns = 2
#' @param sd.flag boolean for standardizing data to sd = 1 (default = F)
#' @return a matrix of zscores
#'
center.matrix<-function(m,dim = 1,sd.flag = F){
  if(dim == 1){
    zscores<-sweep(m,1,rowMeans(m,na.rm = T),FUN = '-')
  }else{
    zscores<-sweep(m,2,colMeans(m,na.rm = T),FUN = '-')
  }
  if(sd.flag){
    zscores<-sweep(zscores,dim,apply(m,dim,function(x) (sd(x,na.rm = T))),FUN = '/')
  }
  return(zscores)
}

#' Discretize data into Quantiles
#'
#' discretizes a vector of values
#'
#' @param v vector of values
#' @param n.cat number of bins
#' @param q1 quantiles (not used anymore?)
#' @return discretized data
#'
discretize.prvt<-function(v,
                          n.cat,
                          q1) {
  q1<-quantile(v,seq(from = (1/n.cat),to = 1,by = (1/n.cat)),na.rm = T)
  u<-matrix(data = 1,nrow = length(v))
  for(i in 2:n.cat){
    u[(v>=q1[i-1])&(v<=q1[i])]<-i
  }
  u<-paste0("Q",u)
  return(u)
}

#' Apply Hierarchical Linear Model with specific formula
#'
#' Wrapper to compute differentially expressed genes with
#' a hierarchical linear model
#'
#' @param r list object of sequencing data
#' @param X fixed effects
#' @param Y outcome
#' @param formula formula for the hierarchical model
#' @param ttest.flag (default = F) is flag for computing ttest?
#' @return pvalues of hlm
#'
apply.formula.HLM<-function(r,X,Y,MARGIN = 1,formula = "y ~ (1 | samples) + x",ttest.flag = F){
  if(is.matrix(Y)){
    if(ttest.flag){
      m1<-t.test.mat(Y,X)
      b<-rowSums(p.adjust.mat(m1[,1:2])<0.1,na.rm = T)>0
      m1<-m1[b,];Y<-Y[b,]
    }
    m<-t(apply(Y,MARGIN = MARGIN,function(y){formula.HLM(y,X,r,formula = formula)}))
  }else{
    m<-t(apply(X,MARGIN = MARGIN,function(x){formula.HLM(Y,x,r,formula = formula)}))
  }
  colnames(m)<-c("Estimate","P")
  m<-cbind.data.frame(Z = get.cor.zscores(m[,"Estimate"],m[,"P"]),m)
  if(ttest.flag){
    m<-cbind.data.frame(m,ttest = m1)
  }
  return(m)
}

#' Execute lmer using formula
#'
#' Wrapper to compute with hierarchical
#' linear model
#'
#' @param y list object of sequencing data
#' @param x fixed effects
#' @param r0 outcome
#' @param formula formula for the hierarchical model
#' @param ttest.flag (default = F) is flag for computing ttest?
#' @return pvalues of hlm
#' @export
formula.HLM<-function(y,x,r0, formula = "y ~ (1 | samples) + x",val = ifelse(is.numeric(x),"","TRUE"),return.all = F){
  r0$x<-x;r0$y<-y
  f<-function(r0){
    M1 <- with(r0, lmer (formula = formula))
    if(return.all){
      c1<-summary(M1)$coef[,c("Estimate","Pr(>|t|)")]
    }else{
      c1<-summary(M1)$coef[paste0("x",val),]
      idx<-match(c("Estimate","Pr(>|t|)"),names(c1))
      c1<-c1[idx]
    }
    return(c1)
  }
  c1<-tryCatch({f(r0)},
               error = function(err){return(c(NA,NA))})
  return(c1)
}

#' Subset a list object based on cell ids
#'
#' Use the cell ids to make a new list object
#' with cells subsets.
#'
#' @param r original list object
#' @param subcells ids of cells to subset
#' @return a new subsetted list object
#' @export
subset_list <- function(r, subcells) {
  n_cells <- length(r$cells)
  n_genes <- length(r$genes)
  q <- lapply(r, function(x){
    if (is.null(dim(x))) {
      # this item is a list
      if(length(x) == n_cells) {
        return (x[r$cells %in% subcells])
      } else {
        return(x)
      }
    } else if (dim(x)[2] == n_cells){
      return(x[,r$cells %in% subcells])
    } else if (dim(x)[1] == n_cells) {
      return(x[r$cells %in% subcells,])
    } else {
      return (x)
    }
  })
  return(q)
}

<<<<<<< HEAD

#' Scale and Center
#'
#' scale and center a matrix
#'
#' @param X a numeric matrix-like object
#' @param MARGIN margin
#' @return a scaled and centered matrix
#' @export
scale_and_center <- function(X, MARGIN){
  X_norm = t(apply(X, MARGIN, function(x){
    loc = mean(x, na.rm =T)
    center = x - loc
    scale = center/sd(x)
    return(scale)
  }))
  return(X_norm)
}


#' Cap Data
#'
#' Cap the values of data (matrix, data.frame, or list)
#'
#' @param X a list-list object, matrix, dataframe, or list
#' @param q quantile to cap at (e.g. 0.05 is capping at the bottom 5th and top 95th quantiles)
#' @return a capped version of your input object.
#' @export
cap_object <- function(X, q){
  ceil_q <- 1 - q
  ceil <- quantile(X, ceil_q)
  floor <- quantile(X, q)
  X[X > ceil] <- ceil
  X[X < floor] <- floor
  return(X)
}



=======
>>>>>>> 8f781af8a6374298c6baaf805667401e5c9a92fb
#' Transfer meta data from a list object to a
#' Seurat object.
#'
#' User can specify which fields to transfer from
#' a list object to a seurat object
#'
#' @param r original list object
#' @param so seurat object
#' @param transfer_list list of meta fields to transfer
#' @return a new subsetted list object
#' @export
transfer_data_list_to_so <- function(r, so, transfer_list = c("coor",
                                                              "samples",
                                                              "TMAs",
                                                              "patients",
                                                              "samples",
                                                              "sites",
                                                              "treatment",
                                                              "cell.types",
                                                              "cell.types2")){
  for (field in transfer_list) {
    if (is.null(dim(r[[field]]))) {
      so@meta.data[[field]] <- r[[field]]
    } else {
      so@meta.data <- cbind(so@meta.data, r[[field]])
    }
  }
  return(so)
}
<<<<<<< HEAD

#' Compute Spearman Correlation on two list objects
#'
#' Wrapper around correlation.
#'
#' @param v1 matrix, rows are samples, columns are genes
#' @param v2 metric to evaluate the correlation with (must be same length as rows of v1)
#' @param method (default = "spearman")
#' @param use (default = "pairwise.complete.obs") type of correlation
#' @param match.flag (default = F)
#' @param alternative (default = "two.sided") alternative hypothesis for correlation test
#' @param upper.tri.flag (default = F)
#' @return a new subsetted list object
#' @export
spearman.cor<-function(v1,v2 = NULL,method = 'spearman',use = "pairwise.complete.obs",
                       match.flag = F,alternative = "two.sided",upper.tri.flag = F){
  if(is.null(v2)){
    v2<-v1
  }
  if(!is.matrix(v1)){v1<-as.matrix(v1)}
  if(!is.matrix(v2)){v2<-as.matrix(v2)}
  if(match.flag){
    n=ncol(v1)
    if(is.null(colnames(v1))){colnames(v1)<-1:ncol(v1)}
    results<-get.mat(m.cols = c("R","P"),m.rows = colnames(v1))
    for(i in 1:ncol(v1)){
      c.i <- cor.test(v1[,i],v2[,i],method = method,use = use, alternative = alternative)
      results[i,1] <- c.i$estimate
      results[i,2] <- c.i$p.value
    }
  }else{
    n1=ncol(v1)
    m<-matrix(nrow = n1,ncol = ncol(v2))
    rownames(m)<-colnames(v1)
    colnames(m)<-colnames(v2)
    results<-list(cor = m, p = m)
    for(i in 1:n1){
      f<-function(x){
        c.i<-cor.test(v1[,i],x,method = method,use = use, alternative = alternative);
        c(c.i$estimate,c.i$p.value)}
      c.i <- apply(v2,2,f)
      results$cor[i,] <- c.i[1,]
      results$p[i,] <- c.i[2,]
    }
    if(ncol(v2)==1){
      results<-cbind(results$cor,results$p)
      colnames(results)<-c('R','P')
    }
  }
  if(upper.tri.flag){
    results$up <- cbind(results$cor[upper.tri(results$cor)],
                        results$p[upper.tri(results$p)])
  }
  return(results)
}

#' Fetch genes with top correlation
#'
#' Wrapper around a rank fetch
#'
#' @param m matrix
#' @param q (default = 100)
#' @param min.ci (default = 0)
#' @param idx (default = NULL)
#' @param add.prefix = (default = "") string to prepend to names
#' @export
get.top.cor<-function(m,q = 100,min.ci = 0,idx = NULL, add.prefix =""){
  m<-as.matrix(m)
  if(is.null(colnames(m))){colnames(m)<-1:ncol(m)}
  m.pos<-(-m);m.neg<-m

  colnames(m.pos)<-paste0(colnames(m.pos),".up")
  colnames(m.neg)<-paste0(colnames(m.neg),".down")
  v<-get.top.elements(cbind(m.pos,m.neg),q,min.ci = (-abs(min.ci)))
  names(v)<-c(colnames(m.pos),colnames(m.neg))
  if(!is.null(idx)){
    v<-v[paste(idx,c("up","down"),sep = ".")]
  }
  names(v)<-paste0(add.prefix,names(v))
  return(v)
}

#' Fetch top elements of a matrix
#'
#' Wrapper around a rank fetch
#'
#' @param m matrix
#' @param q (default = 100)
#' @param min.ci (default = 0)
#' @param idx (default = NULL)
#' @param add.prefix = (default = "") string to prepend to names
#' @export
get.top.elements<-function (m,q = 100,min.ci = NULL,main = ""){
  top.l<-list()
  v<-rownames(m)
  for (i in 1:ncol(m)){
    mi<-m[,i];mi<-mi[!is.na(mi)]
    idx<-order(mi,decreasing = F)
    ci <- mi[idx[min(q,length(mi))]]
    ci <- min(ci,min.ci)
    b <- m[,i]<=ci
    b[is.na(m[,i])]<-F
    top.l[[i]]<-sort(v[b])
  }
  if(main!=""){main<-paste0(main,".")}
  names(top.l)<-paste0(main,colnames(m))
  return(top.l)
}
=======
>>>>>>> 8f781af8a6374298c6baaf805667401e5c9a92fb
