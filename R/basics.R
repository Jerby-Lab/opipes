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
#'
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
#'
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

#' Compute _t_-tests
#'
#' computes _t_-tests for a matrix
#'
#' @param m markers matrix (subsampled)
#' @param b boolean for the right cluster
#' @param two.sided boolean to determine whether you want to run a two-sided t-test
#' @param rankf (default = F)
#' @param fold.changeF (default = F)
#' @return p-value table
#'
t.test.mat<-function(m,b,two.sided=F,rankf = F,fold.changeF = F){
  # m is the markers matrix (subsampled)
  # b is the boolean for the cluster!
  if(length(b)!=ncol(m)){
    print("Error. Inconsistent no. of samples.")
    return()
  }
  if(sum(b)<2||sum(!b)<2){
    return(get.mat(rownames(m),c('more','less',"zscores")))
  }
  if(two.sided){
    p<-as.matrix(apply(m,1,function(x) t.test(x[b],x[!b])$p.value))
  }else{
    p<-t(apply(m,1,function(x) c(t.test(x[b],x[!b],alternative = 'greater')$p.value,
                                 t.test(x[b],x[!b],alternative = 'less')$p.value)))
    colnames(p)<-c('more','less')
    p<-cbind(p,get.p.zscores(p))
    colnames(p)[3]<-"zscores"
  }
  if(rankf){
    p<-cbind(p,rank(p[,1]),rank(p[,2]))
    colnames(p)[4:5]<-c("rank.more","rank.less")
  }
  if(fold.changeF){
    p<-cbind.data.frame(p,pos.mean = rowMeans(m[,b]),neg.mean = rowMeans(m[,!b]))
    p$logFC<-log2(p$pos.mean/p$neg.mean)
  }

  return(p)
}

