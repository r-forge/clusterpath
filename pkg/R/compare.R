### This list defines a certain number of clustering algorithms which
### we can compare. Each element is a list with 2 elements: [[1]] is
### the clustering method name, [[2]] is a function which takes a data
### matrix m and number of clusters k and then returns a vector of
### integer cluster assignments.
clusterings <-
  list(list("old spectral clusterpath",function(m,k,verbose=0,...){
    V <- get.old.best(m,k,exponential.eigenvectors)
    path <- join.clusters2.general(V,verbose=verbose,check.split=0,lambda=1e-6,
                                   join.thresh=1e-3)
    cut.clusterpath.labels(path,k)
  }),list("old spectral exp kmeans",function(m,k,verbose=0,...){
    V <- get.old.best(m,k,exponential.eigenvectors)
    kmeans(V,k)$cluster
  }),list("old spectral kmeans",function(m,k,...){
    V <- get.old.best(m,k,usual.eigenvectors)
    kmeans(V,k)$cluster
  }),list("spectral exp kmeans manual",function(m,k,...){
    N <- nrow(m)
    A <- Matrix(as.matrix(exp(-10*dist(m)^2)))
    sums <- apply(A,1,sum)
    Dm12 <- Diagonal(N,1/sqrt(sums))
    L <- Dm12 %*% A %*% Dm12
    Id <- Diagonal(N,1)
    Lnorm <- Id-L
    edec <- eigen(Lnorm,symmetric=TRUE)
    eii <- exp(-10000*edec$values)
    nonzero <- eii>1e-8
    E <- Diagonal(sum(nonzero),eii[nonzero])
    V <- Matrix(edec$vectors[,nonzero])
    X <- V %*% E
    print(dim(X))
    kmeans(X,k)$cluster
  }),list("clusterpath",function(m,k,path=NULL,gamma=2,verbose=0,...){
    library(clusterpath)
    if(is.null(path))
      {##option to provide your own clusterpath tree for speed
        path <- join.clusters2.general(m,check.splits=0,verbose=verbose,
                                       join.thresh=0.05,opt.thresh=0.1,
                                       lambda.factor=1.05,gamma=gamma,
                                       target.cluster=k,
                                       ...)
      }
    cut.clusterpath.labels(path,k)
  }),list("spectral kmeans",function(m,k,...){
    V <- get.best.eigenvectors(m,k,e.usual.k(k))
    kmeans(V,k)$cluster
  }),list("spectral exp kmeans auto",function(m,k,...){
    V <- get.best.eigenvectors(m,k,e.exponential.auto)
    kmeans(V,k)$cluster
  }),list("average linkage",function(m,k,...){
    tree <- hclust(dist(m),"average")
    cutree(tree,k)
  }),list("kmeans",function(m,k,...){
    kmeans(m,k,100)$cluster
  }),list("gaussian mixture",function(m,k,...){
    require(mclust)
    Mclust(m,k)$classification
  }))

### Call this function with a matrix m and number of clusters k to run
### several clustering algorithms on the matrix.
cluster.points <- structure(function(m,k,...){
  colnames(m) <- alphacolnames(m)
  call.do <- function(args,what)do.call(what,args)
  ## return results of all methods
  cluster.and.time <- function(method,FUN){
    print(method)
    seconds <- system.time({
      guess <- FUN(as.matrix(m),k,...)
    })["elapsed"]
    data.frame(m,row=seq_along(guess),method,guess,seconds,row.names=NULL)
  }
  d <- do.call(rbind,lapply(clusterings,call.do,cluster.and.time))
  ##annotate for plotting
  attr(d,"data") <- m
  attr(d,"alphacolnames") <- alphacolnames(m)
  d
### A data frame with timings and guesses
},ex=function(){
  library(clusterpath)
  iriSc <- scale(as.matrix(iris[,1:4]))
  iclust <- cluster.points(iriSc,3,verbose=1)
  splom(~iclust[1:4]|method,iclust,groups=guess)
  table(iclust$guess,rep(iris$Species,nlevels(iclust$method)),iclust$method)
})

### measure of correspondence between partitions (Hubert and Arabie
### 1985) 1=perfect, 0=completely random, or just the same label for
### everybody.
norm.rand.index <- function(klass,guess){
  n <- table(klass,guess)
  ch2 <- function(x)x*(x-1)/2
  sumi <- sum(ch2(rowSums(n)))
  sumj <- sum(ch2(colSums(n)))
  expected <- sumi*sumj/ch2(sum(n))
  numerator <- sum(ch2(n))-expected
  denominator <- (sumi+sumj)/2-expected
  numerator/denominator
}
