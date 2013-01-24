### calculate some statistics at every point in the path that are
### useful for guessing how many clusters there are.
calcstats <- function(path){
  NAMES <- attr(path,"alphacolnames")
  X <- attr(path,"data")
  w <- attr(path,"wvec")
  do.call(rbind,lapply(unique(path$lambda),function(lambda){
    alpha <- path[path$lambda==lambda,NAMES,drop=FALSE]
    penalty <- sum(dist(alpha)*w)
    loss <- sum((X-alpha)^2)/2
    k <- nrow(unique(alpha))
    data.frame(lambda,k,penalty,loss)
  }))
}

### Given a clusterpath data frame, cut it at a certain number of
### clusters and return the tree up to that point
cut.clusterpath <- function(x,k,...){
  stats <- calcstats(x)
  kstats <- stats[stats$k==k,]
  if(nrow(kstats)==0)stop(sprintf("didnt find %d clusters!",k))
  lcut <- subset(kstats,lambda==min(lambda))$lambda
  x[x$lambda<=lcut,]
}

### given a matrix of optimal alpha values, return a vector of labels
assign.labels <- function(alpha){
  ## alpha starts in good order
  centers <- data.frame(unique(alpha))
  centers$guess <- 1:nrow(centers)
  alpha$row <- 1:nrow(alpha)
  guess <- merge(centers,as.data.frame(alpha))
  guess <- guess[order(guess$row),]
  guess[,"guess"]
}

### given a clusterpath tree, cut it and give the labels
cut.clusterpath.labels <- function(x,k,...){
  tree <- cut.clusterpath(x,k)
  alpha <- subset(tree,lambda==max(lambda))[,attr(x,"alphacolnames")]
  assign.labels(alpha)
}

