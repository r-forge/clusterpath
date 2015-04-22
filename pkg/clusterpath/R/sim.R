### generate a multivariate normal mixture
gendata <- function(N=5,K=2,D=3,SD=0.1){
  means <- matrix(rnorm(K*D),K,D)##K rows (classes), D columns(dimensions)
  mat <- apply(means,2,sapply,function(class.mean)rnorm(N,class.mean,SD))
  list(mat=mat,means=means,class=rep(1:K,each=N))
}

### simulate some data in the shape of 2 half-moons and run some
### clustering algorithms on them.
halfmoon.cluster <- function(seed,N=150,noise=0.5,...){
  set.seed(seed)
  halfcircle <- function(r,center=c(0,0),class,sign){
    angle <- runif(N,0,pi)
    rad <- rnorm(N,r,noise)
    data.frame(x=rad*cos(angle)+center[1],y=sign*rad*sin(angle)+center[2],
               class=factor(class))
  }
  pts <- rbind(halfcircle(4,c(0,0),1,1),
               halfcircle(4,c(4,2),2,-1))
  ##print(xyplot(y~x,pts,groups=class,aspect="iso"))
  ##browser()
  m <- pts[,c("x","y")]
  path <- NULL
  guess <- cluster.points(m,2,path=path,...)
  rows <- guess$method=="clusterpath"
  ##guess[rows,"seconds"] <- guess[rows,"seconds"]+add.seconds
  p <- xyplot(y~x|method,guess,aspect="iso",group=guess,
              main=paste("Clusterpath beats hierarchical and k-means,",
                "and performs similarly to spectral for non-convex clusters"))
  print(p)
  ##stats <- calcstats(path)
  ##l <- min(subset(stats,k==2)$lambda)
  list(data.frame(guess,seed,class=pts$class))##,subset(path,lambda<=l))
}
