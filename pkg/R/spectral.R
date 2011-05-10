### scale the eigenvectors using eigenvalues and the exponential
### function
exponential.eigenvectors <- function
(nonzero.vectors,nonzero.values,bandwidth=1){
  w.exp <- exp(-bandwidth*((nonzero.values/min(nonzero.values))-1))
  w.pos <- w.exp>1e-8
  nonzero.vectors[,w.pos]%*%diag(w.exp[w.pos])
}

### just take the first 2 eigenvectors
usual.eigenvectors <- function(vectors,...){
  vectors[,1:2]
}

### try a bunch of gamma values in spectral clustering and get the one
### that gives the best reconstruction error.
get.old.best <- function(m,k,evec.transform){
  lists <- lapply(c(1:10),function(gamma){
  ##lists <- lapply(c(10^((-7):2)),function(gamma){
    require(Matrix)
    w <- as.matrix(exp(-gamma*dist(m)^2))
    Gm12 <- Diagonal(nrow(m),1/sqrt(apply(w,1,sum)))
    G <- Diagonal(nrow(m),apply(w,1,sum))
    L <- G - w
    Id <- Diagonal(nrow(m),1)
    Lnorm <- Id - Gm12 %*% w %*% Gm12
    norm.dec <- eigen(Lnorm,symmetric=TRUE)
    edec <- eigen(L,symmetric=TRUE)
    chosen <- rev(which(edec$values>1e-8))
    nonzero.vectors <- edec$vectors[,chosen]
    nonzero.values <- edec$values[chosen]
    V <- cbind(1)
    bandwidth <- 10
    while(ncol(V)==1){
      bandwidth <- bandwidth/10
      V <- evec.transform(nonzero.vectors,nonzero.values,bandwidth)
    }
    rownorm.vecs <-
      V/matrix(sqrt(rowSums(V^2)),nrow(V),ncol(V))
    ##browser()
    guess <- kmeans(rownorm.vecs,k,nstart=25)$cluster
    class.error <- sapply(unique(guess),function(g){
      class.pts <- rownorm.vecs[guess==g,]
      sum(scale(class.pts,scale=FALSE)^2)
    })
    list(eigenvectors=rownorm.vecs,stats=data.frame(
         error=sum(class.error),
         gamma=gamma,
         vectors=ncol(rownorm.vecs),
         bandwidth=bandwidth,
         method="old"))
         
  })
  stats <- do.call(rbind,lapply(lists,"[[","stats"))
  ##print(stats)
  chosen <- lists[[which.min(stats$error)]]
  chosen$eigenvectors
}  

### adaptive eigenvector scaling function that takes the biggest
### bandwidth that gives you more than 1 nonzero vector.
e.exponential.auto <- function
(eigenvalues
 ){
  bandwidth <- 1e12
  e <- 1
  while(sum(e>1e-8)==1){
    e <- exp(-bandwidth*eigenvalues)
    bandwidth <- bandwidth/10
  }
  e
}

### eigenvector scaling function with fixed bandwidth.
e.exponential <- function(evals,bw=1e9){
  e <- exp(-bw*evals)
  print(cbind(evals,e))
  e
}

### eigenvector scaling function which just take the first k
### eigenvalues
e.usual <- function
(eigenvalues,
 k=stop("must specify k")
 ){
  increasing <- sort(eigenvalues)
  thresh <- increasing[k]
  ifelse(eigenvalues<=thresh,1,0)
}

### give an eigenvector scaling function for a specific k
e.usual.k <- function(K){
  function(v)e.usual(v,K)
}

### this implements the protocol in Ng Jordan 2001 NIPS for 2
### clusters. weights are calculated for a range of gamma
### values. Kmeans is fit to the top 2 eigenvectors in each case, and
### we pick the eigenvectors the best reconstruction error
get.best.eigenvectors <- structure(function(S,K,evec.transform){
  get.vectors <- function(FUN,edec){
    eii <- FUN(edec$values)
    ##print(cbind(edec$values,eii))
    nonzero <- eii>1e-8
    E <- Diagonal(sum(nonzero),eii[nonzero])
    V <- Matrix(edec$vectors[,nonzero])
    V %*% E
  }
  try.gamma <- function(gamma,FUN){
    require(Matrix)
    N <- nrow(S)
    A <- Matrix(as.matrix(exp(-gamma*dist(S)^2)))
    sums <- apply(A,1,sum)
    if(any(sums==0))return(list(error=NA))
    Dm12 <- Diagonal(N,1/sqrt(sums))
    L <- Dm12 %*% A %*% Dm12
    Id <- Diagonal(N,1)
    Lnorm <- Id-L
    edec <- eigen(Lnorm,symmetric=TRUE)
    X <- get.vectors(FUN,edec)
    ## row normalization:
    eigenvector.l2norms <- sqrt(apply(X^2,1,sum))
    if(any(eigenvector.l2norms==0))return(list(error=NA))
    Y <- as.matrix(X/matrix(eigenvector.l2norms,nrow(X),ncol(X)))
    if(K==2){
      means <- matrix(0,K,ncol(Y))
      means[1,] <- Y[1,]
      dots <- Y %*% Y[1,]
      means[2,] <- Y[which.min(dots),]
      guess <- kmeans(Y,means,nstart=25)$cluster
    }else{
      guess <- kmeans(Y,K,nstart=25)$cluster
    }
    class.error <- sapply(unique(guess),function(g){
      Yk <- Y[guess==g,]
      sum(scale(Yk,scale=FALSE)^2)
    })
    list(eigenvectors=Y,stats=data.frame(
         error=sum(class.error),
         gamma=gamma,
         vectors=ncol(Y),
         method="new"),
         eigen.decomposition=edec)
  }
  lists <- list()
  gamma <- 1e-8
  ## this is a function of the eigenvalues that just returns 1 for the
  ## smallest k values, and 0 elsewhere
  while(!any(is.na(sapply(lists,"[[","error")))){
    lists[[length(lists)+1]] <- try.gamma(gamma,e.usual.k(K))
    gamma <- gamma*10
  }
  stats <- do.call(rbind,lapply(lists[-length(lists)],"[[","stats"))
  print(stats)
  r <- which.min(stats$error)
  gamma.range <- stats[c(r-1,r+1),"gamma"]
  gamma <- seq(gamma.range[1],gamma.range[2],l=20)
  lists <- lapply(gamma,try.gamma,e.usual.k(K))
  stats2 <- do.call(rbind,lapply(lists,"[[","stats"))
  print(stats2)
  gamma <- stats2[which.min(stats2$error),"gamma"]
  chosen <- try.gamma(gamma,evec.transform)
  print(chosen$stats)
  chosen$eigenvectors
},ex=function(){
  set.seed(1)
  sim <- gendata(N=20,D=2,K=2,SD=0.1)
  pts <- data.frame(sim$mat,class=sim$class)
  xyplot(X1~X2,pts,groups=class,aspect="iso")+
    layer_(ltext(x,y,1:length(x)))

  best <- get.best.eigenvectors(sim$mat,2,e.usual.k(2))
  vecs <- data.frame(as.matrix(best),class=sim$class)
  xyplot(X1~X2,vecs,groups=class,aspect="iso")

  best <- get.best.eigenvectors(sim$mat,2,e.exponential)
  vecs <- data.frame(as.matrix(best),class=sim$class)
  xyplot(X1~X2,vecs,groups=class,aspect="iso")
})

