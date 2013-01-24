clusterpath.l2 <- structure(function
### Cluster a matrix of data using the L2 penalty and fused group
### lasso.
(x,
### Matrix of data to cluster.
 ...
### passed to clusterpath.l2.general
 ){
  clusterpath.l2.general(x=x,w=x,check.splits=0,...)
},ex=function(){
  ## cluster the iris data
  path <- clusterpath.l2(iris[,1:4],gamma=1)
  plot(path,groups=iris$Species)
  ## compare with l1 path
  bpts <- clusterpath.l1.id(iris[,1:4])
  plot(bpts)
  bpts4d <- castbreakpoints(bpts)
  solutions <- rbind(bpts4d,path)
  ggplot(solutions,aes(Petal.Length,Petal.Width))+
    geom_path(aes(group=interaction(row,solver),colour=norm))
})

clusterpath.l2.general <- structure(function
### Cluster a matrix of data using the L2 penalty and weights that
### potentially depend on another matrix.
(x,
### \eqn{n\times p} matrix of data to cluster
 w=x,
### \eqn{w_{ij}} weight matrix to use for the coefficients of the
### n*(n-1)/2 pairwise fusion penalties, specified as a "dist"
### object. Use w=as.dist(M) if you have a \eqn{n\times n} matrix
### M. Otherwise, this can be a \eqn{n\times k} matrix Y of data, from
### which the \eqn{w_{ij}} will be calculated using \eqn{w_{ij} =
### \exp(-\gamma||Y_i-Y_j||^2)}, where \eqn{Y_i} is row \eqn{i} of \eqn{Y}.
 gamma=1,
### >= 0. Degree of dependence of the weights on the differences
### between initial points. 0 means identity weights.
 lambda=0.01,
### starting regularization parameter.
 join.thresh=NULL,
### Threshold on \eqn{||\alpha_i-\alpha_j||_2} for fusing points \eqn{i} and
### \eqn{j}. NULL means take a small percentage of the smallest nonzero
### difference in the original data matrix.
 opt.thresh=0.001,
### Threshold on the gradient for deciding when we have found an
### optimal solution.
 lambda.factor=1.05,
### Factor by which we increase lambda after an optimal solution is
### found for the current lambda.
 smooth=0,
### Smoothing parameter \eqn{\epsilon} for calculating the gradient:
### \eqn{||\alpha_i-\alpha_j||_{2,\epsilon} =
### \sqrt{\epsilon + \sum_{k=1}^p \alpha_{ik}-\alpha_{jk}}}
 maxit=1000,
### number of gradient steps to take before abandoning the
### optimization.
 linesearch.freq=2,
### How often to do a line search? linesearch.freq=0 corresponds to
### never doing a line search (always doing decreasing step size), and
### linesearch.freq=k means do a line search k times for every
### gradient step. Look at the examples to see that linesearch speeds
### up the optimization.
 linesearch.points=10,
### On how many points should we calculate the line search?
 check.splits=1,
### 0=do not unfuse clusters after finding an optimal solution for
### this lambda (faster). 1=unfuse clusters (more accurate, allows for
### paths with cluster splits).
 target.cluster=0,
### 0=calculate the entire path. Otherwise, calculate this number of
### clusters and return.
 verbose=0,
### 0=no printout, 1=report every optimal solution found. 2=report
### every gradient and line search step as well (very slow).
 ...
### ignored
 ){
  x <- as.matrix(x)
  wvec <- if(class(w)=="dist"){
    if(missing(gamma))
      warning("gamma unspecified; gamma column in output is nonsense!")
    as.vector(w)
  }else{
    wmat <- exp(-gamma * dist(as.matrix(w))^2)
    as.vector(wmat)
    ##or:
    ##.Call("calcW_convert",w,gamma,PACKAGE="clusterpath")
  }
  if(!is.numeric(x) | !is.numeric(wvec))
    stop("must use numeric data for x and w")
  if(is.null(join.thresh)){ ##default join threshold depends on scale
    dx <- dist(x)
    join.thresh <- min(dx[dx>0])/10
    if(verbose>=1)cat("using join threshold ",join.thresh,"\n",sep="")
  }
  r <- .Call("join_clusters2_restart_convert",x,wvec,
             lambda,join.thresh,opt.thresh,lambda.factor,smooth,
             maxit,linesearch.freq,linesearch.points,check.splits,
             target.cluster,verbose,
             PACKAGE="clusterpath")
  r <- do.call(rbind,r)
  alpha <- r
  colnames(alpha)[1:ncol(x)] <- NAMES <- alphacolnames(x)
  alpha$solved <- alpha$grad<opt.thresh
  alpha$row <- factor(alpha$i+1)
  alpha$gamma <- factor(gamma)
  alpha$norm <- factor(2)
  alpha$solver <-
    if(check.splits)factor("descent.split") else factor("descent.no split")
  if(!is.null(rownames(x)))levels(alpha$row) <- rownames(x)
  ret <- structure(alpha[,c(NAMES,"lambda","row","gamma","norm","solver")],
                   data=x,alphacolnames=NAMES,
                   class=c("l2","clusterpath","data.frame"),
                   wvec=wvec)
  if(class(w)=="dist"){
    attr(ret,"w") <- w
  }else{
    attr(ret,"weight.pts") <- as.matrix(w)
  }
  ret
### Data frame of optimal solutions found, one row for every
### \eqn{\alpha_i} for every lambda. First p columns are the \eqn{alpha_{ik}}
### for k=1,...,p with names taken from the original matrix, if there
### were any. Then rows lambda, row, gamma, norm, solver which permit
### plotting and comparing with other norms, weights and solvers.
},ex=function(){
  library(clusterpath)
  x <- cbind(c(1.1,1.2,4,2),
             c(2.2,2.5,3,0))
  pts <- data.frame(alpha=x,i=1:nrow(x))
  Wmat <- matrix(0,4,4)
  Wmat[2,1] <- Wmat[1,2] <- 9
  Wmat[3,1] <- Wmat[1,3] <- 20
  Wmat[3,2] <- Wmat[2,3] <- 1
  Wmat[4,1] <- Wmat[1,4] <- 1
  Wmat[4,2] <- Wmat[2,4] <- 20
  Wmat[4,3] <- Wmat[3,4] <- 1
  w <- as.dist(Wmat)
  res <- clusterpath.l2.general(x,w,gamma=NA,lambda=0.001,lambda.factor=1.1)
  lvals <- unique(res$lambda)
  if(cvxmod.available()){
    cvx <- cvxcheck(res,lambda=seq(0,max(lvals),l=8))
    library(plyr)
    cvx2 <- ddply(cvx,.(lambda),function(d){
      d$k <- factor(nrow(unique(round(d[,1:2],1))))
      d
    })
    res2 <- ddply(res,.(lambda,solver),function(d){
      d$k <- factor(nrow(unique(d[,1:2])))
      d
    })
    toplot <- rbind(res2,cvx2)
    ggplot(,aes(alpha.1,alpha.2))+
      geom_path(aes(group=row),data=subset(toplot,solver=="descent.split"),
                colour="grey")+
      geom_point(aes(size=lambda,colour=k,shape=solver),data=toplot)+
      scale_shape_manual(values=c(21,20))+
      geom_text(aes(label=i),data=pts)+
      ggtitle("Descent with split checks agrees with cvxmod")+
      coord_equal()
  }
  ## descent algo splitting ex, with matrix and times
  with.timelabels <- function(...){
    ms <- c()
    for(R in 1:20){
      tt <- system.time({
        res <- clusterpath.l2.general(x,w,gamma=NA,
                                      lambda=0.001,
                                      lambda.factor=1.02,
                                      join.thresh=0.01,...)
      })
      ms <- c(ms,tt["elapsed"]*1000)
    }
    ms.mean <- round(mean(ms))
    ms.sd <- round(sd(ms))
    REP <- sprintf("%03d +- %03d ms \\1",ms.mean,ms.sd)
    levels(res$solver) <- gsub("descent.([a-z]+)",REP,levels(res$solver))
    attr(res,"ms.mean") <- ms.mean
    attr(res,"ms.sd") <- ms.sd
    res
  }
  nosplit.res <- with.timelabels(check.splits=0)
  split.res <- with.timelabels(check.splits=1)
  lvals <- split.res$lambda
  p <- ggplot(rbind(split.res,nosplit.res),aes(alpha.1,alpha.2,colour=row))+
    geom_path(aes(group=interaction(row,solver),linetype=solver),
              lwd=1.5)+
    ggtitle(paste("Descent solver (lines) needs to permit splitting",
           "to reach the optimal solutions for general weights"))+
    scale_colour_identity()+
    geom_text(aes(label=i,colour=i),data=pts)+
    coord_equal()
  if(cvxmod.available()){
    cvx <- cvxcheck(split.res,lambda=seq(min(lvals),max(lvals),l=6))
    p <- p+
      geom_point(aes(size=lambda),data=cvx,pch=1)+
      ggtitle(paste("Descent solver (lines) needs to permit splitting",
           "to reach the cvxmod solutions (circles) for general weights"))
  }
  print(p)
  ## compare decreasing step size with mixed decreasing/linesearch
  LAPPLY <- if(require(multicore))mclapply else lapply
  lsres <- LAPPLY(c(0,2,10,20),function(lsval){
    print(lsval)
    r <- with.timelabels(linesearch.freq=lsval,maxit=10e4)
    data.frame(linesearch.freq=lsval,
               ms.mean=attr(r,"ms.mean"),
               ms.sd=attr(r,"ms.sd"))
  })
  (lstab <- do.call(rbind,lsres))
})

clusterpath.l1.general <- structure(function
### Use the l2 descent problem with general weights to solve the
### weighted l1 problem, separately for each dimension. If the
### multicore package is installed each dimension will be processed on
### a separate core. This is useful since the path algorithm can't be
### used for general weights, since it only is proven to find the
### optimal solutions for the identity weights.
##references<< Hocking, Joulin, Bach, Vert ICML 2011.
(x,
### Data matrix to cluster.
 join.thresh=NULL,
### Threshold for fusion. If NULL take a small fraction of the nonzero
### differences found in the original points x. Passed to
### clusterpath.l2.general.
 verbose=0,
### Passed to clusterpath.l2.general.
 gamma=0,
### Used for calculating weights using \eqn{\exp(-\gamma ||x_i-x_j||^2)}.
 ...
### passed to clusterpath.l2.general.
 ){
  if(is.null(join.thresh)){
    dx <- dist(x)
    join.thresh <- min(dx[dx>0])/1000
    if(verbose>=1)cat("using join threshold ",join.thresh,"\n",sep="")
  }
  w <- exp(-gamma*dist(x)^2)
  LAPPLY <- if(require(multicore))mclapply else lapply
  dfs <- LAPPLY(1:ncol(x),function(j){
    clusterpath.l2.general(x[,j],w=w,join.thresh=join.thresh,gamma=gamma,
                           verbose=verbose,...)
  })
  ## need to add rows for lambda values not present in all dimensions
  lambda.lists <- lapply(dfs,function(d)unique(d$lambda))
  lvals <- sort(unique(unlist(lambda.lists)))
  mat <- do.call(cbind,lapply(dfs,function(d){
    alpha <- subset(d,lambda==max(lambda))
    results <- lapply(lvals,function(l){
      if(l%in%d$lambda)d[d$lambda==l,1] else alpha[,1]
    })
    unlist(results)
  }))
  biggest <- which.max(sapply(lambda.lists,length))
  big.df <- dfs[[biggest]]
  result <- transform(data.frame(mat,big.df[,-1]),norm=1)
  attr(result,"alphacolnames") <- names(result)[1:ncol(x)] <- alphacolnames(x)
  class(result) <- class(big.df)
  for(a in c("data","wvec","weight.pts")){
    attr(result,a) <- attr(big.df,a)
  }
  result
### Data frame of optimal solutions, with 1 row for each \eqn{\alpha_i}
### for each lambda.
},ex=function(){
  set.seed(7)
  library(clusterpath)
  sim <- gendata(D=2)
  ## compare with path algorithm
  path <- clusterpath.l1.id(sim$mat)
  bpts <- castbreakpoints(path)
  p <- plot2d(bpts)
  descent.pts <- clusterpath.l1.general(sim$mat)
  p+geom_point(aes(alpha.1,alpha.2,size=lambda/max(lambda)),data=descent.pts)
  desc <- melt(descent.pts[,1:4],measure=c("alpha.1","alpha.2"))
  if(require(latticeExtra)){
    plot(path)+
      xyplot(value~lambda|variable,desc,groups=row)
  }
  library(ggplot2)
  ## compare with cvx
  if(cvxmod.available()){
    gamma <- 0.1
    descent.weights <- clusterpath.l1.general(sim$mat,gamma=gamma,lambda=0.001)
    lvals <- c(seq(min(descent.weights$lambda),0.02,l=4),
               seq(0.02,max(descent.weights$lambda),l=4))
    cvx <- cvxmod.cluster(sim$mat,lvals,norm=1,gamma=gamma)
    both <- rbind(descent.weights,cvx)
    molt <- melt(both,measure=c("alpha.1","alpha.2"))
    ggplot(molt,aes(log(lambda+1),value,colour=solver))+
      geom_point(aes(shape=solver))+
      facet_grid(variable~.,scales="free")
  }
})

clusterpath.l1.id <- structure(function
### Cluster a matrix using the identity weights on each dimension. The
### L1 problem is separable, so we can process each dimension
### separately on each core if the multicore package is available.
(x,
### Matrix of data.
 LAPPLY=if(require(multicore))mclapply else lapply
### Function to use to combine the results of each dimension. Defaults
### to mclapply for parallel processing if the multicore package is
### available, otherwise the standard lapply.
 ){
  x <- as.matrix(x)
  stopifnot(is.numeric(x))
  dfs <- LAPPLY(1:ncol(x),function(k){
    L <- .Call("join_clusters_convert",x[,k],PACKAGE="clusterpath")
    data.frame(L[1:2],
               row=factor(L$i+1),
               col=factor(k),
               gamma=factor(0),
               norm=factor(1),
               solver=factor("path"))
  })
  df <- do.call(rbind,dfs)
  if(!is.null(rownames(x)))levels(df$row) <- rownames(x)
  levels(df$col) <- alphacolnames(x)
  d <- structure(df,data=x,class=c("l1","clusterpath","data.frame"))
  unique(d)
### data frame of optimal solutions at the breakpoints. need unique()
### when there are multiple lines that join at the exact same time
### (only pathological cases).
},ex=function(){
  x <- c(-3,-2,0,3,5)
  df <- clusterpath.l1.id(x)
  head(df)
  mean(x)
  plot(df)
  ## check agreement with cvx
  if(cvxmod.available()){
    cres <- cvxcheck(df,seq(0,max(df$lambda),l=8),verbose=TRUE)
    orig <- data.frame(alpha=x,row=1:length(x),lambda=0)
    p <- ggplot(df,aes(lambda,alpha))+
      geom_line(aes(group=row))+
      geom_point(aes(y=alpha.1),data=cres,pch=21)+
      scale_y_continuous(breaks=x,minor=min(x):max(x))
    print(p)
  }
})

predict.clusterpath <- function
### based on calculated solution path breakpoints, find solutions at
### an individual point lambda.
(object,
### clusterpath data frame of breakpoints returned by clusterpath.l1.id
 lambda=unique(object$lambda),
### lambda values for which we will calculate the optimal solutions.
 ...
### ignored.
 ){
  require(plyr)
  calc1 <- function(d){
    dm <- merge(d,data.frame(lambda),all=TRUE)
    tofill <- transform(dm,row=row[1],col=col[1],gamma=gamma[1],
                        norm=norm[1],solver=solver[1])
    fillin(tofill,"alpha")
  }
  d <- ddply(object,.(col,row),calc1)[,names(object)]
  LL <- lambda ## only return lambdas that were asked for
  d <- subset(d,lambda %in% LL)
  structure(d,data=attr(object,"data"),
            class=c("predictclusterpath",class(object)))
### data frame subclass with the same columns as object, but
### containing all the solutions calculated for lambda.
}

