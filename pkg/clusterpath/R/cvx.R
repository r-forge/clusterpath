python.command <- function
### Return the command to use to run python.
(path=system.file(package="clusterpath")
### If NULL then add nothing to the PYTHONPATH, otherwise add the
### specified path.
 ){
  cmd <- "python"
  if(!is.null(path)){
    stopifnot(is.character(path))
    stopifnot(length(path)==1)
    cmd <- sprintf("PYTHONPATH=%s %s",path,cmd)
  }
  cmd
### A command that should be able to run python with cvxmod on your
### system. e.g. "PYTHONPATH=/path/to/clusterpath python"
}

cvxmod.available <- function
### Test if cvxmod is working on this system.
(python=python.command()
### How should python be invoked?
 ){
  stopifnot(is.character(python))
  stopifnot(length(python)==1)
  cmd <- sprintf("%s -c 'import cvxmod'",python)
  
  system(cmd,ignore.stderr=TRUE)==0
### TRUE if cvxmod is available, FALSE otherwise.
}


cvxmod.cluster <- structure(function
### Perform relaxed convex clustering using cxvmod. This will probably
### only work for small (N<20) matrices! This calls a python interface
### to cvxmod.
(pts,
### matrix of data, observations in rows, variables in columns.
 lambda=NULL,
### tuning parameter
 norm="2",
### use L1 or L2 or Linf norm?
 gamma=1,
### how to calculate weights between positions?
 datafile=tempfile(),
### Backing data file to use (will be read by python program).
 verbose=FALSE,
### show output from cluster.py command-line interface to cvxmod?
 s=NULL,
### 0<=s<=1 is the "natural" regularization parameter for the
### constrained version of the problem. s==0 means no differences
### between any points (the solution will be the mean) and s==1 means
### maximal differences (the solution will be the input data
### matrix). Technically you could use s>= 1 but this will yield the
### same solution as s==1.
 weight.pts=pts,
### other space of points to use for weight calculation.
 W=exp(-as.numeric(as.character(gamma))*dist(weight.pts)^2),
### distance matrix or dist object that represents weights between the
### clusters. If this is specified then we will write this to
### datafile.W, otherwise we will calculate weights based on pts and
### gamma.
 regularization.points=8
 ){
  if(!cvxmod.available()){
    stop("cvxmod not available, try installing cvxopt")
  }
  if(is.null(lambda)){
    param.vals <- if(is.null(s)){
      seq(0,1,l=regularization.points)
    }else s
    param <- "s"
  }else{
    if(is.null(s)){
      param <- "lambda"
      param.vals <- lambda
    }else{
      stop("specify lambda or s, not both")
    }
  }
  if((!missing(W)) && missing(gamma)){
    warning("W specified and gamma unspecified; setting to NA")
    gamma <- NA
  }
  W <- as.matrix(W)
  ##print(W)
  Wfile <- sprintf("%s.W",datafile)
  write.table(W,Wfile,quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table(pts,datafile,quote=FALSE,row.names=FALSE,col.names=FALSE)

  execdir <- system.file("exec",package="clusterpath")
  NAMES <- alphacolnames(pts)
  if(nchar(execdir)==0)stop("clusterpath/exec dir not found")
  clusterpy <- file.path(execdir,"cluster.py")
  cluster <- function(param.val){
    outf <- paste(datafile,param.val,sep=".")
    cmd <- sprintf("%s %s %s %s=%f %s %s",
                   python.command(),
                   shQuote(clusterpy),
                   shQuote(datafile),
                   param,as.numeric(as.character(param.val)),
                   as.character(norm),
                   shQuote(outf))
    if(verbose)cat(cmd,"\n")
    status=system(cmd,ignore.stdout=!verbose)
    if(status!=0)stop(sprintf("error code %d for command\n%s",status,cmd))
    outdata <- scan(outf,quiet=TRUE)
    outdata <- matrix(outdata,nrow=nrow(pts),ncol=ncol(pts))
    if(verbose)print(outdata)
    d <- data.frame(outdata,lambda=param.val,row=factor(1:nrow(pts)),
               gamma=factor(gamma),norm=factor(norm),solver=factor("cvxmod"))
    names(d)[names(d)=="lambda"] <- param
    names(d)[1:length(NAMES)] <- NAMES
    d
  }
  LAP <- if(require(parallel))mclapply else lapply
  ##LAP <- lapply #debug
  dfs <- LAP(param.vals,cluster)
  isdf <- sapply(dfs,is.data.frame)
  if(any(!isdf))warning(dfs[!isdf])
  results <- do.call(rbind,dfs[isdf])
  structure(results,data=pts,class=c("cvxmod","clusterpath","data.frame"),
            w=W,weight.pts=weight.pts,alphacolnames=NAMES)
### data frame of solutions from the cvxmod solver for each lambda/s
### value requested, except those for which the python program failed
### with an error.
},ex=function(){
  if(cvxmod.available()){
    ## by default if you don't specify any regularization parameters, we
    ## use some evenly spaced points in the s-parametrization.
    set.seed(16)
    sim <- gendata(N=5,D=2,K=2,SD=0.5)
    cvx <- cvxmod.cluster(sim$mat,norm=2,gamma=0.5,verbose=TRUE)
    ggplot(cvx,aes(alpha.1,alpha.2))+
      geom_path(aes(group=row),colour="grey")+
      geom_point(aes(size=s),alpha=1/2)+
      ggtitle("cvxmod solutions for the l2 problem using decreasing weights")+
      coord_equal()

    cvx <- data.frame()
    for(norm in c(1,2,"inf"))for(gamma in c(0,0.5)){
      cvx <- rbind(cvx,cvxmod.cluster(sim$mat,norm=norm,gamma=gamma))
    }
    means <- data.frame(t(colMeans(sim$mat)))
    require(grid)
    p <- ggplot(cvx,aes(alpha.2,alpha.1))+
      geom_point(aes(size=s),colour="grey")+
      facet_grid(norm~gamma,labeller=function(var,val)
                 sprintf("%s : %s",var,val))+
      coord_equal()+
      ggtitle("Fused lasso clustering for several norms and weights")+
      geom_point(aes(X2,X1),data=data.frame(sim$mat),pch=21,fill="white")+
      theme(panel.margin=unit(0,"cm"))
    print(p)
    ## otherwise, this is useful for comparing using lambda values
    set.seed(1)
    sim2 <- gendata(N=5,D=1,K=2)
    path <- clusterpath.l1.id(sim2$mat)
    lvals <- seq(0,max(path$lambda),l=8)
    cvx2 <- cvxmod.cluster(sim2$mat,norm=1,gamma=0,lambda=lvals)
    library(latticeExtra)
    plot(path)+xyplot(alpha.1~lambda,cvx2,group=row)
  }
})

cvxcheck <- structure(function
### based on a clustering result, verify using cvxmod
(df,
### data frame of l1 or l2 solutions
 lambda=sort(unique(df$lambda)),
### lambda values on which we will calculate the solutions
 ...
### passed to cvxmod.cluster
 ){
  x <- attr(df,"data")
  arglist <- list(x,lambda=lambda,norm=df$norm[1],gamma=df$gamma[1],...)
  if(!is.null(w <- attr(df,"w",exact=TRUE)))arglist$W <- w
  if(!is.null(w <- attr(df,"weight.pts",exact=TRUE)))arglist$weight.pts <- w
  do.call(cvxmod.cluster,arglist)
},ex=function(){
  sim <- gendata(N <- 5,2,2,0.1)
  colnames(sim$mat) <- c("height","length")
  xyplot(length~height,data.frame(sim$mat,row=1:N),aspect="iso",group=row)
  df <- clusterpath.l1.id(sim$mat)
  if(cvxmod.available()){
    cvx <- cvxcheck(df)
    library(reshape2)
    cvx.melt <- melt(cvx,measure.vars=1:2)

    ## plot each dimension separately using lattice
    library(latticeExtra)
    (p <- plot(df))
    update(p,main="the path algorithm (lines) agrees with cvxmod (points)")+
      xyplot(value~lambda|variable,cvx.melt,groups=row)

    ## plot the 2 dimensions together using ggplot2
    (p <- plot2d(df))
    ## compare with cvx manually
    p+
      geom_point(aes(size=lambda/max(lambda)),data=cvx,shape=17,colour="red")+
      ggtitle(paste("Optimal solutions from path algorithm (black circles)",
                  "agree with cvxmod (red triangles)"))
    ## or use a legend
    p+
      aes(shape=solver,colour=solver)+
      geom_point(aes(size=lambda/max(lambda)),data=cvx)
  }
})

