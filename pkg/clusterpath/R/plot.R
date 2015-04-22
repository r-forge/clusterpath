### plotting method for l2 path for several dimensions.
plotpathmatrix <- structure(function(path,...){
  orig <- attr(path,"data")
  last <- path$lambda==max(path$lambda)
  splom(~orig,panel=function(i,j,...){
    panel.splom(i=i,j=j,...)
    ##i==1 y axis, j==2 x axis first call
    x <- path[,j]
    y <- path[,i]
    clusters <- unique(cbind(x,y)[last,])
    panel.xyplot(x,y,col="grey",type="l",groups=path[,"row"],
                 subscripts=1:nrow(path))
    panel.xyplot(clusters[,1],clusters[,2],col="black",pch=20)
  },...)
},ex=function(){
  ## try to cluster the iris data
  path <- clusterpath.l2.general(iris[,1:4],check.splits=0,gamma=1)
  plot(path,groups=iris$Species)
})

### Get the names of the columns as they will be shown in the
### clustering results data frames.
alphacolnames <- function(x){
  if(is.null(colnames(x)))
     sprintf("alpha.%d",1:ncol(x))
  else colnames(x)
}

### plot the path in rgl 3d as lines. uses only first 3 variables by default
plotpath3d <- structure(function(path,vars=1:3,...){
  for(i in unique(path$row)){
    plot3d(subset(path,row==i)[,vars],add=TRUE,type="l",...)
  }
},ex=function(){
  ## calculate the l1 clusterpath
  breakpoints <- clusterpath.l1.id(iris[,1:4])
  plot(breakpoints)
  ## calculate the weighted l2 clusterpath
  path <- clusterpath.l2(iris[,1:4],gamma=1)
  plot(path,groups=iris$Species)
  require(rgl)
  plot3d(iris[,1:3],type="p",box=FALSE,aspect=TRUE)
  plotpath3d(path,col="blue")
  bpts4d <- castbreakpoints(breakpoints)
  plotpath3d(bpts4d,col="orange")
})

cluster3d <- structure(function
### Cluster a matrix of 3d points using the l1 and l2 algorithms.
(x,
### A numeric matrix of 3d data.
 gamma=0
### gamma value for the weights in the l2 clustering
 ){
  require(rgl)
  stopifnot(ncol(x)==3)
  solutions <- clusterpath.l2(x,gamma=gamma)
  ##3d plotting code
  plot3d(x,type="p",box=FALSE,aspect=TRUE)
  plotpath3d(solutions,col="blue",lwd=2)
  ## add l1 path
  path <- clusterpath.l1.id(x)
  bpts <- castbreakpoints(path)
  plotpath3d(bpts,col="orange",lwd=2)
  invisible(list(l1.breakpoints.dims=path,
                 l1.breakpoints.3d=bpts,
                 l2.solutions=solutions))
### An invisble list of the results of the clustering algorithms.
},ex=function(){
  sim <- gendata(D=3,K=5,N=10,SD=0.1)
  L <- cluster3d(sim$mat)
  L2 <- cluster3d(sim$mat,gamma=3)

  iris3 <- iris[,1:3]
  cluster3d(iris3)
  cluster3d(iris3,gamma=1)
})

plot.clusterpath <- structure(function
### Plot method for clusterpath data frames. These data frames contain
### all the breakpoints for each dimension, so we can plot lines to
### represent exactly the entire path of optimal solutions.
(x,
### data frame returned by clusterpath.l1.id
 type="l",
### plot type, see ?xyplot
 main="The entire regularization path of optimal solutions for each variable",
### main title
 xlab=expression(paste("location in the regularization path  ",lambda)),
### annotation for the horizontal axis
 ylab=expression(paste("optimal coefficient  ",alpha)),
### annotation for the vertical axis
 strip=strip.custom(strip.names=TRUE),
### strip function for annotation in the strips
 layout=c(1,nlevels(x$col)),
### layout for the plot, by default the panels are stacked vertically
 ...
### Other arguments for xyplot
 ){
  xyplot(alpha~lambda|col,x,group=row,type=type,layout=layout,
         xlab=xlab,ylab=ylab,main=main,strip=strip,...)
### lattice xyplot showing the entire path of solutions.
},ex=function(){
  sim <- gendata(N=20,D=5,K=5)
  df <- clusterpath.l1.id(sim$mat)
  plot(df)

  df.iris <- clusterpath.l1.id(as.matrix(iris[,1:4]))
  plot(df.iris)
})

### plot method for l2 solution path
plot.l2 <- function(x,...){
if(length(attr(x,"alphacolnames"))==2)plot2d(x)
else plotpathmatrix(x,...)
}

plot2d <- structure(function
### Plot the result of clusterpath for a 2d data set. In this case,
### the breakpoints in the clusterpath data frame are not sufficient
### to visualize the correct path of solutions in 2d. You need
### additional breakpoints, which are calculated by castbreakpoints
### before plotting using ggplot2.
(x
### clusterpath data frame as returned by clusterpath.l1.id.
 ){
  if(!is.null(x$col)){
    if(nlevels(x$col)==2){
      d <- castbreakpoints(x)
      colnames <- levels(x$col)
    }else{
      stop(sprintf("dont know how to plot %d variables in 2d",length(colnames)))
    }
  }else{
    if(is.null(colnames <- attr(x,"alphacolnames"))){
      stop("dont know how to plot this thing")
    }else{
      d <- x
    }
  }
  a <- eval(as.call(lapply(c("aes",colnames),as.name)))
  ggplot(d,a)+
    geom_path(aes(group=row),colour="grey")+
    coord_equal()+
    geom_point(aes(size=lambda/max(lambda)),pch=21)+
    scale_size(expression(lambda/lambda[max]))+
    ggtitle("Path of optimal solutions for data in 2d")
},ex=function(){
  sim <- gendata(N=5,D=2,K=2)
  df <- clusterpath.l1.id(sim$mat)
  plot2d(df)
})

fillin <- function
### Use linear interpolation to complete the data frame of optimal
### alpha values in a given column. For the linear interpolation, we
### use the \code{\link{approx}} function, with x-value lambda and
### y-value from col.
(alpha,
### data frame containing columns lambda and col
 col
### Column to complete
 ){
  y <- alpha[,col]
  na <- is.na(y)
  x <- alpha$lambda
  y[na] <- approx(x[!na],y[!na],x[na])$y
  y[is.na(y)] <- rev(y[!na])[1] ## grab the last one
  alpha[,col] <- y
  alpha
### Same as original alpha, but with NAs in column col replaced by
### their linear interpolation.
}

castbreakpoints <- structure(function
### Cast the variables in a clusterpath l1 data frame into the columns
### of a new data frame, to facilitate plotting in [23]d. In practice
### this is quite a bit faster than brute-force plotting all points
### returned by predict.clusterpath.
(df
### clusterpath data frame.
 ){
  NAMES <- levels(df$col)
  d <- dcast(rename(df,c(alpha="value")),...~col,mean)
  for(col in NAMES)d <- ddply(d,.(row),fillin,col)
  structure(d[,c(NAMES,names(d)[!names(d)%in%NAMES])],
            data=attr(df,"data"),
            alphacolnames=NAMES)
### Data frame with a row for every breakpoint in the dimensionality
### of the problem, suitable for plotting with lines in the high
### dimension to represent the entire linear path of optimal solutions
### for the l1 problem.
},ex=function(){
  sim <- gendata(N=10,D=2,K=2)
  df <- clusterpath.l1.id(sim$mat)
  system.time(str(castbreakpoints(df)))
  system.time(str(predict(df)))
  system.time(str(castbreakpoints(predict(df))))  
})

### Plot the result of a predictclusterpath data frame as points. 
plot.predictclusterpath <- function(x,type="p",...){
  plot.clusterpath(x,type=type,...)
}

