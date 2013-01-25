\name{cvxmod.cluster}
\alias{cvxmod.cluster}
\title{cvxmod cluster}
\description{Perform relaxed convex clustering using cxvmod. This will probably
only work for small (N<20) matrices! This calls a python interface
to cvxmod.}
\usage{cvxmod.cluster(pts, lambda = NULL, norm = "2", gamma = 1, datafile = tempfile(), 
    verbose = FALSE, s = NULL, weight.pts = pts, W = exp(-as.numeric(as.character(gamma)) * 
        dist(weight.pts)^2), regularization.points = 8)}
\arguments{
  \item{pts}{matrix of data, observations in rows, variables in columns.}
  \item{lambda}{tuning parameter}
  \item{norm}{use L1 or L2 or Linf norm?}
  \item{gamma}{how to calculate weights between positions?}
  \item{datafile}{Backing data file to use (will be read by python program).}
  \item{verbose}{show output from cluster.py command-line interface to cvxmod?}
  \item{s}{0<=s<=1 is the "natural" regularization parameter for the
constrained version of the problem. s==0 means no differences
between any points (the solution will be the mean) and s==1 means
maximal differences (the solution will be the input data
matrix). Technically you could use s>= 1 but this will yield the
same solution as s==1.}
  \item{weight.pts}{other space of points to use for weight calculation.}
  \item{W}{distance matrix or dist object that represents weights between the
clusters. If this is specified then we will write this to
datafile.W, otherwise we will calculate weights based on pts and
gamma.}
  \item{regularization.points}{
}
}

\value{data frame of solutions from the cvxmod solver for each lambda/s
value requested, except those for which the python program failed
with an error.}

\author{Toby Dylan Hocking}




\examples{
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
               sprintf("\%s : \%s",var,val))+
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
}
