\name{clusterpath.l2.general}
\alias{clusterpath.l2.general}
\title{clusterpath l2 general}
\description{Cluster a matrix of data using the L2 penalty and weights that
potentially depend on another matrix.}
\usage{clusterpath.l2.general(x, w = x, gamma = 1, lambda = 0.01, join.thresh = NULL, 
    opt.thresh = 0.001, lambda.factor = 1.05, smooth = 0, maxit = 1000, 
    linesearch.freq = 2, linesearch.points = 10, check.splits = 1, 
    target.cluster = 0, verbose = 0, ...)}
\arguments{
  \item{x}{\eqn{n\times p} matrix of data to cluster}
  \item{w}{\eqn{w_{ij}} weight matrix to use for the coefficients of the
n*(n-1)/2 pairwise fusion penalties, specified as a "dist"
object. Use w=as.dist(M) if you have a \eqn{n\times n} matrix
M. Otherwise, this can be a \eqn{n\times k} matrix Y of data, from
which the \eqn{w_{ij}} will be calculated using \eqn{w_{ij} =
\exp(-\gamma||Y_i-Y_j||^2)}, where \eqn{Y_i} is row \eqn{i} of \eqn{Y}.}
  \item{gamma}{>= 0. Degree of dependence of the weights on the differences
between initial points. 0 means identity weights.}
  \item{lambda}{starting regularization parameter.}
  \item{join.thresh}{Threshold on \eqn{||\alpha_i-\alpha_j||_2} for fusing points \eqn{i} and
\eqn{j}. NULL means take a small percentage of the smallest nonzero
difference in the original data matrix.}
  \item{opt.thresh}{Threshold on the gradient for deciding when we have found an
optimal solution.}
  \item{lambda.factor}{Factor by which we increase lambda after an optimal solution is
found for the current lambda.}
  \item{smooth}{Smoothing parameter \eqn{\epsilon} for calculating the gradient:
\eqn{||\alpha_i-\alpha_j||_{2,\epsilon} =
\sqrt{\epsilon + \sum_{k=1}^p \alpha_{ik}-\alpha_{jk}}}}
  \item{maxit}{number of gradient steps to take before abandoning the
optimization.}
  \item{linesearch.freq}{How often to do a line search? linesearch.freq=0 corresponds to
never doing a line search (always doing decreasing step size), and
linesearch.freq=k means do a line search k times for every
gradient step. Look at the examples to see that linesearch speeds
up the optimization.}
  \item{linesearch.points}{On how many points should we calculate the line search?}
  \item{check.splits}{0=do not unfuse clusters after finding an optimal solution for
this lambda (faster). 1=unfuse clusters (more accurate, allows for
paths with cluster splits).}
  \item{target.cluster}{0=calculate the entire path. Otherwise, calculate this number of
clusters and return.}
  \item{verbose}{0=no printout, 1=report every optimal solution found. 2=report
every gradient and line search step as well (very slow).}
  \item{\dots}{ignored}
}

\value{Data frame of optimal solutions found, one row for every
\eqn{\alpha_i} for every lambda. First p columns are the \eqn{alpha_{ik}}
for k=1,...,p with names taken from the original matrix, if there
were any. Then rows lambda, row, gamma, norm, solver which permit
plotting and comparing with other norms, weights and solvers.}

\author{Toby Dylan Hocking}




\examples{
library(clusterpath)
x <- cbind(c(1.1,1.2,4,2),
           c(2.2,2.5,3,0))
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
pts <- data.frame(alpha=x,i=1:nrow(x))
toplot <- rbind(res2,cvx2)
ggplot(,aes(alpha.1,alpha.2))+
  geom_path(aes(group=row),data=subset(toplot,solver=="descent.split"),
            colour="grey")+
  geom_point(aes(size=lambda,colour=k,shape=solver),data=toplot)+
  scale_shape_manual(values=c(21,20))+
  geom_text(aes(label=i),data=pts)+
  ggtitle("Descent with split checks agrees with cvxmod")+
  coord_equal()
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
  REP <- sprintf("\%03d +- \%03d ms \\1",ms.mean,ms.sd)
  levels(res$solver) <- gsub("descent.([a-z]+)",REP,levels(res$solver))
  attr(res,"ms.mean") <- ms.mean
  attr(res,"ms.sd") <- ms.sd
  res
}
nosplit.res <- with.timelabels(check.splits=0)
split.res <- with.timelabels(check.splits=1)
lvals <- split.res$lambda
cvx <- cvxcheck(split.res,lambda=seq(min(lvals),max(lvals),l=6))
p <- ggplot(rbind(split.res,nosplit.res),aes(alpha.1,alpha.2,colour=row))+
  geom_path(aes(group=interaction(row,solver),linetype=solver),
            lwd=1.5)+
  ggtitle(paste("Descent solver (lines) needs to permit splitting",
         "to reach the cvxmod solutions (circles) for general weights"))+
  geom_point(aes(size=lambda),data=cvx,pch=1)+
  scale_colour_identity()+
  geom_text(aes(label=i,colour=i),data=pts)+
  coord_equal()
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
}
