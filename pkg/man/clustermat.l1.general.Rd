\name{clustermat.l1.general}
\alias{clustermat.l1.general}
\title{clustermat l1 general}
\description{Use the l2 descent problem with general weights to solve the
weighted l1 problem, separately for each dimension. If the
multicore package is installed each dimension will be processed on
a separate core. This is useful since the path algorithm can't be
used for general weights, since it only is proven to find the
optimal solutions for the identity weights.}
\usage{
clustermat.l1.general(x, join.thresh = NULL, verbose = 0, gamma = 0, ...)
}
\arguments{
  \item{x}{Data matrix to cluster.}
  \item{join.thresh}{Threshold for fusion. If NULL take a small fraction of the nonzero
differences found in the original points x. Passed to
join.clusters2.general.}
  \item{verbose}{Passed to join.clusters2.general.}
  \item{gamma}{Used for calculating weights using \eqn{\exp(-\gamma ||x_i-x_j||^2)}.}
  \item{\dots}{passed to join.clusters2.general.}
}

\value{Data frame of optimal solutions, with 1 row for each \eqn{\alpha_i}
for each lambda.}
\references{Hocking, Joulin, Bach, Vert ICML 2011.}
\author{Toby Dylan Hocking <toby.hocking@inria.fr>}




\examples{
set.seed(7)
library(clusterpath)
sim <- gendata(D=2)
## compare with path algorithm
path <- clustermat.l1.id(sim$mat)
bpts <- castbreakpoints(path)
p <- plot2d(bpts)
descent.pts <- clustermat.l1.general(sim$mat)
p+geom_point(aes(alpha.1,alpha.2,size=lambda/max(lambda)),data=descent.pts)
desc <- melt(descent.pts[,1:4],measure=c("alpha.1","alpha.2"))
if(require(latticeExtra)){
  plot(path)+
    xyplot(value~lambda|variable,desc,groups=row)
}
library(ggplot2)
## compare with cvx
gamma <- 0.1
descent.weights <- clustermat.l1.general(sim$mat,gamma=gamma,lambda=0.001)
lvals <- c(seq(min(descent.weights$lambda),0.02,l=4),
           seq(0.02,max(descent.weights$lambda),l=4))
cvx <- cvxmod.cluster(sim$mat,lvals,norm=1,gamma=gamma)
both <- rbind(descent.weights,cvx)
molt <- melt(both,measure=c("alpha.1","alpha.2"))
ggplot(molt,aes(log(lambda+1),value,colour=solver))+
  geom_point(aes(shape=solver))+
  facet_grid(variable~.,scales="free")
}
