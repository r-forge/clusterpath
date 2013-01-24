\name{cluster3d}
\alias{cluster3d}
\title{cluster3d}
\description{Cluster a matrix of 3d points using the l1 and l2 algorithms.}
\usage{cluster3d(x, gamma = 0)}
\arguments{
  \item{x}{A numeric matrix of 3d data.}
  \item{gamma}{gamma value for the weights in the l2 clustering}
}

\value{An invisble list of the results of the clustering algorithms.}

\author{Toby Dylan Hocking}




\examples{
sim <- gendata(D=3,K=5,N=10,SD=0.1)
L <- cluster3d(sim$mat)
L2 <- cluster3d(sim$mat,gamma=3)

iris3 <- iris[,1:3]
cluster3d(iris3)
cluster3d(iris3,gamma=1)
}
