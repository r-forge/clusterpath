\name{join.clusters2}
\alias{join.clusters2}
\title{join clusters2}
\description{Cluster a matrix of data using the L2 penalty and fused group
lasso.}
\usage{
join.clusters2(x, ...)
}
\arguments{
  \item{x}{Matrix of data to cluster.}
  \item{\dots}{passed to join.clusters2.general}
}



\author{Toby Dylan Hocking <toby.hocking@inria.fr>}




\examples{
## cluster the iris data
path <- join.clusters2.general(scale(iris[,1:4]),gamma=1)
plot(path,groups=iris$Species)
}
