\name{plotpathmatrix}
\alias{plotpathmatrix}
\title{plotpathmatrix}
\description{plotting method for l2 path for several dimensions.}
\usage{
plotpathmatrix(path, ...)
}
\arguments{
  \item{path}{
}
  \item{\dots}{
}
}



\author{Toby Dylan Hocking <toby.hocking@inria.fr>}




\examples{
## try to cluster the iris data
path <- clusterpath.l2.general(iris[,1:4],check.splits=0,gamma=1)
plot(path,groups=iris$Species)
}
