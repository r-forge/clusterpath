\name{plotpath3d}
\alias{plotpath3d}
\title{plotpath3d}
\description{plot the path in rgl 3d as lines. uses only first 3 variables by default}
\usage{
plotpath3d(path, vars = 1:3, ...)
}
\arguments{
  \item{path}{
}
  \item{vars}{
}
  \item{\dots}{
}
}



\author{Toby Dylan Hocking <toby.hocking@inria.fr>}




\examples{
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
}
