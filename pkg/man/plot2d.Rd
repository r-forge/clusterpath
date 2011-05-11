\name{plot2d}
\alias{plot2d}
\title{plot2d}
\description{Plot the result of clusterpath for a 2d data set. In this case,
the breakpoints in the clusterpath data frame are not sufficient
to visualize the correct path of solutions in 2d. You need
additional breakpoints, which are calculated by castbreakpoints
before plotting using ggplot2.}
\usage{
plot2d(x)
}
\arguments{
  \item{x}{clusterpath data frame as returned by clusterpath.l1.id.}
}



\author{Toby Dylan Hocking <toby.hocking@inria.fr>}




\examples{
sim <- gendata(N=5,D=2,K=2)
df <- clusterpath.l1.id(sim$mat)
plot2d(df)
}
