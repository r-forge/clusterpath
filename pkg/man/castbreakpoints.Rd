\name{castbreakpoints}
\alias{castbreakpoints}
\title{castbreakpoints}
\description{Cast the variables in a clusterpath l1 data frame into the columns
of a new data frame, to facilitate plotting in [23]d. In practice
this is quite a bit faster than brute-force plotting all points
returned by predict.clusterpath.}
\usage{
castbreakpoints(df)
}
\arguments{
  \item{df}{clusterpath data frame.}
}

\value{Data frame with a row for every breakpoint in the dimensionality
of the problem, suitable for plotting with lines in the high
dimension to represent the entire linear path of optimal solutions
for the l1 problem.}

\author{Toby Dylan Hocking <toby.hocking@inria.fr>}




\examples{
sim <- gendata(N=10,D=2,K=2)
df <- clustermat.l1.id(sim$mat)
system.time(str(castbreakpoints(df)))
system.time(str(predict(df)))
system.time(str(castbreakpoints(predict(df))))  
}
