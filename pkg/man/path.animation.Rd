\name{path.animation}
\alias{path.animation}
\title{path animation}
\description{Create an animation that summarizes the path calculation in
2d. This works best for the path algorithm, since each animation
frame will represent a breakpoint in the path.}
\usage{
path.animation(df, outdir = tempfile(), title = "Calculation of the path in 2d")
}
\arguments{
  \item{df}{Result of clustermat.l1.id.}
  \item{outdir}{Subdirectory for plot files, to be created.}
  \item{title}{Title for the animation.}
}



\author{Toby Dylan Hocking <toby.hocking@inria.fr>}




\examples{
## several clusters
sim <- gendata(N <- 5,2,2,0.1)
colnames(sim$mat) <- c("height","length")
xyplot(length~height,data.frame(sim$mat,row=1:N),aspect="iso",group=row)
df <- clustermat.l1.id(sim$mat)
##path.animation(df)
}
