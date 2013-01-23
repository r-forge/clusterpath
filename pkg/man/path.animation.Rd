\name{path.animation}
\alias{path.animation}
\title{path animation}
\description{Create an animation that summarizes the path calculation in
2d. This works best for the path algorithm, since each animation
frame will represent a breakpoint in the path.}
\usage{path.animation(df, outdir = tempfile(), title = "Calculation of the path in 2d")}
\arguments{
  \item{df}{Result of clusterpath.l1.id.}
  \item{outdir}{Subdirectory for plot files, to be created.}
  \item{title}{Title for the animation.}
}



\author{Toby Dylan Hocking}




\examples{
## several clusters
sim <- gendata(5,2,2,0.1)
colnames(sim$mat) <- c("height","length")
plot(sim$mat)
df <- clusterpath.l1.id(sim$mat)
path.animation(df)

iris2 <- iris[,c("Petal.Length","Petal.Width")]
iris.path <- clusterpath.l1.id(iris2)
path.animation(iris.path)
}
