\name{plot.clusterpath}
\alias{plot.clusterpath}
\title{plot clusterpath}
\description{Plot method for clusterpath data frames. These data frames contain
all the breakpoints for each dimension, so we can plot lines to
represent exactly the entire path of optimal solutions.}
\usage{\method{plot}{clusterpath}(x, type = "l", main = "The entire regularization path of optimal solutions for each variable", 
    xlab = expression(paste("location in the regularization path  ", 
        lambda)), ylab = expression(paste("optimal coefficient  ", 
        alpha)), strip = strip.custom(strip.names = TRUE), layout = c(1, 
        nlevels(x$col)), ...)}
\arguments{
  \item{x}{data frame returned by clusterpath.l1.id}
  \item{type}{plot type, see ?xyplot}
  \item{main}{main title}
  \item{xlab}{annotation for the horizontal axis}
  \item{ylab}{annotation for the vertical axis}
  \item{strip}{strip function for annotation in the strips}
  \item{layout}{layout for the plot, by default the panels are stacked vertically}
  \item{\dots}{Other arguments for xyplot}
}

\value{lattice xyplot showing the entire path of solutions.}

\author{Toby Dylan Hocking}




\examples{
sim <- gendata(N=20,D=5,K=5)
df <- clusterpath.l1.id(sim$mat)
plot(df)

df.iris <- clusterpath.l1.id(as.matrix(iris[,1:4]))
plot(df.iris)
}
