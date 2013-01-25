\name{cluster.points}
\alias{cluster.points}
\title{cluster points}
\description{Call this function with a matrix m and number of clusters k to run
several clustering algorithms on the matrix.}
\usage{cluster.points(m, k, ...)}
\arguments{
  \item{m}{
}
  \item{k}{
}
  \item{\dots}{
}
}

\value{A data frame with timings and guesses}

\author{Toby Dylan Hocking}




\examples{
iriSc <- scale(as.matrix(iris[,1:4]))
iclust <- cluster.points(iriSc,3,verbose=1)
splom(~iclust[1:4]|method,iclust,groups=guess)
table(iclust$guess,rep(iris$Species,nlevels(iclust$method)),iclust$method)
}
