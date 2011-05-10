\name{predict.clusterpath}
\alias{predict.clusterpath}
\title{predict clusterpath}
\description{based on calculated solution path breakpoints, find solutions at
an individual point lambda.}
\usage{
\method{predict}{clusterpath}(object, lambda = unique(object$lambda), ...)
}
\arguments{
  \item{object}{clusterpath data frame of breakpoints returned by clustermat.l1.id}
  \item{lambda}{lambda values for which we will calculate the optimal solutions.}
  \item{\dots}{ignored.}
}

\value{data frame subclass with the same columns as object, but
containing all the solutions calculated for lambda.}

\author{Toby Dylan Hocking <toby.hocking@inria.fr>}





