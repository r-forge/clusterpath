\name{fillin}
\alias{fillin}
\title{fillin}
\description{Use linear interpolation to complete the data frame of optimal
alpha values in a given column. For the linear interpolation, we
use the \code{\link{approx}} function, with x-value lambda and
y-value from col.}
\usage{fillin(alpha, col)}
\arguments{
  \item{alpha}{data frame containing columns lambda and col}
  \item{col}{Column to complete}
}

\value{Same as original alpha, but with NAs in column col replaced by
their linear interpolation.}

\author{Toby Dylan Hocking}





