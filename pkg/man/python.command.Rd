\name{python.command}
\alias{python.command}
\title{python command}
\description{Return the command to use to run python.}
\usage{python.command(path = system.file(package = "clusterpath"))}
\arguments{
  \item{path}{If NULL then add nothing to the PYTHONPATH, otherwise add the
specified path.}
}

\value{A command that should be able to run python with cvxmod on your
system. e.g. "PYTHONPATH=/path/to/clusterpath python"}

\author{Toby Dylan Hocking}





