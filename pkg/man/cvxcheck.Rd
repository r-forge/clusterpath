\name{cvxcheck}
\alias{cvxcheck}
\title{cvxcheck}
\description{based on a clustering result, verify using cvxmod}
\usage{
cvxcheck(df, lambda = sort(unique(df$lambda)), ...)
}
\arguments{
  \item{df}{data frame of l1 or l2 solutions}
  \item{lambda}{lambda values on which we will calculate the solutions}
  \item{\dots}{passed to cvxmod.cluster}
}



\author{Toby Dylan Hocking <toby.hocking@inria.fr>}




\examples{
sim <- gendata(N <- 5,2,2,0.1)
colnames(sim$mat) <- c("height","length")
xyplot(length~height,data.frame(sim$mat,row=1:N),aspect="iso",group=row)
df <- clusterpath.l1.id(sim$mat)
cvx <- cvxcheck(df)
library(reshape)
cvx.melt <- melt(cvx,measure.vars=1:2)

## plot each dimension separately using lattice
library(latticeExtra)
(p <- plot(df))
update(p,main="the path algorithm (lines) agrees with cvxmod (points)")+
  xyplot(value~lambda|variable,cvx.melt,groups=row)

## plot the 2 dimensions together using ggplot2
(p <- plot2d(df))
## compare with cvx manually
p+
  geom_point(aes(size=lambda/max(lambda)),data=cvx,shape=17,colour="red")+
    opts(title=paste("Optimal solutions from path algorithm (black circles)",
           "agree with cvxmod (red triangles)"))
## or use a legend
p+
  aes(shape=solver,colour=solver)+
    geom_point(aes(size=lambda/max(lambda)),data=cvx)
}
