\name{clusterpath.l2}
\alias{clusterpath.l2}
\title{clusterpath l2}
\description{Cluster a matrix of data using the L2 penalty and fused group
lasso.}
\usage{clusterpath.l2(x, ...)}
\arguments{
  \item{x}{Matrix of data to cluster.}
  \item{\dots}{passed to clusterpath.l2.general}
}



\author{Toby Dylan Hocking}




\examples{
## cluster the iris data
path <- clusterpath.l2(iris[,1:4],gamma=1)
plot(path,groups=iris$Species)
## compare with l1 path
bpts <- clusterpath.l1.id(iris[,1:4])
plot(bpts)
bpts4d <- castbreakpoints(bpts)
solutions <- rbind(bpts4d,path)
ggplot(solutions,aes(Petal.Length,Petal.Width))+
  geom_path(aes(group=interaction(row,solver),colour=norm))
}
