\name{get.best.eigenvectors}
\alias{get.best.eigenvectors}
\title{get best eigenvectors}
\description{this implements the protocol in Ng Jordan 2001 NIPS for 2
clusters. weights are calculated for a range of gamma
values. Kmeans is fit to the top 2 eigenvectors in each case, and
we pick the eigenvectors the best reconstruction error}
\usage{get.best.eigenvectors(S, K, evec.transform)}
\arguments{
  \item{S}{
}
  \item{K}{
}
  \item{evec.transform}{
}
}



\author{Toby Dylan Hocking}




\examples{
set.seed(1)
sim <- gendata(N=20,D=2,K=2,SD=0.1)
pts <- data.frame(sim$mat,class=sim$class)
xyplot(X1~X2,pts,groups=class,aspect="iso")+
  layer_(ltext(x,y,1:length(x)))

best <- get.best.eigenvectors(sim$mat,2,e.usual.k(2))
vecs <- data.frame(as.matrix(best),class=sim$class)
xyplot(X1~X2,vecs,groups=class,aspect="iso")

best <- get.best.eigenvectors(sim$mat,2,e.exponential)
vecs <- data.frame(as.matrix(best),class=sim$class)
xyplot(X1~X2,vecs,groups=class,aspect="iso")
}
