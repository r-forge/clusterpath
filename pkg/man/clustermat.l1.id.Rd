\name{clustermat.l1.id}
\alias{clustermat.l1.id}
\title{clustermat l1 id}
\description{Cluster a matrix using the identity weights on each dimension. The
L1 problem is separable, so we can process each dimension
separately on each core if the multicore package is available.}
\usage{
clustermat.l1.id(x)
}
\arguments{
  \item{x}{Matrix of data.}
}

\value{data frame of optimal solutions at the breakpoints. need unique()
when there are multiple lines that join at the exact same time
(only pathological cases).}

\author{Toby Dylan Hocking <toby.hocking@inria.fr>}




\examples{
x <- c(-3,-2,0,3,5)
df <- clustermat.l1.id(x)
head(df)
mean(x)
plot(df)
## check agreement with cvx
cres <- cvxcheck(df,seq(0,max(df$lambda),l=8),verbose=TRUE)
orig <- data.frame(alpha=x,row=1:length(x),lambda=0)
p <- ggplot(df,aes(lambda,alpha))+
  geom_line(aes(group=row))+
  geom_point(aes(y=alpha.1),data=cres,pch=21)+
  scale_y_continuous(breaks=x,minor=min(x):max(x))
print(p)
}
