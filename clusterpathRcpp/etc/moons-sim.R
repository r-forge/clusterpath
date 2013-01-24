sim.two.moons <- function(seed,N=150,noise=0.5){
  set.seed(seed)
  halfcircle <- function(r,center=c(0,0),class,sign){
    angle <- runif(N,0,pi)
    rad <- rnorm(N,r,noise)
    data.frame(x=rad*cos(angle)+center[1],y=sign*rad*sin(angle)+center[2],
               class=factor(class))
  }
  pts <- rbind(halfcircle(4,c(0,0),1,1),
               halfcircle(4,c(4,2),2,-1))
  data.frame(pts,seed=factor(seed))
}
alldata <- do.call(rbind,lapply(1:20,sim.two.moons))
library(lattice)
xyplot(y~x|seed,alldata,groups=class,aspect="iso")
write.csv(alldata,"moons-sim-clusterpath.csv",row.names=FALSE,quote=FALSE)
