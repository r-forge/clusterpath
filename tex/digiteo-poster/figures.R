library(tikzDevice)
options( tikzDocumentDeclaration="\\documentclass[24pt]{beamer}")
library(clusterpath)
library(ggplot2)
if(require(RColorBrewer)){
  require(grDevices)
  custom.pal <- brewer.pal(9,"Set1")
  custom.pal[6] <- "#DDDD33"
  trellis.par.set(theme=simpleTheme(col=custom.pal))
  ##show.settings()
}

centers <- expand.grid(x=c(1,6),y=c(0,5,10))
centers$cluster <- 1:nrow(centers)
N <- 100
set.seed(9)
library(plyr)
pts <- ddply(centers,.(cluster),with,data.frame(
      x=rnorm(N*length(x),x),
      y=rnorm(N*length(y),y),
      cluster=cluster))
pts$kmeans <- kmeans(as.matrix(pts[,1:2]),nrow(centers))$cluster
library(lattice)
p <- xyplot(x~y,pts,groups=kmeans,aspect="iso",main="k-means clustering, k=6",
            xlab="",ylab="",scales=list())
#tikz("kmeans.tex",height=5,width=8)
pdf("kmeans.pdf",height=5,width=8)
print(p)
dev.off()



### Figure 2. 10 points to demonstrate different path geometries
set.seed(19)
sim <- gendata(N=5,D=2,K=2,SD=0.6)
plot(sim$mat,asp=1)
cvx <- data.frame()
## for(norm in c(1,2,"inf"))for(gamma in c(0,1)){
##   cvx <- rbind(cvx,cvxmod.cluster(sim$mat,norm=norm,gamma=gamma,
##                                   regularization.points=50))
## }
load(url("http://cbio.ensmp.fr/~thocking/clusterpath-figure-2.RData"))
cutline <- 0.16
extra <- cvxmod.cluster(sim$mat,cutline,norm=1,gamma=0)
means <- data.frame(alpha=t(colMeans(sim$mat)))
normweights <- function(var,val){
  val <- as.character(val)
  if(var=="gamma")var <- "\\gamma"
  else var <- sprintf("\\textrm{%s}",var)
  val[val=="inf"] <- "\\infty"
  sprintf("$%s=%s$",var,val)
}
p <- ggplot(cvx,aes(alpha.1,alpha.2))+
  geom_text(data=means,label="$\\overline X$",col="grey")+
  ##geom_point(aes(size=s),alpha=1/4)+
  geom_path(aes(group=row),lwd=1)+
  geom_point(col="red",data=extra,size=3)+
  facet_grid(norm~gamma,labeller=normweights)+
  coord_equal()+
  scale_size("$s$")+
  geom_point(data=data.frame(alpha=sim$mat),pch=21,fill="white")+
  scale_x_continuous("$\\alpha^1$")+
  scale_y_continuous("$\\alpha^2$")+
  opts(title="Norm and weights control the clusterpath")+
  theme_bw()
print(p)
tikz("cvx-allnorms.tex",height=8,width=8.5)
print(p+theme_bw())
dev.off()

pts <- data.frame(alpha=sim$mat)
pts$row <- 1:nrow(pts)
path <- clusterpath.l1.id(pts[,1:2])
path2d <- castbreakpoints(path)
pred <- subset(castbreakpoints(predict(path,c(0,cutline))),lambda>0)
toplot <- path
levels(toplot$col) <- gsub("alpha.([12])","$\\\\alpha^\\1$",levels(path$col))
path.plot <- ggplot(toplot,aes(lambda,alpha))+
  geom_line(aes(group=row))+
  geom_point(subset=.(lambda==0),pch=21,fill="white")+
  geom_vline(xintercept=cutline,colour="red",lwd=2)+
  facet_grid(col~.,scales="free")+
  xlab("Location in the regularization path $\\lambda$")+
  ylab("Optimal value of $\\ell_1$ clusterpath")+
  opts(axis.title.y=theme_text(vjust=-0.1,angle=90),
       title="$O(pn\\log n)$ path algorithm for $\\ell_1$ clusterpath")
tikz("l1-dims.tex",height=8,width=6)
print(path.plot)
dev.off()



### Figure 1. geometric interpretation
set.seed(3)
x <- replicate(2,rnorm(3))
pts <- data.frame(alpha=x,row=1:nrow(x))
getlines <- function(x){
  require(foreach)
  N <- nrow(x)
  foreach(i=1:(N-1),.combine=rbind)%do%{
    foreach(j=(i+1):N,.combine=rbind)%do%{
      start <- x[i,]
      end <- x[j,]
      data.frame(start.1=start[1],end.1=end[1],start.2=start[2],end.2=end[2],
                 norm=factor(2))
    }
  }
}
getpolys <- function(path,l=0){
  gamma <- as.numeric(as.character(path$gamma[1]))
  W <- as.matrix(exp(-gamma*dist(attr(path,"weight.pts"))))
  alpha <- subset(path,lambda==l)[,c("row",attr(path,"alphacolnames"))]
  clusters <- unique(alpha[,-1])
  N <- nrow(clusters)
  equal <- function(u)apply(alpha[,-1],1,function(v)all(v==u))
  plot(as.matrix(clusters),type="n",asp=1,xlim=c(-2,1),ylim=c(-2,1))
  with(alpha,text(alpha.1,alpha.2,row))
  polys <- data.frame()
  ws <- data.frame()
  width <- data.frame()
  for(i in 1:(N-1))for(j in (i+1):N){
    ci <- as.matrix(clusters[i,])
    cj <- as.matrix(clusters[j,])
    m <- (ci[2]-cj[2])/(ci[1]-cj[1])
    rects <- expand.grid(i=alpha[equal(ci),"row"],j=alpha[equal(cj),"row"])
    rects$wij <- W[as.matrix(rects)]
    total <- sum(rects$wij)
    ## this is the horizontal displacement from ci
    intercept <- ci[2]-ci[1]*m
    abline(intercept,m)
    int2 <- ci[1]/m+ci[2]
    m2 <- -1/m
    iperp <- function(x)m2*(x-ci[1])+ci[2]
    jperp <- function(x)m2*(x-cj[1])+cj[2]
    abline(int2,m2)
    xfound <- ci[1]-sqrt(((total/2)^2)/(m2^2+1))
    xoff <- abs(ci[1]-xfound)
    xvals <- xoff*c(-1,1)+ci[1]
    points(xvals,iperp(xvals),pch="+")
    ## this is the angle of the rectangle from the horizontal
    wsum <- c(0,cumsum(rects$wij))
    rate <- xoff/total*2
    ix <- function(w)w*rate+ci[1]-xoff
    jx <- function(w)w*rate+cj[1]-xoff
    for(k in 1:(length(wsum)-1)){
      IX <- c(ix(wsum[k]),ix(wsum[k+1]))
      JX <- c(jx(wsum[k+1]),jx(wsum[k]))
      m <- cbind(c(IX,JX),c(iperp(IX),jperp(JX)))
      row <- t(c(start=colMeans(m[c(1,4),]),end=colMeans(m[c(2,3),])))
      width <- rbind(width,data.frame(row,rect.num=k))
      ij <- with(rects,sprintf("%d%d",i[k],j[k]))
      ws <- rbind(ws,data.frame(alpha=t(colMeans(m)),
                                label=sprintf("$w_{%s}$",ij)))
      polys <- rbind(polys,data.frame(alpha=m,rect.num=factor(k),
                    pair=ij,row.names=NULL))
    }
  }
  list(polys=polys,clusters=clusters,ws=ws,width=width)
}
idname <- "Identity weights, $t=\\Omega(X)$"
line.df <- data.frame(getlines(x),figure=factor(idname))
path <- clusterpath.l2(x,verbose=1,join.thresh=0.01,gamma=1.3)
plot(path)
lval <- path$lambda[300]
poly.list <- getpolys(path,lval)
jointitle <- sprintf("Decreasing weights after join, $t<\\Omega(X)$")
add_polygon <- function(p,L,figure){
  p+
  geom_polygon(aes(group=interaction(pair,rect.num)),
               data=data.frame(L$polys,figure),
               fill="grey",colour="white",lwd=1.2)+
  geom_segment(aes(start1,start2,xend=end1,yend=end2),
               data=data.frame(L$width,figure),
               colour="white")+
  geom_text(aes(label=label),
            data=data.frame(L$ws,figure))
}
p <- ggplot(poly.list$clusters,aes(alpha.1,alpha.2))+
  coord_equal()+
  geom_point(pch=21,fill="white")
add_polygon(p,poly.list,jointitle)
apart <- "\\small Decreasing weights, $t=\\Omega(X)$"
apart.list <- getpolys(path,path$lambda[1])
outside <- data.frame(start=rbind(x,x),figure=idname,
                      end.1=c(x[1,1],x[1,1],x[3,1],x[3,1],x[3,1],x[3,1]),
                      end.2=c(x[2,2],x[2,2],x[2,2],x[1,2],x[2,2],x[1,2]))
segs <- rbind(outside,line.df[,names(outside)])
getcenter <- function(d,ann,xoff=0,yoff=0){
  with(d,data.frame(x=(start.1+end.1)/2+xoff,y=(start.2+end.2)/2+yoff,
                    label=sprintf("$\\ell_%s$",ann),
                    figure=idname))
}
inf <- getcenter(outside,"\\infty")[c(1,4,5),]
inf$y[1] <- inf$y[1]+0.2
inf$x[-1] <- inf$x[-1]+0.2
inf$y[3] <- inf$y[3]+0.1
norms <- rbind(getcenter(line.df,"2"),getcenter(outside,"1"),inf)
norms$y[6] <- norms$y[6]+0.05
norms$y[8] <- norms$y[8]+0.1
p <- ggplot(NULL,aes(alpha.1,alpha.2))+
  geom_segment(aes(start.1,start.2,xend=end.1,yend=end.2),
               data=segs,lwd=2,colour="grey")+
  geom_text(aes(x,y,label=label),data=norms)
p2 <- add_polygon(p,apart.list,apart)+
  geom_path(aes(group=row),data=data.frame(path,figure=jointitle),
            subset=.(lambda<=lval))
p3 <- add_polygon(p2,poly.list,jointitle)+
  geom_point(data=data.frame(poly.list$clusters,figure=jointitle))+
  coord_equal()+
  geom_text(aes(label=label),
            data=data.frame(alpha.1=c(-1,x[2,1]-0.05,0.35),
              alpha.2=c(-1.25,0.35,x[3,2]),
              label=sprintf("$X_%d$",1:3),
              figure=rep(c(jointitle,idname,apart),each=3)))+
  facet_grid(~figure)+
  geom_text(aes(label=label),
            data.frame(alpha.1=c(-0.8,-0.5),alpha.2=c(-1.05,0.1),
                       figure=jointitle,
                       label=c("$\\alpha_1$",
                         "$\\alpha_{C}=\\alpha_2=\\alpha_3$")))
levs <- unique(unlist(lapply(p3$layers,function(l)levels(l$data$figure))))
allpts <- do.call(rbind,lapply(levs,function(l)data.frame(pts,figure=l)))
p4 <- p3+geom_point(data=allpts,pch=21,fill="white")+
 scale_x_continuous("",breaks=c(-10))+
 scale_y_continuous("",breaks=c(-10))+
 theme_bw()+
  opts(title="Geometric interpretation: constrain area between points")
print(p4)
tikz("geometry.tex",height=7,width=14)
print(p4)
dev.off()
