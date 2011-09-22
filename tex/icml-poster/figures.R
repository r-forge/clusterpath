library(tikzDevice)
options( tikzDocumentDeclaration="\\documentclass[24pt]{beamer}")
library(clusterpath)
library(ggplot2)
theme_set(theme_grey())
set.seed(1 )
sim <- gendata(5,2,2,0.3)
pts <- data.frame(alpha=sim$mat)
pts$row <- 1:nrow(pts)
path <- clusterpath.l1.id(pts[,1:2])
path2d <- castbreakpoints(path)
cutline <- 0.18
pred <- subset(castbreakpoints(predict(path,c(0,cutline))),lambda>0)
mytheme <-   opts(axis.title.x=theme_text(vjust=-0.5),
  axis.text.x=theme_text(vjust=-0.5),
  axis.ticks.margin=unit(0.55,"lines"))
p <- ggplot(pts,aes(alpha.2,alpha.1))+
  geom_path(aes(group=row),data=path2d)+
  geom_point(subset=.(row==10))+
  geom_point(pch=21)+
  geom_point(data=pred,colour="red")+
  coord_equal()+
  opts(title="$\\ell_1$ clusterpath of 10 points in 2d")+
  geom_text(aes(x,y,label=label),data.frame(x=0.7,y=-0.6,label="Joins the left cluster on $\\alpha^1$\\ before joining right cluster."))+
  geom_text(aes(x,y,label=label),data.frame(x=0.3,y=-0.1,label="Solution at $\\lambda=0.18$\\ yields 2 clusters."),colour="red")+
  xlab("$\\alpha^2$\\hspace{2in}")+
  ylab("$\\alpha^1$")+
  mytheme
tikz("l1-2d.tex",width=10,height=6,pointsize=24)
print(p)
dev.off()


library(ggplot2)
toplot <- path
levels(toplot$col) <- gsub("alpha.([12])","$\\\\alpha^\\1$",levels(path$col))
path.plot <- ggplot(toplot,aes(lambda,alpha))+
  geom_line(aes(group=row))+
  geom_point(subset=.(lambda==0 & row==10))+
  geom_vline(xintercept=cutline,colour="red")+
  facet_grid(col~.,scales="free")+
  xlab("Location in the regularization path $\\lambda$")+
  ylab("Optimal value of $\\ell_1$ clusterpath")+
  mytheme+opts(axis.title.y=theme_text(vjust=-0.1,angle=90))
tikz("l1-dims.tex",height=6,width=4)
print(path.plot)
dev.off()
 

library(clusterpath)
library(multicore)
agreement.matrix <- function(l){
  1-dist(sapply(unique(l),function(x)l==x),"max")
}
calc.curves <- function(x,true.labels,kmax=12,gamma.vals=c(0.5,2,10)){
  names(gamma.vals) <- sprintf("clusterpath $\\gamma=%.1f$",gamma.vals)
  FUNS <- c(mclapply(gamma.vals,function(gamma){
    path <-
      clusterpath.l2(x,gamma=gamma,join.thresh=0.05,opt.thresh=0.1,
                             verbose=1,lambda.factor=1.01)
    stats <- calcstats(path)
    substats <- subset(stats,k<=kmax)
    function(kval){
      l <- subset(substats,k==kval)[1,"lambda"]
      if(is.na(l)){
        rep(NA,length(true.labels))
      }else{
        alpha <- subset(path,lambda==l)[,attr(path,"alphacolnames")]
        assign.labels(alpha)
      }
    }
  }),list("gaussian mixture"=function(k){
    require(mclust)
    Mclust(x,k)$classification
  },"average linkage"=function(k){
    tree <- hclust(dist(x),"average")
    cutree(tree,k)
  },kmeans=function(k){
    sapply(1:100,function(seed){
      set.seed(seed)
      kmeans(x,k,100)$cluster
    })
  }))
  agree <- agreement.matrix(true.labels)
  do.call(rbind,lapply(names(FUNS),function(method){
    FUN <- FUNS[[method]]
    do.call(rbind,lapply(1:kmax,function(k){
      res <- as.matrix(FUN(k))
      rand <- apply(res,2,function(klass){
        norm.rand.index(klass,true.labels)
      })
      agreement <- apply(res,2,function(klass){
        sum((agreement.matrix(klass)-agree)^2)
      })
      data.frame(method,k,
                 mean.rand=mean(rand),sd.rand=sd(rand),
                 mean.agree=mean(agreement),sd.agree=sd(agreement))
    }))
  }))
}
halfcircle <- function(r,center=c(0,0),class,sign,N=150,noise=0.5){
  angle <- runif(N,0,pi)
  rad <- rnorm(N,r,noise)
  data.frame(x=rad*cos(angle)+center[1],y=sign*rad*sin(angle)+center[2],
             class=factor(class))
}
set.seed(1)
pts <- rbind(halfcircle(4,c(0,0),1,1),
             halfcircle(4,c(4,2),2,-1))
print(xyplot(y~x,pts,groups=class,aspect="iso"))
m <- pts[,c("x","y")]
moons.results <- calc.curves(m,pts$class)

## iris
irisSc <- scale(iris[,-5])
iris.results <- calc.curves(irisSc,iris[,5])

all.res <- rbind(
                 data.frame(iris.results,data="iris"),
                 data.frame(moons.results,data="moons"))
toplot <- subset(all.res,k<12 & k>1 & !is.nan(mean.rand))
levels(toplot$method) <- gsub(" ","\n",levels(toplot$method))
levels(toplot$method) <-
  c("$\\gamma=0.5$","$\\gamma=2$","$\\gamma=10$","GMM","HC","kmeans")
toplot$method <- reorder(toplot$method,-toplot$mean.rand,mean,na.rm=TRUE)
p <- ggplot(toplot,aes(k,mean.rand,colour=method,linetype=method))+
  geom_line(size=2)+
  geom_linerange(aes(ymin=mean.rand-sd.rand,ymax=mean.rand+sd.rand),size=2)+
  scale_y_continuous(paste("Normalized Rand index (bigger means",
                           "better agreement with known clusters)"),
                     breaks=c(0.4,0.6,0.8,1))+
  scale_x_continuous("Number of clusters",breaks=c(2,3,5,7,9,11))+
  coord_cartesian(xlim=c(2,11),ylim=c(0.2,1))+
  facet_grid(data~.,
             labeller=function(var,val)sprintf("%s: %s",var,val))+
  opts(title="Performance for several model sizes")
##library(directlabels)
##p <- xyplot(mean.rand~k|data,toplot,groups=method,type="l",
##            par.settings=simpleTheme(lwd=1.5,lty=1:5),xlim=c(2,14))
##p2 <- direct.label(p,list("\\gamma=10"))
print(p)
tikz("moons-iris.tex",height=9,width=7)
print(p+mytheme+opts(axis.title.y=theme_text(vjust=-0.1,angle=90)))
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
means <- data.frame(alpha=t(colMeans(sim$mat)))
p <- ggplot(cvx,aes(alpha.1,alpha.2))+
  ##geom_text(data=means,label="$\\overline X$",col="grey")+
  ##geom_point(aes(size=s),alpha=1/4)+
  geom_path(aes(colour=gamma,linetype=gamma,group=interaction(row,gamma)),lwd=1)+
  facet_grid_label(~norm)+
  coord_equal()+
  scale_size("$s$")+
  geom_point(data=data.frame(alpha=sim$mat),pch=21,fill="white")+
  scale_x_continuous("",breaks=-10)+
  scale_y_continuous("",breaks=-10)+
  opts(title="Norm and weights control the clusterpath")
print(p)
tikz("cvx-allnorms.tex",height=8,width=7)
print(p+theme_bw())
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
                                label=sprintf("\\tiny$w_{%s}$",ij)))
      polys <- rbind(polys,data.frame(alpha=m,rect.num=factor(k),
                    pair=ij,row.names=NULL))
    }
  }
  list(polys=polys,clusters=clusters,ws=ws,width=width)
}
idname <- "\\small Identity weights, $t=\\Omega(X)$"
line.df <- data.frame(getlines(x),figure=factor(idname))
path <- clusterpath.l2(x,verbose=1,join.thresh=0.01,gamma=1.3)
plot(path)
lval <- path$lambda[300]
poly.list <- getpolys(path,lval)
jointitle <- sprintf("\\scriptsize Decreasing weights after join, $t<\\Omega(X)$")
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
                    label=sprintf("\\scriptsize$\\ell_%s$",ann),
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
              label=sprintf("\\scriptsize$X_%d$",1:3),
              figure=rep(c(jointitle,idname,apart),each=3)))+
  facet_wrap("figure")+
  geom_text(aes(label=label),
            data.frame(alpha.1=c(-0.8,-0.5),alpha.2=c(-1.05,0.1),
                       figure=jointitle,
                       label=c("\\scriptsize$\\alpha_1$",
                         "\\scriptsize$\\alpha_{C}=\\alpha_2=\\alpha_3$")))
levs <- unique(unlist(lapply(p3$layers,function(l)levels(l$data$figure))))
allpts <- do.call(rbind,lapply(levs,function(l)data.frame(pts,figure=l)))
p4 <- p3+geom_point(data=allpts,pch=21,fill="white")+
 scale_x_continuous("",breaks=c(-10))+
 scale_y_continuous("",breaks=c(-10))+
 theme_bw()+
  opts(title="Geometric interpretation: constrain area between points")
print(p4)
tikz("geometry.tex",height=8,width=7.25)
print(p4)
dev.off()
