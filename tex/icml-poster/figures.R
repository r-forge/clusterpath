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
 

