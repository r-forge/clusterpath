path.animation <- structure(function
### Create an animation that summarizes the path calculation in
### 2d. This works best for the path algorithm, since each animation
### frame will represent a breakpoint in the path.
(df,
### Result of clusterpath.l1.id.
 ##outdir=tempfile(),
### Subdirectory for plot files, to be created.
 title="Calculation of the path in 2d"
### Title for the animation.
 ){
  require(animation)
  require(reshape2)
  ##dir.create(outdir)
  ## convert relative path to full path (animation package bug)
  ##outdir <- tools::file_path_as_absolute(outdir)
  lambdas <- sort(unique(df$lambda))
  lmax <- max(lambdas)
  toomany <- dcast(rename(predict(df),c(alpha="value")),...~col)
  bpts <- transform(castbreakpoints(df),lratio=lambda/lmax)
  ani.start(nmax=length(lambdas),
            title=title,
            ##outdir=outdir,
            ani.width=1000,
            ani.height=800,
            loop=FALSE,
            description="The clusterpath algorithm is calculated using identity weights and the l1 norm for a 2d data matrix.")
  on.exit(ani.stop())
  a <- do.call("aes_string",as.list(levels(df$col)))
  for(l in lambdas){
    ##cat(outdir,l,"\n")
    p <- ggplot()+
      a+
        coord_equal()+
        geom_path(aes(group=row),data=subset(toomany,lambda<=l),colour="grey")+
        geom_point(aes(size=lratio),data=subset(bpts,lambda<=l))+
        scale_size(expression(lambda/lambda[max]),limits=c(0,1))
    print(p)
  }
  ani.stop()
},ex=function(){
  ## several clusters
  sim <- gendata(5,2,2,0.1)
  colnames(sim$mat) <- c("height","length")
  plot(sim$mat)
  df <- clusterpath.l1.id(sim$mat)
  path.animation(df)

  iris2 <- iris[,c("Petal.Length","Petal.Width")]
  iris.path <- clusterpath.l1.id(iris2)
  path.animation(iris.path)
})
