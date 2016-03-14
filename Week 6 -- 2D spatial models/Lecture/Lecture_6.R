
setwd( "C:/Users/James.Thorson/Desktop/UW Hideaway/Course plan 2016 -- spatiotemporal models/Week 6 -- 2D spatial models/Lecture" )

############################
# Show SPDE approximation
############################

library( INLA )  # FROM: http://www.r-inla.org/download

# Simulate locations
loc_xy = cbind( "x"=runif(10), "y"=runif(10))

# create mesh
mesh = inla.mesh.create( loc_xy, plot.delay=NULL, refine=FALSE)

# Create matrices in INLA
spde <- inla.spde2.matern(mesh, alpha=2)

png(file="SPDE_explanation_pt1.png", width=8, height=4, res=200, units="in")
  par( mfrow=c(1,2), mar=c(0,0,2,0), mgp=c(2,0.5,0), tck=-0.02)
  # Plot samples
  plot( loc_xy, xlim=range(mesh$loc[,1]), ylim=range(mesh$loc[,2]), main="Sample locations")
  # Plot mesh
  plot(mesh, main="Mesh composed of triangles")
  text( x=mesh$loc[,1], y=mesh$loc[,2], labels=1:mesh$n, col=ifelse(1:mesh$n%in%mesh$idx$loc,"blue","black"))
  title("Mesh composed of triangles")
dev.off()

png(file="SPDE_explanation_pt2.png", width=9, height=3, res=200, units="in")
  par( mfrow=c(1,3), mar=c(2,2,2,0), mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs='i')
  # Visualize SPDE approx.
  Col = colorRampPalette(colors=c("blue","white","red"))
  Points = function(X,Y,Z){
    DF = cbind( expand.grid('X'=as.vector(X), 'Y'=as.vector(Y)), 'Z'=as.vector(Z) )
    DF = DF[which(DF[,'Z']!=0),]
    points( x=DF[,'X'], y=DF[,'Y'], pch=20)
  }
  image(z=as.matrix(spde$param.inla$M0), x=1:mesh$n, y=1:mesh$n, main="M0", zlim=c(-1,1)*max(abs(spde$param.inla$M0)), col=Col(11), xlab="", ylab=""); box()
  Points( X=1:mesh$n, Y=1:mesh$n, Z=as.matrix(spde$param.inla$M0))
  image(z=as.matrix(spde$param.inla$M1), x=1:mesh$n, y=1:mesh$n, main="M1", zlim=c(-1,1)*max(abs(spde$param.inla$M1)), col=Col(11), xlab="", ylab=""); box()
  Points( X=1:mesh$n, Y=1:mesh$n, Z=as.matrix(spde$param.inla$M1))
  image(z=as.matrix(spde$param.inla$M2), x=1:mesh$n, y=1:mesh$n, main="M2", zlim=c(-1,1)*max(abs(spde$param.inla$M2)), col=Col(11), xlab="", ylab=""); box()
  Points( X=1:mesh$n, Y=1:mesh$n, Z=as.matrix(spde$param.inla$M2))
dev.off()

