## phylomorphospace3d: projection of a tree into three dimensional morphospace
## written by Liam J. Revell 2012, 2013, 2014, 2016, 2018

phylomorphospace3d<-function(tree,X,A=NULL,label=TRUE,control=list(),method=c("dynamic","static"),...){
	method<-method[1]
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	# control list
	con=list(spin=TRUE,axes=TRUE,box=TRUE,simple.axes=FALSE,lwd=1,ftype="reg",
		col.edge=rep("black",nrow(tree$edge)))
	con[(namc<-names(control))]<-control
	if(con$simple.axes) con$box<-con$axes<-FALSE
	con$ftype<-which(c("off","reg","b","i","bi")==con$ftype)-1
	if(is.null(A)) A<-apply(X,2,function(x,tree) fastAnc(tree,x),tree=tree)
	else A<-A[as.character(1:tree$Nnode+length(tree$tip)),]
	x<-y<-z<-matrix(NA,nrow(tree$edge),2)
	X<-X[tree$tip.label,]
	for(i in 1:length(tree$tip)){
		x[tree$edge[,2]==i,2]<-X[i,1]
		y[tree$edge[,2]==i,2]<-X[i,2]
		z[tree$edge[,2]==i,2]<-X[i,3]
	}
	for(i in length(tree$tip)+1:tree$Nnode){
		x[tree$edge[,1]==i,1]<-x[tree$edge[,2]==i,2]<-A[as.character(i),1]
		y[tree$edge[,1]==i,1]<-y[tree$edge[,2]==i,2]<-A[as.character(i),2]
		z[tree$edge[,1]==i,1]<-z[tree$edge[,2]==i,2]<-A[as.character(i),3]
	}
	if(is.null(colnames(X))) colnames(X)<-c("x","y","z")
	if(method=="dynamic"){
		chk<-.check.pkg("rgl")
		if(!chk){
			cat("  method = \"dynamic\" requires the package \"rgl\"\n  Defaulting to method = \"static\"\n\n")
			method<-"static"
			lines3d<-play3d<-plot3d<-segments3d<-spheres3d<-spin3d<-text3d<-function(...) NULL
		}
	}
	if(method=="dynamic"){
		params<-get("r3dDefaults")
		plot3d(rbind(X,A),xlab=colnames(X)[1],ylab=colnames(X)[2],zlab=colnames(X)[3],
			axes=con$axes,box=con$box,params=params)
		spheres3d(X,radius=0.02*mean(apply(X,2,max)-apply(X,2,min)))
		for(i in 1:nrow(tree$edge)) segments3d(x[i,],y[i,],z[i,],lwd=con$lwd,
			col=con$col.edge[i])
		ms<-colMeans(X)
		rs<-apply(rbind(X,A),2,range)
		if(con$simple.axes){
			lines3d(x=rs[,1],y=c(rs[1,2],rs[1,2]),z=c(rs[1,3],rs[1,3]))
			lines3d(x=c(rs[1,1],rs[1,1]),y=rs[,2],z=c(rs[1,3],rs[1,3]))
			lines3d(x=c(rs[1,1],rs[1,1]),y=c(rs[1,2],rs[1,2]),z=rs[,3])
		}
		rs<-rs[2,]-rs[1,]
		for(i in 1:length(tree$tip)){
			adj<-0.03*rs*(2*(X[i,]>ms)-1)
			if(con$ftype) text3d(X[i,]+adj,texts=tree$tip.label[i],font=con$ftype)
		}
		if(con$spin){
			xx<-spin3d(axis=c(0,0,1),rpm=10)
			play3d(xx,duration=5)
			invisible(xx)
		}
		else invisible(NULL)
	} else if(method=="static"){
		if(hasArg(angle)) angle<-list(...)$angle
		else angle<-30
		if(hasArg(xlim)) xlim<-list(...)$xlim
		else xlim<-NULL
		if(hasArg(ylim)) ylim<-list(...)$ylim
		else ylim<-NULL
		if(hasArg(zlim)) zlim<-list(...)$zlim
		else zlim=NULL
		xx<-scatterplot3d(X,xlab=colnames(X)[1],zlab=colnames(X)[3],pch=19,angle=angle,
			ylab=colnames(X)[2],cex.symbols=1.3,xlim=xlim,ylim=ylim,zlim=zlim)
		aa<-xx$xyz.convert(A)
		points(aa$x,aa$y,pch=19,cex=0.8)
		for(i in 1:nrow(tree$edge)){
			aa<-xx$xyz.convert(x[i,],y[i,],z[i,])
			lines(aa$x,aa$y,lwd=2,col=con$col.edge[i])
		}
		for(i in 1:length(tree$tip.label)){
			aa<-xx$xyz.convert(x[which(tree$edge[,2]==i),2],y[which(tree$edge[,2]==i),2],
				z[which(tree$edge[,2]==i),2])
			if(con$ftype) text(tree$tip.label[i],x=aa$x,y=aa$y,pos=2,font=con$ftype)
		}
		invisible(xx)
	}
}

## check for a package (modified from 'geiger')
## primarily used for rgl
## written by Liam J. Revell 2014, 2022
.check.pkg<-function(pkg){
	if(pkg%in%rownames(installed.packages())){
		require(pkg,character.only=TRUE)
		return(TRUE)
	} else return(FALSE)
}
