## function depends on phytools (& dependencies) and maps (& dependencies)
## written by Liam J. Revell 2013, 2017, 2019

phylo.to.map<-function(tree,coords,rotate=TRUE,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	# optional arguments
	if(hasArg(database)) database<-list(...)$database
	else database<-"world"
	if(hasArg(regions)) regions<-list(...)$regions
	else regions<-"."
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-c(-180,180)
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-c(-90,90)
	# create a map
	map<-map(database,regions,xlim=xlim,ylim=ylim,plot=FALSE,fill=TRUE,resolution=0)
	# if rotate
	if(hasArg(type)) type<-list(...)$type else type<-"phylogram"
	if(hasArg(direction)) direction<-list(...)$direction else direction<-"downwards"
	if(is.data.frame(coords)) coords<-as.matrix(coords)
	if(rotate&&type=="phylogram"){
		cc<-aggregate(coords,by=list(rownames(coords)),mean)
		cc<-matrix(as.matrix(cc[,2:3]),nrow(cc),2,dimnames=list(cc[,1],colnames(cc)[2:3]))
		tree<-minRotate(tree,cc[,if(direction=="rightwards") 1 else 2])
	} else direction<-"unoptimized"
	x<-list(tree=tree,map=map,coords=coords,direction=direction)
	class(x)<-"phylo.to.map"
	if(plot) plot.phylo.to.map(x,...)
	invisible(x)
}

## S3 method to plot object of class "phylo.to.map"
## written by Liam J. Revell 2013, 2014, 2016, 2019

plot.phylo.to.map<-function(x,type=c("phylogram","direct"),...){
	type<-type[1]
	if(class(x)=="phylo.to.map"){
		tree<-x$tree
		map<-x$map
		coords<-x$coords
	} else stop("x should be an object of class \"phylo.to.map\"")
	# get optional arguments
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-map$range[1:2]
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-map$range[3:4]
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-1.0
	if(hasArg(split)) split<-list(...)$split
	else split<-c(0.4,0.6)
	if(hasArg(psize)) psize<-list(...)$psize
	else psize<-1.0
	if(hasArg(cex.points)){ 
		cex.points<-list(...)$cex.points
		if(length(cex.points)==1) cex.points<-c(0.6*cex.points,cex.points)
	} else cex.points<-c(0.6*psize,psize)
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-rep(0,4)
	if(hasArg(asp)) asp<-list(...)$asp
	else asp<-1.0
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-"reg"
	ftype<-which(c("off","reg","b","i","bi")==ftype)-1
	if(!ftype) fsize=0 
	if(hasArg(from.tip)) from.tip<-list(...)$from.tip
	else from.tip<-FALSE
	if(hasArg(colors)) colors<-list(...)$colors
	else colors<-"red"
	if(length(colors)==1) rep(colors[1],2)->colors
	if(length(colors)==2&&type=="phylogram"){
		colors<-matrix(rep(colors,nrow(coords)),nrow(coords),2,byrow=TRUE)
		rownames(colors)<-rownames(coords)
	} else if(is.vector(colors)&&(length(colors)==Ntip(tree))) {
		COLS<-matrix("red",nrow(coords),2,dimnames=list(rownames(coords)))
		for(i in 1:length(colors)) COLS[which(rownames(COLS)==names(colors)[i]),1:2]<-colors[i]
		colors<-COLS
	}
	if(hasArg(direction)) direction<-list(...)$direction
	else direction<-"downwards"
	if(hasArg(pch)) pch<-list(...)$pch
	else pch<-21
	if(length(pch)==1) pch<-rep(pch,2)
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-c(2,1)
	if(length(lwd)==1) lwd<-rep(lwd,2)
	if(hasArg(lty)) lty<-list(...)$lty
	else lty<-"dashed"
	if(hasArg(pts)) pts<-list(...)$pts
	else pts<-TRUE
	if(type=="phylogram"){
		if(x$direction=="downwards"&&direction=="rightwards"){
			cat("\"phylo.to.map\" direction is \"downwards\" but plot direction has been given as \"rightwards\".\n")
			cat("Re-optimizing object....\n")
			cc<-aggregate(coords,by=list(rownames(coords)),mean)
			cc<-matrix(as.matrix(cc[,2:3]),nrow(cc),2,dimnames=list(cc[,1],colnames(cc)[2:3]))
			tree<-minRotate(tree,cc[,1])
		} else if(x$direction=="rightwards"&&direction=="downwards"){
			cat("\"phylo.to.map\" direction is \"rightwards\" but plot direction has been given as \"downwards\".\n")
			cat("Re-optimizing object....\n")
			cc<-aggregate(coords,by=list(rownames(coords)),mean)
			cc<-matrix(as.matrix(cc[,2:3]),nrow(cc),2,dimnames=list(cc[,1],colnames(cc)[2:3]))
			tree<-minRotate(tree,cc[,2])
		}
	}
	# recompute ylim or xlim to leave space for the tree
	if(type=="phylogram"){
		if(direction=="downwards"){
			if(!ftype) ylim<-c(ylim[1],ylim[2]+0.03*diff(ylim))
			ylim<-c(ylim[1],ylim[2]+
				split[1]/split[2]*(ylim[2]-ylim[1]))
		} else if(direction=="rightwards"){
			if(!ftype) xlim<-c(xlim[1]-0.03*diff(xlim),xlim[2])
			xlim<-c(xlim[1]-split[1]/split[2]*(xlim[2]-xlim[1]),xlim[2])
		}
	}
	# open & size a new plot
	#par(mar=mar)
	#plot.new()
	#plot.window(xlim=c(-1,1),ylim=c(-1,1),asp=1)
	#abline(h=0,col=make.transparent("blue",0.5))
	#stop("to here?")
	if(all(mar==0)) mar<-mar+0.01
	plot.new()
	par(mar=mar)
	plot.window(xlim=xlim,ylim=ylim,asp=asp)
	map(map,add=TRUE,fill=TRUE,col="gray95",mar=rep(0,4))
	if(type=="phylogram"){
		## preliminaries
		cw<-reorder(tree,"cladewise")
		n<-Ntip(cw)
		if(direction=="downwards"){
			# plot a white rectangle
			dx<-abs(diff(xlim))
			rect(xlim[1]-1.04*dx,ylim[2]-split[1]*(ylim[2]-ylim[1]),
				xlim[2]+1.04*dx,ylim[2],col="white",border="white")
			# rescale tree so it fits in the upper half of the plot
			# with enough space for labels
			pdin<-par()$din[2]
			sh<-(fsize*strwidth(paste(" ",cw$tip.label,sep=""))+
				0.3*fsize*strwidth("W"))*(par()$din[1]/par()$din[2])*
				(diff(par()$usr[3:4])/diff(par()$usr[1:2]))
			cw$edge.length<-cw$edge.length/max(nodeHeights(cw))*
				(split[1]*(ylim[2]-ylim[1])-max(sh))
			pw<-reorder(cw,"postorder") ## post-order
			x<-vector(length=n+cw$Nnode)
			x[cw$edge[cw$edge[,2]<=n,2]]<-0:(n-1)/(n-1)*(xlim[2]-xlim[1])+xlim[1]
			nn<-unique(pw$edge[,1])
			# compute horizontal position of each edge
			for(i in 1:length(nn)){
				xx<-x[pw$edge[which(pw$edge[,1]==nn[i]),2]]
				x[nn[i]]<-mean(range(xx))
			}
			# compute start & end points of each edge
			Y<-ylim[2]-nodeHeights(cw)
			# plot coordinates & linking lines
			coords<-coords[,2:1]
			for(i in 1:nrow(coords)){ 
				tip.i<-which(cw$tip.label==rownames(coords)[i])
				lines(c(x[tip.i],coords[i,1]),c(Y[which(cw$edge[,2]==tip.i),2]-
					if(from.tip) 0 else sh[tip.i],coords[i,2]),
					col=colors[i,1],lty=lty,lwd=lwd[2])
			}
			points(coords,pch=pch,cex=cex.points[2],bg=colors[,2])
			# plot vertical edges
			for(i in 1:nrow(Y)) lines(rep(x[cw$edge[i,2]],2),Y[i,],
				lwd=lwd[1],lend=2)
			# plot horizontal relationships
			for(i in 1:cw$Nnode+n) 
				lines(range(x[cw$edge[which(cw$edge[,1]==i),2]]),
				Y[which(cw$edge[,1]==i),1],lwd=lwd[1],lend=2)
			# plot tip labels
			for(i in 1:n){ 
				if(ftype) text(x[i],Y[which(cw$edge[,2]==i),2],
					paste(" ",sub("_"," ",cw$tip.label[i]),sep=""),
						pos=4,offset=c(0,1),
					srt=-90,cex=fsize,font=ftype)
				if(pts) points(x[i],Y[which(cw$edge[,2]==i),2],pch=21,
					bg=colors[cw$tip.label,][i,2],
					cex=cex.points[1])
			}
			PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
				show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
				font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,
				label.offset=fsize*strwidth(" ")/(par()$usr[2]-
				par()$usr[1])*(par()$usr[4]-par()$usr[3]),
				x.lim=par()$usr[1:2],y.lim=par()$usr[3:4],
				direction=direction,tip.color="black",Ntip=Ntip(cw),
				Nnode=cw$Nnode,edge=cw$edge,xx=x,yy=sapply(1:(Ntip(cw)+
				cw$Nnode),function(x,y,z) y[match(x,z)],y=Y,z=cw$edge))
		} else {
			dy<-abs(diff(ylim))
			rect(xlim[1],ylim[1],xlim[1]+split[1]*(xlim[2]-
				xlim[1]),ylim[2],col="white",border="white")
			sh<-fsize*strwidth(paste(" ",cw$tip.label,sep=""))+
				0.2*fsize*strwidth("W")
			cw$edge.length<-cw$edge.length/max(nodeHeights(cw))*
				(split[1]*(xlim[2]-xlim[1])-max(sh))
			pw<-reorder(cw,"postorder")
			y<-vector(length=n+cw$Nnode)
			y[cw$edge[cw$edge[,2]<=n,2]]<-0:(n-1)/(n-1)*(ylim[2]-ylim[1])+ylim[1]
			nn<-unique(pw$edge[,1])
			for(i in 1:length(nn)){
				yy<-y[pw$edge[which(pw$edge[,1]==nn[i]),2]]
				y[nn[i]]<-mean(range(yy))
			}
			H<-nodeHeights(cw)
			X<-xlim[1]+H
			coords<-coords[,2:1]
			for(i in 1:nrow(coords)){
				tip.i<-which(cw$tip.label==rownames(coords)[i])
				lines(c(X[which(cw$edge[,2]==tip.i),2]+if(from.tip) 0 else sh[tip.i],coords[i,1]),
					c(y[tip.i],coords[i,2]),col=colors[i,1],lty=lty,lwd=lwd[2])
			}
			points(coords,pch=pch,cex=cex.points[2],bg=colors[,2])
			for(i in 1:nrow(X)) lines(X[i,],rep(y[cw$edge[i,2]],2),lwd=lwd[1],lend=2)
			for(i in 1:cw$Nnode+n) lines(X[which(cw$edge[,1]==i),1],
				range(y[cw$edge[which(cw$edge[,1]==i),2]]),lwd=lwd[1],lend=2)
			for(i in 1:n){
				if(ftype) text(X[which(cw$edge[,2]==i),2],y[i],
					paste(" ",sub("_"," ",cw$tip.label[i]),sep=""),
					pos=4,offset=0.1,cex=fsize,font=ftype)
				if(pts) points(X[which(cw$edge[,2]==i),2],y[i],
					pch=21,bg=colors[cw$tip.label,][i,2],cex=cex.points[1])
			}
			PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
				show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
				font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=0.1,
				x.lim=par()$usr[1:2],y.lim=par()$usr[3:4],
				direction=direction,tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
				edge=cw$edge,xx=sapply(1:(Ntip(cw)+cw$Nnode),
				function(x,y,z) y[match(x,z)],y=X,z=cw$edge),yy=y)
		}
		assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
	} else if(type=="direct"){
		phylomorphospace(tree,coords[,2:1],add=TRUE,label="horizontal",
			node.size=c(0,psize),lwd=lwd[2],control=list(col.node=
			setNames(rep(colors[2],max(tree$edge)),1:max(tree$edge)),
			col.edge=setNames(rep(colors[1],nrow(tree$edge)),tree$edge[,2])),
			ftype=c("off","reg","b","i","bi")[ftype+1],fsize=fsize)
	}
}


## rotates all nodes to try and match tip an ordering 
## written by Liam J. Revell 2013, 2015
minRotate<-function(tree,x,...){
	if(hasArg(print)) print<-list(...)$print
	else print<-TRUE
	tree<-reorder(tree)
	nn<-1:tree$Nnode+Ntip(tree)
	x<-x[tree$tip.label]
	for(i in 1:tree$Nnode){
		tt<-read.tree(text=write.tree(rotate(tree,nn[i])))
		oo<-sum(abs(order(x[tree$tip.label])-1:length(tree$tip.label)))
		pp<-sum(abs(order(x[tt$tip.label])-1:length(tt$tip.label)))
		if(oo>pp) tree<-tt
		if(print) message(paste("objective:",min(oo,pp)))
	}
	attr(tree,"minRotate")<-min(oo,pp)
	return(tree)
}

print.phylo.to.map<-function(x,...){
	cat("Object of class \"phylo.to.map\" containing:\n\n")
	cat(paste("(1) A phylogenetic tree with",Ntip(x$tree),"tips and",x$tree$Nnode,"internal nodes.\n\n",sep=" "))
	cat("(2) A geographic map with range:\n")
	cat(paste("     ",paste(round(x$map$range[3:4],2),collapse="N, "),"N\n",sep=""))
	cat(paste("     ",paste(round(x$map$range[1:2],2),collapse="W, "),"W.\n\n",sep=""))
	cat(paste("(3) A table containing ",nrow(x$coords)," geographic coordinates (may include\n",
			"    more than one set per species).\n\n",sep=""))
	if(x$direction%in%c("downwards","rightwards")) 
		cat(paste("If optimized, tree nodes have been rotated to maximize alignment\n",
			"with the map when the tree is plotted in a ",x$direction," direction.\n\n",sep=""))
	else
		cat("The nodes of the tree may not have yet been rotated to maximize\nalignment between the phylogeny & the map.\n\n")
}
