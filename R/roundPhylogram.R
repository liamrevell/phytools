## function plots a round phylogram
## written by Liam J. Revell 2014, 2015

roundPhylogram<-function(tree,fsize=1.0,ftype="reg",lwd=2,mar=NULL,offset=NULL,direction="rightwards",type="phylogram",xlim=NULL,ylim=NULL){
	if(inherits(tree,"multiPhylo")){
		par(ask=TRUE)
		tt<-lapply(tree,roundPhylogram,fsize=fsize,ftype=ftype,lwd=lwd,mar=mar,offset=offset, direction=direction,type=type,xlim=xlim,ylim=ylim)
	} else {
		if(type=="cladogram"||is.null(tree$edge.length)) tree<-compute.brlen(tree)
		ftype<-which(c("off","reg","b","i","bi")==ftype)-1
		if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
		# swap out "_" character for spaces (assumes _ is a place holder)
		tree$tip.label<-gsub("_"," ",tree$tip.label)
		if(is.null(mar)) mar=rep(0.1,4)
		n<-length(tree$tip.label)
		# set offset fudge (empirically determined)
		offsetFudge<-1.37
		# reorder cladewise to assign tip positions
		cw<-reorder(tree,"cladewise")
		y<-vector(length=n+cw$Nnode)
		y[cw$edge[cw$edge[,2]<=n,2]]<-1:n
		# reorder pruningwise for post-order traversal
		pw<-reorder(tree,"pruningwise")
		nn<-unique(pw$edge[,1])
		# compute vertical position of each edge
		for(i in 1:length(nn)){
			yy<-y[pw$edge[which(pw$edge[,1]==nn[i]),2]]
			y[nn[i]]<-mean(range(yy))
		}
		# compute start & end points of each edge
		X<-nodeHeights(cw)
		if(is.null(xlim)){
			pp<-par("pin")[1]
			sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+offsetFudge*fsize*strwidth("W",units="inches")
			alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=X,sw=sw,pp=pp,interval=c(0,1e6))$minimum
			xlim<-c(min(X),max(X)+sw/alp)
		}
		if(is.null(ylim)) ylim=c(1,max(y))
		## end preliminaries
		# open & size a new plot
		plot.new()
		par(mar=mar)
		if(is.null(offset)) offset<-0.2*lwd/3+0.2/3
		plot.window(xlim=xlim,ylim=ylim)
		# plot edges
		for(i in 1:nrow(X)){
			x<-NULL
			b<-y[cw$edge[i,1]]
			c<-X[i,1]
			d<-if(y[cw$edge[i,2]]>y[cw$edge[i,1]]) 1 else -1
			xx<-X[i,2]
			yy<-y[cw$edge[i,2]]
			a<-(xx-c)/(yy-b)^2
			curve(d*sqrt((x-c)/a)+b,from=X[i,1],to=X[i,2],add=TRUE,lwd=lwd)
		}
		# plot tip labels
		for(i in 1:n)
			if(ftype) text(X[which(cw$edge[,2]==i),2],y[i],tree$tip.label[i],pos=4,offset=offset,font=ftype,cex=fsize)
		PP<-list(type=type,use.edge.length=if(type=="phylogram") TRUE else FALSE,node.pos=1,
			show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,font=ftype,cex=fsize,
			adj=0,srt=0,no.margin=FALSE,label.offset=offset,x.lim=par()$usr[1:2],y.lim=par()$usr[3:4],
			direction=direction,tip.color="black",Ntip=length(cw$tip.label),Nnode=cw$Nnode,edge=cw$edge,
			xx=sapply(1:(length(cw$tip.label)+cw$Nnode),function(x,y,z) y[match(x,z)],y=X,z=cw$edge),
			yy=y)
		assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
	}
}
