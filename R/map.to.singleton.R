## convert stochastic map style tree to a tree with singleton nodes
## written by Liam J. Revell 2013, 2015

map.to.singleton<-function(tree){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	tree<-reorderSimmap(tree)
	Nedge<-nrow(tree$edge)+sum(sapply(tree$maps,length)-1)
	Ntips<-length(tree$tip.label)
	edge<-tree$edge
	edge[edge>Ntips]<--edge[edge>Ntips]+Ntips
	xx<-vector(); edge.length<-vector(); ii<-1
	for(i in 1:nrow(tree$edge)){
		if(length(tree$maps[[i]])==1){ 
			xx[ii]<-names(tree$maps[[i]])
			edge.length[ii]<-tree$maps[[i]]
			ii<-ii+1
		} else {
			nn<-length(tree$maps[[i]])
			new<-matrix(NA,nn,2)
			new[1,1]<-edge[ii,1]
			nextnode<--1+if(i>1) min(c(edge[1:(ii-1),],edge[ii,1])) else edge[ii,1]
			new[,2]<-nextnode-0:(nn-1)
			for(j in 2:nn) new[j,1]<-new[j-1,2]
			if(edge[ii,2]>0) new[nrow(new),2]<-edge[ii,2]
			if(i==nrow(tree$edge)) edge<-rbind(edge[1:(ii-1),],new)
			else {
				ee<-edge[(ii+1):nrow(edge),]
				ee[ee<edge[ii,1]]<-ee[ee<edge[ii,1]]-nrow(new)+1
				if(ii==1) edge<-rbind(new,ee)
				else edge<-rbind(edge[1:(ii-1),],new,ee)
			}
			xx[1:length(tree$maps[[i]])+ii-1]<-names(tree$maps[[i]])
			edge.length[1:length(tree$maps[[i]])+ii-1]<-tree$maps[[i]]
			ii<-ii+length(tree$maps[[i]])
		}
	}
	tip.label<-tree$tip.label
	Nnode<--min(edge)
	edge[edge<0]<-length(tip.label)-edge[edge<0]
	tree<-list(edge=edge,Nnode=Nnode,tip.label=tip.label,edge.length=setNames(edge.length,xx))
	class(tree)<-c("singleton","phylo")
	attr(tree,"order")<-"cladewise"
	return(tree)
}

## function plots tree with singleton nodes
## written by Liam J. Revell 2013, 2015, 2017
plotTree.singletons<-function(tree){
	## preliminaries
	n<-length(tree$tip.label)
	if(is.null(names(tree$edge.length))) names(tree$edge.length)<-rep(1,length(tree$edge.length))
	colors<-setNames(palette()[1:length(unique(names(tree$edge.length)))],sort(unique(names(tree$edge.length))))
	colors<-colors[names(tree$edge.length)]
	# assume order is cladewise (otherwise we're in trouble!)
	if(attr(tree,"order")=="cladewise") cw<-tree
	else stop("tree must be in \"cladewise\" order.")
	y<-vector(length=n+cw$Nnode)
	y[cw$edge[cw$edge[,2]<=n,2]]<-1:n
	# reorder pruningwise for post-order traversal
	pw<-reorderPhylo(cw,"pruningwise")
	nn<-unique(pw$edge[,1])
	# compute vertical position of each edge
	for(i in 1:length(nn)){
		yy<-y[pw$edge[which(pw$edge[,1]==nn[i]),2]]
		y[nn[i]]<-mean(range(yy))
	}
	# compute start & end points of each edge
	X<-nodeHeights(cw)
	# open & size a new plot
	plot.new(); par(mar=c(1.1,1.1,0.1,0.1))
	pp<-par("pin")[1]
	sw<-par("cex")*(max(strwidth(cw$tip.label,units="inches")))+1.37*par("cex")*strwidth("W",units="inches")
	alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=X,sw=sw,pp=pp,interval=c(0,1e6))$minimum
	xlim<-c(min(X),max(X)+sw/alp)
	plot.window(xlim=xlim,ylim=c(1,max(y)))
	axis(1,at=c(0,max(X)),labels=FALSE); axis(2,at=1:max(y),labels=FALSE)
	# plot horizontal edges
	for(i in 1:nrow(X)) lines(X[i,],rep(y[cw$edge[i,2]],2),lwd=2,lend=2,col=colors[i])
	# plot vertical relationships
	for(i in 1:tree$Nnode+n){
		xx<-X[which(cw$edge[,1]==i),1]
		yy<-y[cw$edge[which(cw$edge[,1]==i),2]]
		if(length(xx)>1) lines(xx,yy,lwd=2,lend=2,col=colors[which(cw$edge[,1]==i)])
	}
	# plot points
	for(i in 1:nrow(X)) points(X[i,],rep(y[cw$edge[i,2]],2),pch=21,bg="gray")
	# plot tip labels
	for(i in 1:n) text(X[which(cw$edge[,2]==i),2],y[i],tree$tip.label[i],pos=4,offset=0.3)
	PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
		show.tip.label=TRUE,show.node.label=FALSE,
		font=par()$font,cex=par()$cex,adj=0,srt=0,no.margin=FALSE,label.offset=0.3,
		x.lim=xlim,y.lim=c(1,max(y)),
		direction="rightward",tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
		edge=cw$edge,xx=sapply(1:(Ntip(cw)+cw$Nnode),
		function(x,y,z) y[match(x,z)],y=X,z=cw$edge),yy=y)
	assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
}

## function to reorder the edges of the tree for postorder traversal should *only* be used 
## internally to plotTree.singletons (i.e., is not designed for general use).
## written by Liam J. Revell 2013
reorderPhylo<-function(x,order="pruningwise",index.only=FALSE,...){
	if(index.only) stop("index.only=TRUE not permitted")
	if(order!="pruningwise") stop("function only returns trees in modified pruningwise format for plotTree.singletons")
	if(attr(x,"order")!="cladewise") stop("input tree should be in cladewise order.")
	aa<-lapply(x$edge[,1],getDescendants,tree=x)
	ll<-sapply(aa,length)
	ii<-order(ll)
	x$edge<-x$edge[ii,]
	x$edge.length<-x$edge.length[ii]
	attr(x,"order")<-"pruningwise"
	return(x)
}

## function converts a tree with a root edge to a tree with a singleton node instead
## written by Liam J. Revell 2016, re-written 2019
rootedge.to.singleton<-function(tree){
	if(!inherits(tree,"phylo"))
		stop("tree should be object of class \"phylo\".")
	if(!is.null(tree$root.edge)){
		tree$edge[tree$edge>Ntip(tree)]<-
			tree$edge[tree$edge>Ntip(tree)]+1
		if(attr(tree,"order")%in%c("postorder","pruningwise")){
			tree$edge<-rbind(tree$edge,c(1,2)+Ntip(tree))
			tree$edge.length<-c(tree$edge.length,tree$root.edge)
		} else {
			tree$edge<-rbind(c(1,2)+Ntip(tree),tree$edge)
			tree$edge.length<-c(tree$root.edge,tree$edge.length)
		}
		tree$root.edge<-NULL
		tree$Nnode<-tree$Nnode+1
		if(!is.null(tree$node.label)) 
			tree$node.label<-c("",tree$node.label)
	}
	tree
}



