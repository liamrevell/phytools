## function to manipulate & plot objects of class "backbonePhylo"
## written by Liam J. Revell 2013

## convert "phylo" to "backbonePhylo"
phylo.toBackbone<-function(x,trans,...){
	if(!inherits(x,"phylo")) stop("tree should be an object of class \"phylo\".")
	x$tip.clade<-list()
	for(i in 1:length(x$tip.label)){
		if(x$tip.label[i]%in%trans$tip.label){
			ii<-which(trans$tip.label==x$tip.label[i])
			x$tip.clade[[i]]<-list()
			x$tip.clade[[i]]$label<-trans$clade.label[ii]
			x$tip.clade[[i]]$N<-trans$N[ii]
			x$tip.clade[[i]]$depth<-trans$depth[ii]
		} else {
			x$tip.clade[[i]]<-list()
			x$tip.clade[[i]]$label<-x$tip.label[i]
			x$tip.clade[[i]]$N<-1
			x$tip.clade[[i]]$depth<-x$edge.length[which(x$edge[,2]==i)]
		}
	}
	x$tip.label<-NULL
	class(x)<-"backbonePhylo"
	x
}		

## convert to object of class "phylo"
backbone.toPhylo<-function(x){
	if(!inherits(x,"backbonePhylo")) stop("x not an object of class \"backbonePhylo\"")
	x$tip.label<-sapply(x$tip.clade,function(x) x$label)
	x$tip.clade<-NULL
	class(x)<-"phylo"
	x
}

## reorder backbone phylogeny
reorder.backbonePhylo<-function(x,order="cladewise",...){
	ii<-reorder(backbone.toPhylo(x),order,index.only=TRUE)
	x$edge<-x$edge[ii,]
	x$edge.length<-x$edge.length[ii]
	attr(x,"order")<-order
	x
}

## print method
print.backbonePhylo<-function(x,...){
	cat(paste("\nBackbone phylogenetic tree with",length(x$tip.clade),"subtrees and",x$Nnode,"resolved internal nodes.\n"))
	n<-min(length(x$tip.clade),5)
	cat("\nLabels: ")
	cat(paste(sapply(x$tip.clade[1:n],function(y) y$label),collapse=", "))
	cat(", ...\nDiversities: ")
	cat(paste(sapply(x$tip.clade[1:n],function(y) y$N),collapse=", "))
	cat(", ...\n\n")
}

## scale N
scaleN<-function(x,k){
	for(i in 1:length(x$tip.clade)) if(x$tip.clade[[i]]$N>1) x$tip.clade[[i]]$N<-x$tip.clade[[i]]$N*k
	x
}

## plot backbone phylogeny with triangles
plot.backbonePhylo<-function(x,...){
	if(!inherits(x,"backbonePhylo")) stop("x not an object of class \"backbonePhylo\"")
	if(hasArg(vscale)) vscale<-list(...)$vscale
	else vscale<-1
	x<-scaleN(x,vscale)
	tt<-backbone.toPhylo(x)
	n<-sum(sapply(x$tip.clade,function(x) x$N))
	cw<-reorder.backbonePhylo(x,"cladewise")
	y<-vector(length=length(cw$tip.clade)+cw$Nnode)
	z<-c(0,cumsum(sapply(x$tip.clade,function(x) x$N)))
	nn<-sapply(2:length(z),function(i,x) (x[i]-x[i-1])/2+x[i-1],x=z)
	y[cw$edge[cw$edge[,2]<=length(cw$tip.clade),2]]<-nn[1:length(cw$tip.clade)]
	pw<-reorder.backbonePhylo(x,"pruningwise")
	nn<-unique(pw$edge[,1])
	for(i in 1:length(nn)){
		yy<-y[pw$edge[which(pw$edge[,1]==nn[i]),2]]
		y[nn[i]]<-mean(range(yy))
	}
	# compute start & end points of each edge
	X<-nodeHeights(tt)
	# open & size a new plot
	plot.new(); par(mar=rep(0.1,4))
	pp<-par("pin")[1]
	sw<-par("cex")*(max(strwidth(sapply(cw$tip.clade,function(x) x$label),units="inches")))+1.37*par("cex")*strwidth("W",units="inches")
	alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=X,sw=sw,pp=pp,interval=c(0,1e6))$minimum
	xlim<-c(min(X),max(X)+sw/alp)
	plot.window(xlim=xlim,ylim=c(0,n))
	# plot horizontal edges
	for(i in 1:nrow(X)){
		if(cw$edge[i,2]>length(cw$tip.clade)) lines(X[i,],rep(y[cw$edge[i,2]],2),lwd=2,lend=2)
		else lines(c(X[i,1],X[i,2]-cw$tip.clade[[cw$edge[i,2]]]$depth),rep(y[cw$edge[i,2]],2),lwd=2,lend=2)
	}
	# plot vertical relationships
	for(i in 1:x$Nnode+length(x$tip.clade)) lines(X[which(cw$edge[,1]==i),1],range(y[cw$edge[which(cw$edge[,1]==i),2]]),lwd=2,lend=2)
	for(i in 1:length(x$tip.clade)){
		xx<-c(X[which(cw$edge[,2]==i),2]-cw$tip.clade[[i]]$depth,X[which(cw$edge[,2]==i),2],X[which(cw$edge[,2]==i),2])
		yy<-c(y[cw$edge[which(cw$edge[,2]==i),2]],y[cw$edge[which(cw$edge[,2]==i),2]]+cw$tip.clade[[i]]$N/2-0.5,y[cw$edge[which(cw$edge[,2]==i),2]]-cw$tip.clade[[i]]$N/2+0.5)
		if(yy[2]<yy[3]) yy[2]<-yy[3]<-yy[1]
		polygon(x=xx,y=yy,col="grey",lwd=2)
	}
	for(i in 1:length(cw$tip.clade)) text(X[which(cw$edge[,2]==i),2],y[i],cw$tip.clade[[i]]$label,pos=4,offset=0.1)
}



	

