## creates an object of class "cophylo"
## written by Liam J. Revell 2015, 2016
cophylo<-function(tr1,tr2,assoc=NULL,rotate=TRUE,...){
	if(!inherits(tr1,"phylo")||!inherits(tr2,"phylo")) 
		stop("tr1 & tr2 should be objects of class \"phylo\".")
	## hack to make sure tip labels of each tree are in cladewise order
	tr1<-untangle(tr1,"read.tree")
	tr2<-untangle(tr2,"read.tree")
	## if no association matrix check for exact matches
	if(is.null(assoc)){
		assoc<-intersect(tr1$tip.label,tr2$tip.label)
		assoc<-if(length(assoc)>0) cbind(assoc,assoc) else NULL
		if(is.null(assoc)){ 
			cat("No associations provided or found.\n")
			rotate<-FALSE
		}
	}
	## check to verify that all taxa in assoc are in tree
	ii<-sapply(assoc[,1],"%in%",tr1$tip.label)
	if(any(!ii)){ 
		assoc<-assoc[ii,]
		cat("Some species in assoc[,1] not in tr1. Removing species & links.\n")
	}
	ii<-sapply(assoc[,2],"%in%",tr2$tip.label)
	if(any(!ii)){ 
		assoc<-assoc[ii,]
		cat("Some species in assoc[,2] not in tr2. Removing species & links.\n")
	}
	## now check if rotation is to be performed
	if(rotate){
		cat("Rotating nodes to optimize matching...\n")
		flush.console()
		x<-setNames(sapply(assoc[,2],match,table=tr2$tip.label),assoc[,1])
		tr1<-tipRotate(tr1,x*Ntip(tr1)/Ntip(tr2),...)
		best.tr1<-Inf
		x<-setNames(sapply(assoc[,1],match,table=tr1$tip.label),assoc[,2])
		tr2<-tipRotate(tr2,x*Ntip(tr2)/Ntip(tr1),...)
		best.tr2<-Inf
		while((best.tr2-attr(tr2,"minRotate"))>0||(best.tr1-attr(tr1,"minRotate"))>0){
			best.tr1<-attr(tr1,"minRotate")
			best.tr2<-attr(tr2,"minRotate")
			x<-setNames(sapply(assoc[,2],match,table=tr2$tip.label),assoc[,1])
			tr1<-tipRotate(tr1,x*Ntip(tr1)/Ntip(tr2),...)
			x<-setNames(sapply(assoc[,1],match,table=tr1$tip.label),assoc[,2])
			tr2<-tipRotate(tr2,x*Ntip(tr2)/Ntip(tr1),...)
		}
		cat("Done.\n")
	}
	tt<-list(tr1,tr2)
	class(tt)<-"multiPhylo"
	obj<-list(trees=tt,assoc=assoc)
	class(obj)<-"cophylo"
	obj
}

## called internally by plot.cophylo to plot a phylogram
## written by Liam J. Revell
phylogram<-function(tree,part=1,direction="rightwards",fsize=1,ftype="i",lwd=1,...){
	if(hasArg(pts)) pts<-list(...)$pts
	else pts<-TRUE
	d<-if(direction=="rightwards") 1 else -1
	## check if edge lenths
	if(is.null(tree$edge.length)) tree<-compute.brlen(tree)
	## rescale tree so it fits in one half of the plot
	## with enough space for labels
	if(ftype=="off") fsize<-0
	sh<-fsize*strwidth(tree$tip.label)
	tree$edge.length<-tree$edge.length/max(nodeHeights(tree))*(part-max(sh))
	n<-Ntip(tree)
	## reorder cladewise to assign tip positions
	cw<-reorder(tree,"cladewise")
	y<-vector(length=n+cw$Nnode)
	y[cw$edge[cw$edge[,2]<=n,2]]<-0:(n-1)/(n-1)
	## reorder pruningwise for post-order traversal
	pw<-reorder(tree,"pruningwise")
	nn<-unique(pw$edge[,1])
	## compute vertical position of each edge
	for(i in 1:length(nn)){
		yy<-y[pw$edge[which(pw$edge[,1]==nn[i]),2]]
		y[nn[i]]<-mean(range(yy))
	}
	## compute start & end points of each edge
	X<-nodeHeights(cw)-0.5
	## plot horizontal edges
	for(i in 1:nrow(X)) lines(d*X[i,],rep(y[cw$edge[i,2]],2),lwd=lwd,lend=2)
	## plot vertical relationships
	for(i in 1:tree$Nnode+n) lines(d*X[which(cw$edge[,1]==i),1],
		sort(y[cw$edge[which(cw$edge[,1]==i),2]]),lwd=lwd,lend=2)
	## plot links to tips
	h<-max(X)+0.1*(max(X)-min(X))+max(fsize*strwidth(tree$tip.label))-
		fsize*strwidth(tree$tip.label)
	for(i in 1:n){ 
		lines(d*c(X[which(cw$edge[,2]==i),2],h[i]),rep(y[i],2),lwd=1,lty="dotted")
		if(pts) points(d*X[which(cw$edge[,2]==i),2],y[i],pch=16,cex=0.7*sqrt(lwd))
	}
	## plot tip labels
	font<-which(c("off","reg","b","i","bi")==ftype)-1
	if(font>0){
		for(i in 1:n) text(d*max(h+fsize*strwidth(tree$tip.label)),y[i],
			sub("_"," ",tree$tip.label[i]), pos=if(d<0) 4 else 2,offset=0,
			cex=fsize,font=font)
	}
	PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
		show.tip.label=if(ftype!="off") TRUE else FALSE,show.node.label=FALSE,
		font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=0,
		x.lim=par()$usr[1:2],y.lim=par()$usr[3:4],
		direction=direction,tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
		edge=cw$edge,xx=d*sapply(1:(Ntip(cw)+cw$Nnode),
		function(x,y,z) y[match(x,z)],y=X,z=cw$edge),yy=y)
	assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
	## return rightmost or leftmost edge of tip labels
	invisible(d*max(h+fsize*strwidth(tree$tip.label)))
}

## plot links between tip taxa according to assoc
## written by Liam J. Revell 2015
makelinks<-function(obj,x){
	for(i in 1:nrow(obj$assoc)){
		ii<-which(obj$trees[[1]]$tip.label==obj$assoc[i,1])
		jj<-which(obj$trees[[2]]$tip.label==obj$assoc[i,2])
		y<-c((ii-1)/(Ntip(obj$trees[[1]])-1),(jj-1)/(Ntip(obj$trees[[2]])-1))
		lines(x,y,lty="dashed")
	}
}

## plot an object of class "cophylo"
## written by Liam J. Revell 2015, 2016
plot.cophylo<-function(x,...){
	plot.new()
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-c(0.1,0.1,0.1,0.1)
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-c(-0.5,0.5)
	if(hasArg(scale.bar)) scale.bar<-list(...)$scale.bar
	else scale.bar<-rep(0,2)
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-if(any(scale.bar>0)) c(-0.1,1) else c(0,1)
	obj<-list(...)
	par(mar=mar)
	plot.window(xlim=xlim,ylim=ylim)
	leftArgs<-rightArgs<-obj
	if(!is.null(obj$fsize)){
		if(length(obj$fsize)>1){
			leftArgs$fsize<-obj$fsize[1]
			rightArgs$fsize<-obj$fsize[2]
			sb.fsize<- if(length(obj$fsize)>2) obj$fsize[3] else 1
		} else sb.fsize<-1
	} else sb.fsize<-1
	x1<-do.call("phylogram",c(list(tree=x$trees[[1]],part=0.4),leftArgs))
	left<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	x2<-do.call("phylogram",c(list(tree=x$trees[[2]],part=0.4,
		direction="leftwards"),rightArgs))
	right<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(!is.null(x$assoc)) makelinks(x,c(x1,x2))
	else cat("No associations provided.\n")
	if(any(scale.bar>0)) add.scalebar(x,scale.bar,sb.fsize)
	assign("last_plot.cophylo",list(left=left,right=right),envir=.PlotPhyloEnv)
}

## add scale bar
## written by Liam J. Revell 2015
add.scalebar<-function(obj,scale.bar,fsize){
	if(scale.bar[1]>0){
		s1<-(0.4-max(fsize*strwidth(obj$trees[[1]]$tip.label)))/max(nodeHeights(obj$trees[[1]]))
		lines(c(-0.5,-0.5+scale.bar[1]*s1),rep(-0.05,2))
		lines(rep(-0.5,2),c(-0.05,-0.06))
		lines(rep(-0.5+scale.bar[1]*s1,2),c(-0.05,-0.06))
		text(mean(c(-0.5,-0.5+scale.bar[1]*s1)),rep(-0.05,2),scale.bar[1],pos=1)
	}
	if(scale.bar[2]>0){
		s2<-(0.4-max(fsize*strwidth(obj$trees[[2]]$tip.label)))/max(nodeHeights(obj$trees[[2]]))
		lines(c(0.5-scale.bar[2]*s2,0.5),rep(-0.05,2))	
		lines(rep(0.5-scale.bar[2]*s2,2),c(-0.05,-0.06))
		lines(rep(0.5,2),c(-0.05,-0.06))
		text(mean(c(0.5-scale.bar[2]*s2,0.5)),rep(-0.05,2),scale.bar[2],pos=1)
	}
}

## print an object of class "cophylo"
## written by Liam J. Revell 2015
print.cophylo<-function(x, ...){
    cat("Object of class \"cophylo\" containing:\n\n")
    cat("(1) 2 (possibly rotated) phylogenetic trees in an object of class \"multiPhylo\".\n\n")
    cat("(2) A table of associations between the tips of both trees.\n\n")
}

## written by Liam J. Revell 2015
tipRotate<-function(tree,x,...){
	if(hasArg(fn)) fn<-list(...)$fn
	else fn<-function(x) x^2
	if(hasArg(methods)) methods<-list(...)$methods
	else methods<-"pre"
	if(hasArg(print)) print<-list(...)$print
	else print<-FALSE
	tree<-reorder(tree)
	nn<-1:tree$Nnode+length(tree$tip.label)
	if("pre"%in%methods){
		for(i in 1:tree$Nnode){
			tt<-read.tree(text=write.tree(rotate(tree,nn[i])))
			oo<-sum(fn(x-setNames(1:length(tree$tip.label),tree$tip.label)[names(x)]))
			pp<-sum(fn(x-setNames(1:length(tt$tip.label),tt$tip.label)[names(x)]))
			if(oo>pp) tree<-tt
			if(print) message(paste("objective:",min(oo,pp)))
		}
	}
	if("post"%in%methods){
		for(i in tree$Nnode:1){
			tt<-read.tree(text=write.tree(rotate(tree,nn[i])))
			oo<-sum(fn(x-setNames(1:length(tree$tip.label),tree$tip.label)[names(x)]))
			pp<-sum(fn(x-setNames(1:length(tt$tip.label),tt$tip.label)[names(x)]))
			if(oo>pp) tree<-tt
			if(print) message(paste("objective:",min(oo,pp)))
		}
	}
	attr(tree,"minRotate")<-min(oo,pp)
	tree
}

## labeling methods for plotted "cophylo" object
## written by Liam J. Revell 2015

nodelabels.cophylo<-function(...,which=c("left","right")){
	obj<-get("last_plot.cophylo",envir=.PlotPhyloEnv)
	if(which[1]=="left") assign("last_plot.phylo",obj[[1]],envir=.PlotPhyloEnv)
	else if(which[1]=="right") assign("last_plot.phylo",obj[[2]],envir=.PlotPhyloEnv)
	nodelabels(...)
}

edgelabels.cophylo<-function(...,which=c("left","right")){
	obj<-get("last_plot.cophylo",envir=.PlotPhyloEnv)
	if(which[1]=="left") assign("last_plot.phylo",obj[[1]],envir=.PlotPhyloEnv)
	else if(which[1]=="right") assign("last_plot.phylo",obj[[2]],envir=.PlotPhyloEnv)
	edgelabels(...)
}

tiplabels.cophylo<-function(...,which=c("left","right")){
	obj<-get("last_plot.cophylo",envir=.PlotPhyloEnv)
	if(which[1]=="left") assign("last_plot.phylo",obj[[1]],envir=.PlotPhyloEnv)
	else if(which[1]=="right") assign("last_plot.phylo",obj[[2]],envir=.PlotPhyloEnv)
	tiplabels(...)
}
