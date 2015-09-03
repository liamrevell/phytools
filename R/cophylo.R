## creates an object of class "cophylo"
## written by Liam J. Revell 2015
cophylo<-function(tr1,tr2,assoc=NULL,rotate=TRUE,...){
	if(!inherits(tr1,"phylo")||!inherits(tr2,"phylo")) 
		stop("tr1 & tr2 should be objects of class \"phylo\".")
	## hack to make sure tip labels of each tree are in cladewise order
	tr1<-read.tree(text=write.tree(tr1))
	tr2<-read.tree(text=write.tree(tr2))
	## if no association matrix check for exact matches
	if(is.null(assoc)){
		assoc<-intersect(tr1$tip.label,tr2$tip.label)
		assoc<-if(length(assoc)>0) cbind(assoc,assoc) else NULL
		if(is.null(assoc)){ 
			cat("No associations provided or found.\n")
			rotate<-FALSE
		}
	}
	## now check if rotation is to be performed
	if(rotate){
		cat("Rotating nodes to optimize matching...\n")
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
phylogram<-function(tree,part=1,direction="right",fsize=1,ftype="i",lwd=1){
	d<-if(direction=="right") 1 else -1
	## rescale tree so it fits in one half of the plot
	## with enough space for labels
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
		points(d*X[which(cw$edge[,2]==i),2],y[i],pch=16,cex=0.7*sqrt(lwd))
	}
	## plot tip labels
	for(i in 1:n) text(d*max(h+fsize*strwidth(tree$tip.label)),y[i],
		sub("_"," ",tree$tip.label[i]), pos=if(d<0) 4 else 2,offset=0,
		cex=fsize,font=which(c("off","reg","b","i","bi")==ftype)-1)
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
## written by Liam J. Revell 2015
plot.cophylo<-function(x,...){
	plot.new()
	par(mar=c(0.1,0.1,0.1,0.1))
	plot.window(xlim=c(-0.5,0.5),ylim=c(0,1))
	x1<-phylogram(x$trees[[1]],part=0.4,...)
	x2<-phylogram(x$trees[[2]],part=0.4,direction="left",...)
	if(!is.null(x$assoc)) makelinks(x,c(x1,x2))
	else cat("No associations provided.\n")
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
