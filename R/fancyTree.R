# function to plot special types of phylogeny visualizations
# so far the implemented types are:
# "extinction" in which all branches leading to extinct taxa (or prior to the MRCA of extant species) are plotted with red dashed lines;
# "traitgram3d" which creates a 3D graph projecting the tree into two-dimensional morphospace (with time as the third axis)
# "droptip" creates a two panel plot with the tips to be pruned marked (panel 1) and then removed, and returns the pruned tree
# "xkcd" creates an xkcd-comic style phylogeny [no longer an option]
# "densitymap" maps the posterior density of a binary stochastic character mapping
# "contmap" maps reconstructed trait evolution for a continuous character on the tree
# "phenogram95" plots a 95% CI phenogram
# "scattergram" plots a phylogenetic scatterplot matrix
# written by Liam J. Revell 2012, 2013, 2014, 2015

fancyTree<-function(tree,type=c("extinction","traitgram3d","droptip","densitymap","contmap","phenogram95","scattergram"),...,control=list()){
	type<-matchType(type,c("extinction","traitgram3d","droptip","densitymap","contmap","phenogram95","scattergram"))
	if(!inherits(tree,"phylo")&&type%in%c("extinction","traitgram3d","droptip")) stop("tree should be an object of class \"phylo\".")
	else if(!inherits(tree,"multiPhylo")&&type=="densitymap") stop("for type='densitymap' tree should be an object of class \"multiPhylo\".")
	if(type=="extinction") extinctionTree(tree)
	else if(type=="traitgram3d") invisible(traitgram3d(tree,...,control=control))
	else if(type=="droptip") return(droptipTree(tree,...))
	else if(type=="densitymap") plotDensityMap(tree,...)
	else if(type=="contmap") plotContMap(tree,...)
	else if(type=="phenogram95") phenogram95(tree,...)
	else if(type=="scattergram") phyloScattergram(tree,...)
	else stop(paste("do not recognize type = \"",type,"\"",sep=""))
}

# phyloScattergram internal function
# written by Liam J. Revell 2013, 2014

phyloScattergram<-function(tree,...){
	if(hasArg(X)) X<-list(...)$X
	else stop("phenotypic data should be provided in the matrix X")
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-0.7
	if(hasArg(colors)) colors<-list(...)$colors
	else if(!is.null(tree$maps)) colors<-setNames(palette()[1:ncol(tree$mapped.edge)],sort(colnames(tree$mapped.edge)))
	if(hasArg(label)) label<-list(...)$label
	else label<-"radial"
	if(hasArg(hold)) hold<-list(...)$hold
	else hold<-TRUE
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	m<-ncol(X)
	if(hold) null<-dev.hold()
	if(!quiet&&hold){ 
		cat("Computing multidimensional phylogenetic scatterplot matrix...\n")
		flush.console()
	}
	par(mfrow=c(m,m))
	par(cex=fsize)
	par(mar=c(0,0,0,0))
	par(oma=c(5,5,3,3))
	m<-ncol(X)
	for(i in 1:m) for(j in 1:m){
		if(i==j) contMap(tree,X[,i],legend=FALSE,lwd=2,outline=F,fsize=fsize)
		else { 
			phylomorphospace(tree,X[,c(j,i)],lwd=1,node.by.map=TRUE,axes=FALSE,node.size=c(0,1),colors=colors,label=label,xlab="",ylab="")
			if(i==1) axis(side=3) # top row
			if(i==m) axis(side=1) # first column
			if(j==1) axis(side=2) # bottom row
			if(j==m) axis(side=4) # last column
		}
	}
	par(cex=0.9)
	if(is.null(colnames(X))) colnames(X)<-paste("V",1:m,sep="")
	invisible(mapply(title,xlab=colnames(X),adj=seq(0,(m-1)/m,1/m)+1/(2*m),MoreArgs=list(outer=TRUE,cex=0.9)))
	invisible(mapply(title,ylab=colnames(X)[m:1],adj=seq(0,(m-1)/m,1/m)+1/(2*m),MoreArgs=list(outer=TRUE,cex=0.9)))
	if(hold) null<-dev.flush()
}

# phenogram95 internal function
# written by Liam J. Revell 2013, 2014

phenogram95<-function(tree,...){
	if(hasArg(x)) x<-list(...)$x
	else stop("no phenotypic data provided")
	if(hasArg(spread.labels)) spread.labels<-list(...)$spread.labels
	else spread.labels<-TRUE
	if(hasArg(link)) link<-list(...)$link
	else link<-0.05*max(nodeHeights(tree))
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-0
	if(hasArg(hold)) hold<-list(...)$hold
	else hold<-TRUE
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	# get ancestral states
	A<-fastAnc(tree,x,CI=TRUE)
	# compute transparencies
	if(hasArg(tlim)) tlim<-list(...)$tlim
	else tlim<-c(0,25)
	trans<-as.character(floor(seq(tlim[1],tlim[2],length.out=51)))
	trans[as.numeric(trans)<10]<-paste("0", trans[as.numeric(trans)<10],sep="")
	# now get the arguments for phenogram
	args<-list(...)
	args$tree<-tree
	args$lwd<-1
	args$link<-0.05*max(nodeHeights(tree))
	args$offset<-0
	if(hold) null<-dev.hold()
	if(!quiet&&hold){ 
		cat("Computing density traitgram...\n")
		flush.console()
	}
	for(i in 0:50){
  		p<-i/length(trans)
		args$add<-i>0
		args$spread.labels<-if(i==0) spread.labels else FALSE
		args$link<-if(i==0) link else 0
		args$offset<-if(i==0) offset else offset+link
		args$x<-c(x,(1-p)*A$CI95[,1]+p*A$ace)
		args$colors<-paste("#0000ff",trans[i+1],sep="")
		do.call(phenogram,args)
		args$x<-c(x,(1-p)*A$CI95[,2]+p*A$ace)
		args$add<-TRUE
		args$spread.labels<-FALSE
		args$link<-0
		args$offset<-offset+link
		do.call(phenogram,args)
	}
	args$x<-c(x,A$ace)
	args$add<-TRUE
	args$colors<-"white"
	args$lwd<-2
	args$offset<-offset+link
	do.call(phenogram,args)
	null<-dev.flush()
}

# extinctionTree internal function
# written by Liam J. Revell 2012

extinctionTree<-function(tree){
	edges<-rep(0,nrow(tree$edge))
	names(edges)<-tree$edge[,2]
	extant<-getExtant(tree)
	ca<-findMRCA(tree,extant)
	root.node<-length(tree$tip)+1
	if(ca!=root.node){
		z<-setdiff(getDescendants(tree,root.node),getDescendants(tree,ca))
		edges[as.character(z)]<-1
	}
	z<-getDescendants(tree,ca)
	y<-lapply(z,getDescendants,tree=tree)
	for(i in 1:length(z)) if(!any(tree$tip.label[y[[i]]]%in%extant)) edges[as.character(z[i])]<-1
	plot.phylo(tree,edge.color=edges+1,edge.lty=edges+1,edge.width=2,no.margin=TRUE)
}

# traitgram3d internal function
# written by Liam J. Revell 2012, 2013

traitgram3d<-function(tree,...,control){
	if(hasArg(X)) X<-list(...)$X
	else stop("no phenotypic data provided")
	if(!hasArg(A)){
		if(is.null(control$maxit)) maxit<-2000
		else maxit<-control$maxit
		Y<-apply(X,2,function(x,tree) anc.ML(tree,x,maxit),tree=tree)
		convergence<-sapply(Y,function(x) x$convergence)
		if(any(convergence!=0)) warning("anc.ML may not have converged; consider increasing maxit.")
		A<-sapply(Y,function(x) x$ace)
	} else { 
		A<-list(...)$A
		A<-A[as.character(1:tree$Nnode+length(tree$tip)),]
	}
	if(is.null(colnames(X))) colnames(X)<-c("x","y")
	X<-cbind(X,diag(vcv(tree))[rownames(X)])
	A<-cbind(A,nodeHeights(tree)[match(rownames(A)[1:nrow(A)],tree$edge)])
	colnames(X)[3]<-colnames(A)[3]<-"time"
	# other optional arguments
	if(hasArg(method)) method<-list(...)$method
	else method<-"dynamic"
	if(hasArg(angle)) angle<-list(...)$angle
	else angle<-30
	# done other optional arguments
	xx<-phylomorphospace3d(tree,X,A,control=control,method=method,angle=angle,zlim=range(nodeHeights(tree)))
	return(xx)
}

# droptipTree internal function
# written by Liam J. Revell 2012

droptipTree<-function(tree,...){
	if(hasArg(tip)) tip<-list(...)$tip
	else stop("need to provide tip or tips to drop")
	edges<-rep(0,nrow(tree$edge))
	names(edges)<-tree$edge[,2]
	keep<-setdiff(tree$tip.label,tip)
	ca<-findMRCA(tree,keep)
	root.node<-length(tree$tip)+1
	if(ca!=root.node){
		z<-setdiff(getDescendants(tree,root.node),getDescendants(tree,ca))
		edges[as.character(z)]<-1
	}
	z<-getDescendants(tree,ca)
	foo<-function(x,tree){
		n<-length(tree$tip.label)
		y<-getDescendants(tree,x)
		y<-y[y<=n]
		return(y)
	}
	y<-lapply(z,foo,tree=tree)
	for(i in 1:length(z)) if(!any(tree$tip.label[y[[i]]]%in%keep)) edges[as.character(z[i])]<-1
	par(mfrow=c(2,1))
	plot.phylo(tree,edge.color=edges+1,edge.lty=edges+1,edge.width=2,no.margin=TRUE)
	dtree<-drop.tip(tree,tip); dtree$root.edge<-max(nodeHeights(tree))-max(nodeHeights(dtree))
	plot.phylo(dtree,edge.width=2,no.margin=TRUE,root.edge=TRUE)
	return(dtree)
}

# plotDensityMap internal function
# written by Liam J. Revell 2012

plotDensityMap<-function(trees,...){
	if(hasArg(res)) res<-list(...)$res
	else res<-100
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-NULL
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-NULL
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-3
	if(hasArg(check)) check<-list(...)$check
	else check<-FALSE
	if(hasArg(legend)) legend<-list(...)$legend
	else legend<-NULL
	if(hasArg(outline)) outline<-list(...)$outline
	else outline<-FALSE
	densityMap(trees,res,fsize,check,legend,outline)
}

# plotContMap internal function
# written by Liam J. Revell 2012

plotContMap<-function(tree,...){
	if(hasArg(x)) x<-list(...)$x
	else stop("need to provide vector 'x' of phenotypic trait values")
	if(hasArg(res)) res<-list(...)$res

	else res<-100
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-NULL
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-NULL
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-4
	if(hasArg(legend)) legend<-list(...)$legend
	else legend<-NULL
	if(hasArg(lims)) lims<-list(...)$lims
	else lims<-NULL
	if(hasArg(outline)) outline<-list(...)$outline
	else outline<-TRUE
	if(hasArg(sig)) sig<-list(...)$sig
	else sig<-3
	contMap(tree,x,res,fsize,ftype,lwd,legend,lims,outline,sig)
}


