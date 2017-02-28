## function plots reconstructed values for ancestral characters along the edges of the tree
## written by Liam J. Revell 2012-2017
contMap<-function(tree,x,res=100,fsize=NULL,ftype=NULL,lwd=4,legend=NULL,
	lims=NULL,outline=TRUE,sig=3,type="phylogram",direction="rightwards",
	plot=TRUE,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-rep(0.3,4)
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-NULL
	if(hasArg(method)) method<-list(...)$method
	else method<-"fastAnc"
	if(hasArg(hold)) hold<-list(...)$hold
	else hold<-TRUE
	if(hasArg(leg.txt)) leg.txt<-list(...)$leg.txt
	else leg.txt<-"trait value"
	h<-max(nodeHeights(tree))
	steps<-0:res/res*max(h)
	H<-nodeHeights(tree)
	if(method=="fastAnc") a<-fastAnc(tree,x) 
	else if(method=="anc.ML") { 
		fit<-anc.ML(tree,x)
		a<-fit$ace
		if(!is.null(fit$missing.x)) x<-c(x,fit$missing.x)
	} else if(method=="user"){
		if(hasArg(anc.states)) anc.states<-list(...)$anc.states
		else {
			cat("No ancestral states have been provided. Using states estimated with fastAnc.\n\n")
			a<-fastAnc(tree,x)
		}
		if(length(anc.states)<tree$Nnode){
			nodes<-as.numeric(names(anc.states))
			tt<-tree
			for(i in 1:length(nodes)){
				M<-matchNodes(tt,tree,method="distances")
				ii<-M[which(M[,2]==nodes[i]),1]
				tt<-bind.tip(tt,nodes[i],edge.length=0,where=ii)
			}
			a<-fastAnc(tt,c(x,anc.states))
			M<-matchNodes(tree,tt,method="distances")
			a<-a[as.character(M[,2])]
			names(a)<-M[,1]
		} else { 
			if(is.null(names(anc.states))) names(anc.states)<-1:tree$Nnode+Ntip(tree)
			a<-anc.states[as.character(1:tree$Nnode+Ntip(tree))]
		}
	}
	y<-c(a,x[tree$tip.label])
	names(y)[1:Ntip(tree)+tree$Nnode]<-1:Ntip(tree)
	A<-matrix(y[as.character(tree$edge)],nrow(tree$edge),ncol(tree$edge))
	cols<-rainbow(1001,start=0,end=0.7)
	names(cols)<-0:1000
	if(is.null(lims)) lims<-c(min(y),max(y))
	trans<-0:1000/1000*(lims[2]-lims[1])+lims[1]
	names(trans)<-0:1000
	tree$maps<-vector(mode="list",length=nrow(tree$edge))
	for(i in 1:nrow(tree$edge)){
		XX<-cbind(c(H[i,1],steps[intersect(which(steps>H[i,1]),which(steps<H[i,2]))]),
			c(steps[intersect(which(steps>H[i,1]),which(steps<H[i,2]))],H[i,2]))-H[i,1]
		YY<-rowMeans(XX)
		if(!all(YY==0)){
			b<-vector()
			for(j in 1:length(YY))
				b[j]<-(A[i,1]/YY[j]+A[i,2]/(max(XX)-YY[j]))/(1/YY[j]+1/(max(XX)-YY[j]))
		} else b<-A[i,1]
		d<-sapply(b,getState,trans=trans)
		tree$maps[[i]]<-XX[,2]-XX[,1]
		names(tree$maps[[i]])<-d
	}
	tree$mapped.edge<-makeMappedEdge(tree$edge,tree$maps)
	tree$mapped.edge<-tree$mapped.edge[,order(as.numeric(colnames(tree$mapped.edge)))]
	class(tree)<-c("simmap",setdiff(class(tree),"simmap"))
	xx<-list(tree=tree,cols=cols,lims=lims)
	class(xx)<-"contMap"
	if(plot) plot.contMap(xx,fsize=fsize,ftype=ftype,lwd=lwd,legend=legend,outline=outline,
		sig=sig,type=type,mar=mar,direction=direction,offset=offset,hold=hold,leg.txt=leg.txt)
	invisible(xx)
}

# function
# written by Liam J. Revell 2012
getState<-function(x,trans){
	i<-1
	while(x>trans[i]){
		state<-names(trans)[i]
		i<-i+1
	}
	return(state)
}

## S3 print method for objects of class "contMap"
## uses print.densityMap internally
## written by Liam J. Revell 2012, 2013, 2014, 2015, 2016

plot.contMap<-function(x,...){
	if(inherits(x,"contMap")){
		lims<-x$lims
		x<-list(tree=x$tree,cols=x$cols)
		class(x)<-"densityMap"
	} else stop("x should be an object of class \"contMap\"")
	H<-nodeHeights(x$tree)
	# get & set optional arguments
	if(hasArg(legend)) legend<-list(...)$legend
	else legend<-NULL
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-NULL
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-NULL
	if(hasArg(outline)) outline<-list(...)$outline
	else outline<-TRUE
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-4
	if(hasArg(sig)) sig<-list(...)$sig
	else sig<-3
	if(hasArg(type)) type<-list(...)$type
	else type<-"phylogram"
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-rep(0.3,4)
	if(hasArg(direction)) direction<-list(...)$direction
	else direction<-"rightwards"
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-NULL
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-NULL
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-NULL
	if(hasArg(hold)) hold<-list(...)$hold
	else hold<-TRUE
	if(is.null(legend)) legend<-0.5*max(H)
	if(is.null(fsize)) fsize<-c(1,1)
	if(length(fsize)==1) fsize<-rep(fsize,2)
	if(is.null(ftype)) ftype<-c("i","reg")
	if(length(ftype)==1) ftype<-c(ftype,"reg")
	if(hasArg(leg.txt)) leg.txt<-list(...)$leg.txt
	else leg.txt<-"trait value"
	# done optional arguments
	leg.txt<-c(round(lims[1],sig),leg.txt,round(lims[2],sig))
	plot(x,fsize=fsize,ftype=ftype,lwd=lwd,legend=legend,outline=outline,leg.txt=leg.txt,
		type=type,mar=mar,direction=direction,offset=offset,xlim=xlim,ylim=ylim,hold=hold)
}

## S3 print method for object of class 'contMap'
## written by Liam J. Revell 2013
print.contMap<-function(x,digits=6,...){
	cat("Object of class \"contMap\" containing:\n\n")
	cat(paste("(1) A phylogenetic tree with ",length(x$tree$tip.label)," tips and ",x$tree$Nnode," internal nodes.\n\n",sep=""))
	cat(paste("(2) A mapped continuous trait on the range (",round(x$lims[1],digits),", ",round(x$lims[2],digits),").\n\n",sep="")) 
}

## drop tips from an object of class 'contMap'
## written by Liam J. Revell 2014
drop.tip.contMap<-function(x,tip){
	if(!inherits(x,"contMap")) cat("x should be an object of class \"contMap\"\n")
	else {
		x$tree<-drop.tip.simmap(x$tree,tip)
		return(x)
	}
}

## add error bars to contMap plot
## written by Liam J. Revell 2017
errorbar.contMap<-function(obj,...){
	if(hasArg(x)) x<-list(...)$x
	else x<-setNames(sapply(1:Ntip(obj$tree),function(x,obj){
		ii<-which(obj$tree$edge[,2]==x)
		ss<-names(obj$tree$maps[[ii]][length(obj$tree$maps[[ii]])])
		obj$lims[1]+as.numeric(ss)/(length(obj$cols)-1)*diff(obj$lims)
		},obj=obj),obj$tree$tip.label)
	if(hasArg(scale.by.ci)) scale.by.ci<-list(...)$scale.by.ci
	else scale.by.ci<-TRUE
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-14
	tree<-obj$tree
	aa<-fastAnc(tree,x,CI=TRUE)
	xlim<-range(aa$CI95)
	if(xlim[2]>obj$lims[2]||xlim[1]<obj$lims[1]){
		cat(paste("  -----\n  The range of the contMap object, presently (",
			round(obj$lims[1],4),",",
			round(obj$lims[2],4),
			"), should be equal to\n  or greater than the range of the CIs on ancestral states: (",
			round(xlim[1],4),",",round(xlim[2],4),").\n",sep=""))
		cat(paste("  To ensure that your error bars are correctly plotted, please recompute your\n", 
			"  contMap object and increase lims.\n  -----\n",sep=""))
	}
	d<-diff(obj$lims)
	if(scale.by.ci){
		v<-aa$CI95[,2]-aa$CI95[,1]
		v<-v/max(v)
	} else v<-rep(0.5,tree$Nnode)
	n<-length(obj$cols)-1
	lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	h<-max(nodeHeights(tree))
	for(i in 1:tree$Nnode){
		ii<-round((aa$CI95[i,1]-obj$lims[1])/d*n)
		jj<-round((aa$CI95[i,2]-obj$lims[1])/d*(n+1))
		cols<-obj$cols[ii:jj]   
		add.color.bar(leg=0.1*h*v[i],cols=cols,prompt=FALSE,
			x=lastPP$xx[i+Ntip(tree)]-0.05*h*v[i],
			y=lastPP$yy[i+Ntip(tree)],title="",subtitle="",lims=NULL,
			lwd=lwd)
	}
}
