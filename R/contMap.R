## function plots reconstructed values for ancestral characters along the edges of the tree
## written by Liam J. Revell 2012-2023
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
	steps<-c(0:(res-1)/(res-1)*h,h+h/(res-1)) ## 0:res/res*(h+h/res)
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
	attr(tree,"map.order")<-"right-to-left"
	xx<-list(tree=tree,cols=cols,lims=lims)
	class(xx)<-"contMap"
	if(plot) plot.contMap(xx,fsize=fsize,ftype=ftype,lwd=lwd,legend=legend,outline=outline,
		sig=sig,type=type,mar=mar,direction=direction,offset=offset,hold=hold,leg.txt=leg.txt)
	invisible(xx)
}

## internally used function
## written by Liam J. Revell 2012, 2017
getState<-function(x,trans){
	if(x<=trans[1]) state<-names(trans)[1]
	else if(x>=trans[length(trans)]) state<-names(trans)[length(trans)]
	else {
		i<-1
		while(x>trans[i]){
			state<-names(trans)[i]
			i<-i+1
		}
	}
	state
}

## S3 print method for objects of class "contMap"
## uses print.densityMap internally
## written by Liam J. Revell 2012, 2013, 2014, 2015, 2016, 2023

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
	if(hasArg(sig)) sig<-list(...)$sig
	else sig<-3
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-NULL
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-NULL
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
	if(hasArg(leg.txt)) leg.txt<-list(...)$leg.txt
	else leg.txt<-"trait value"
	if(hasArg(underscore)) underscore<-list(...)$underscore
	else underscore<-FALSE
	if(hasArg(outline)) outline<-list(...)$outline
	else outline<-TRUE
	if(hasArg(nodes_only)) nodes_only<-list(...)$nodes_only
	else nodes_only<-FALSE
	if(hasArg(arc_height)) arc_height<-list(...)$arc_height
	else arc_height<-2
	if(is.null(legend)) legend<-if(type=="arc") max(H) else 0.5*max(H)
	if(is.null(fsize)) fsize<-c(1,1)
	if(length(fsize)==1) fsize<-rep(fsize,2)
	if(is.null(ftype)) ftype<-c("i","reg")
	if(length(ftype)==1) ftype<-c(ftype,"reg")
	if(nodes_only){
		if(hasArg(cex)) cex<-list(...)$cex
		else cex<-c(1.5,1)
		if(hasArg(lwd)) lwd<-list(...)$lwd
		else lwd<-c(1,4)
		obj<-x$tree
		N<-Ntip(obj)
		node.cols<-setNames(x$cols[
			c(names(obj$maps[[which(obj$edge[1,]==(N+1))]][1]),
			sapply(obj$maps,function(x) names(x[length(x)])))],
			c(N+1,obj$edge[,2]))
		node.cols<-node.cols[as.character(1:(N+obj$Nnode))]
		if(legend&&is.null(ylim)&&type%in%c("phylogram","cladogram")){
			if(direction%in%c("rightwards","leftwards")) ylim<-c(1-0.12*(N-1),N)
			else if(direction%in%c("upwards","downwards")) {
				pp<-par("pin")[2]
				sw<-(fsize*(max(strwidth(obj$tip.label,units="inches")))+
					1.37*fsize*strwidth("W",units="inches"))[1]
				alp<-optimize(function(a,H,sw,pp) (a*1.2*max(H)+sw-pp)^2,H=H,sw=sw,pp=pp,
					interval=c(0,1e6))$minimum
				ylim<-if(direction=="downwards") c(min(H)-sw/alp-0.16*max(H),max(H)) else 
					c(min(H)-0.16*max(H),max(H)+sw/alp)
			}
		} else if(is.null(ylim)) ylim<-NULL
		if(is.null(offset)) {
			if(type%in%c("cladogram","phylogram")) offset<-0.2*lwd[1]/3+0.2/3
			else if(type%in%c("fan","arc")) offset<-1
		}
		args<-list(...)
		args$ylim<-ylim
		args$tree<-as.phylo(obj)
		args$offset<-offset
		args$fsize<-fsize[1]
		args$lwd<-lwd[1]
		args$arc_height<-arc_height
		do.call(plotTree,args)
		pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		xx<-pp$xx
		yy<-pp$yy
		if(type%in%c("phylogram","cladogram")){
			if(direction=="rightwards")
				xx[1:N]<-xx[1:N]+strwidth(paste(obj$tip.label,"__",sep=""),cex=fsize[1])+offset
		}
		DROP<-if(type=="arc") DROP<-Ntip(x$tree)+1 else DROP<-NULL
		points(xx[-DROP],yy[-DROP],pch=if(outline) 21 else 16,
			col=if(outline) par()$fg else node.cols[-DROP],
			bg=if(outline) node.cols[-DROP] else NULL,
			cex=c(rep(cex[2],N),rep(cex[1],obj$Nnode-length(DROP))))
		if(legend){
			if(is.logical(legend)) legend<-0.5*max(H)
			if(length(leg.txt)==1) 
				leg.txt<-c(round(lims[1],sig),leg.txt,round(lims[2],sig))
			if(type%in%c("phylogram","cladogram")){
				if(direction%in%c("rightwards","leftwards")){
					add.color.bar(legend,x$cols,title=leg.txt[2],
						as.numeric(leg.txt[c(1,3)]),
						digits=sig,prompt=FALSE,
						x=if(direction=="leftwards") max(H)-legend else 0,
						y=1-0.08*(N-1),lwd=lwd[2],
						fsize=fsize[2],outline=outline,
						direction=if(!is.null(xlim)) 
						if(xlim[2]<xlim[1]) "leftwards" else 
						"rightwards" else "rightwards")
				} else if(direction%in%c("upwards","downwards")){
					add.color.bar(legend,x$cols,title=leg.txt[2],
						as.numeric(leg.txt[c(1,3)]),
						digits=sig,prompt=FALSE,x=1,
						y=ylim[1]+0.04*max(nodeHeights(x$tree)),
						lwd=lwd[2],outline=outline,
						fsize=fsize[2],direction="rightwards",
						subtitle=paste("length=",round(legend,3),
						sep=""))
				}
			} else if(type=="arc"){
				add.color.bar(legend,x$cols,
					title=leg.txt[2],
					as.numeric(leg.txt[c(1,3)]),
					digits=sig,prompt=FALSE,
					outline=outline,
					x=mean(par()$usr[1:2])-0.5*legend,
					y=par()$usr[3]+0.1*diff(par()$usr[3:4]),
					lwd=lwd[2],
					fsize=fsize[2])
			} else if(type=="fan"){
				add.color.bar(legend,x$cols,
					title=leg.txt[2],
					as.numeric(leg.txt[c(1,3)]),
					digits=sig,prompt=FALSE,
					outline=outline,
					x=0.9*par()$usr[1],
					y=0.9*par()$usr[3],
					lwd=lwd[2],
					fsize=fsize[2])
			}
		}
	} else {
		if(hasArg(lwd)) lwd<-list(...)$lwd
		else lwd<-4
		if(hasArg(hold)) hold<-list(...)$hold
		else hold<-TRUE
		# done optional arguments
		leg.txt<-c(round(lims[1],sig),leg.txt,round(lims[2],sig))
		plot(x,fsize=fsize,ftype=ftype,lwd=lwd,legend=legend,outline=outline,leg.txt=leg.txt,
			type=type,mar=mar,direction=direction,offset=offset,xlim=xlim,ylim=ylim,hold=hold,
			underscore=underscore,arc_height=arc_height)
	}
}

## S3 print method for object of class 'contMap'
## written by Liam J. Revell 2013
print.contMap<-function(x,digits=6,...){
	cat("Object of class \"contMap\" containing:\n\n")
	cat(paste("(1) A phylogenetic tree with ",length(x$tree$tip.label)," tips and ",x$tree$Nnode," internal nodes.\n\n",sep=""))
	cat(paste("(2) A mapped continuous trait on the range (",round(x$lims[1],digits),", ",round(x$lims[2],digits),").\n\n",sep="")) 
}

## drop tips from an object of class 'contMap'
## written by Liam J. Revell 2014, 2023
drop.tip.contMap<-function(phy,tip,...){
	if(!inherits(phy,"contMap")) cat("phy should be an object of class \"contMap\"\n")
	else {
		phy$tree<-drop.tip.simmap(phy$tree,tip,...)
		return(phy)
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

keep.tip.contMap<-function(phy,tip,...){
	tips<-setdiff(phy$tree$tip.label,tip)
	drop.tip.contMap(phy,tip=tips,...)
}
