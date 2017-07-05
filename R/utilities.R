## some utility functions
## written by Liam J. Revell 2011, 2012, 2013, 2014, 2015, 2016, 2017

## function to expand clades in a plot by a given factor
## written by Liam J. Revell 2017
expand.clade<-function(tree,node,factor=5){
	cw<-reorder(tree)
	tips<-setNames(rep(1,Ntip(tree)),cw$tip.label)
	get.tips<-function(node,tree){
    		dd<-getDescendants(tree,node)
    		tree$tip.label[dd[dd<=Ntip(tree)]]
	}
	desc<-unlist(lapply(node,get.tips,tree=cw))
	for(i in 2:Ntip(cw)){
		tips[i]<-tips[i-1]+
			if(names(tips)[i]%in%desc){
				1 
			} else if(names(tips)[i-1]%in%desc){
				1
			} else 1/factor
	}
	obj<-list(tree=tree,tips=tips)
	class(obj)<-"expand.clade"
	obj
}

## S3 print method for the object class "expand.clade"
print.expand.clade<-function(x,...){
	cat("An object of class \"expand.clade\" consisting of:\n")
	cat(paste("(1) A phylogenetic tree (x$tree) with",Ntip(x$tree),
		"tips and\n   ",x$tree$Nnode,"internal nodes.\n"))
	cat("(2) A vector (x$tips) containing the desired tip-spacing.\n\n")
}

## S3 plot method for the object class "expand.clade"
plot.expand.clade<-function(x,...){
	args<-list(...)
	args$tree<-x$tree
	args$tips<-x$tips
	if(inherits(args$tree,"simmap")) do.call(plotSimmap,args)
	else do.call(plotTree,args)
}

## function to add a geological or other temporal legend to a plotted tree
## written by Liam J. Revell 2017
geo.legend<-function(leg=NULL,colors=NULL,alpha=0.2,...){
	if(hasArg(cex)) cex<-list(...)$cex
	else cex<-par()$cex
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	if(hasArg(show.lines)) show.lines<-list(...)$show.lines
	else show.lines<-TRUE
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(is.null(colors)){
		colors<-setNames(c(
			rgb(255,242,127,255,maxColorValue=255),
			rgb(255,230,25,255,maxColorValue=255),
			rgb(253,154,82,255,maxColorValue=255),
			rgb(127,198,78,255,maxColorValue=255),
			rgb(52,178,201,255,maxColorValue=255),
			rgb(129,43,146,255,maxColorValue=255),
			rgb(240,64,40,255,maxColorValue=255),
			rgb(103,165,153,255,maxColorValue=255),
			rgb(203,140,55,255,maxColorValue=255),
			rgb(179,225,182,255,maxColorValue=255),
			rgb(0,146,112,255,maxColorValue=255),
			rgb(127,160,86,255,maxColorValue=255),
			rgb(247,67,112,255,maxColorValue=255)),
			c("Quaternary","Neogene","Paleogene",
			"Cretaceous","Jurassic","Triassic",
			"Permian","Carboniferous","Devonian",
			"Silurian","Ordovician","Cambrian",
			"Precambrian"))
	}
	if(is.null(leg)){
		leg<-rbind(c(2.588,0),
			c(23.03,2.588),
			c(66.0,23.03),
			c(145.0,66.0),
			c(201.3,145.0),
			c(252.17,201.3),
			c(298.9,252.17),
			c(358.9,298.9),
			c(419.2,358.9),
			c(443.8,419.2),
			c(485.4,443.8),
			c(541.0,485.4),
			c(4600,541.0))
		rownames(leg)<-c("Quaternary","Neogene","Paleogene",
			"Cretaceous","Jurassic","Triassic",
			"Permian","Carboniferous","Devonian",
			"Silurian","Ordovician","Cambrian",
			"Precambrian")
		t.max<-max(obj$xx)
		ii<-which(leg[,2]<=t.max)
		leg<-leg[ii,]
		leg[max(ii),1]<-t.max
	}
	colors<-sapply(colors,make.transparent,alpha=alpha)
	if(plot){	
		y<-c(rep(0,2),rep(par()$usr[4],2))
		ylabel<--1/25*obj$Ntip
		for(i in 1:nrow(leg)){
			strh<-strheight(rownames(leg)[i])
			polygon(c(leg[i,1:2],leg[i,2:1]),y,
				col=colors[rownames(leg)[i]],border=NA)
			if(show.lines){
				lines(x=rep(leg[i,1],2),y=c(0,par()$usr[4]),
					lty="dotted",col="grey")
				lines(x=c(leg[i,1],mean(leg[i,])-0.8*cex*
					get.asp()*strheight(rownames(leg)[i])),
					y=c(0,ylabel),lty="dotted",col="grey")
				lines(x=c(leg[i,2],mean(leg[i,])+0.8*cex*
					get.asp()*strheight(rownames(leg)[i])),
					y=c(0,ylabel),lty="dotted",col="grey")
				lines(x=rep(mean(leg[i,])-0.8*cex*
					get.asp()*strheight(rownames(leg)[i]),2),
					y=c(ylabel,par()$usr[3]),lty="dotted",col="grey")
				lines(x=rep(mean(leg[i,])+0.8*cex*
					get.asp()*strheight(rownames(leg)[i]),2),
					y=c(ylabel,par()$usr[3]),lty="dotted",col="grey")
			}
			polygon(x=c(leg[i,1],
				mean(leg[i,])-0.8*cex*get.asp()*strh,
				mean(leg[i,])-0.8*cex*get.asp()*strh,
				mean(leg[i,])+0.8*cex*get.asp()*strh,
				mean(leg[i,])+0.8*cex*get.asp()*strh,
				leg[i,2]),y=c(0,ylabel,par()$usr[3],
				par()$usr[3],ylabel,0),
				col=colors[rownames(leg)[i]],border=NA)
			text(x=mean(leg[i,])+
				if(obj$direction=="leftwards") 0.12*strh else -0.12*strh,
				y=ylabel,labels=rownames(leg)[i],
				srt=90,adj=c(1,0.5),cex=cex)
		}
	}
	invisible(list(leg=leg,colors=colors))
}

## borrowed from mapplots
get.asp<-function(){
	pin<-par("pin")
	usr<-par("usr")
	asp<-(pin[2]/(usr[4]-usr[3]))/(pin[1]/(usr[2]-usr[1]))
	asp
}

round.polygon<-function(x,y,col="transparent"){
	## just space holding for now	
}

## draw a box around a clade
## written by Liam J. Revell 2017

cladebox<-function(tree,node,color=NULL,...){
	if(is.null(color)) color<-make.transparent("yellow",0.2)
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	h<-max(nodeHeights(tree))
	parent<-tree$edge[which(tree$edge[,2]==node),1]
	x0<-max(c(obj$xx[node]+obj$xx[parent])/2,obj$xx[node]-0.05*h)
	x1<-obj$x.lim[2]
	dd<-getDescendants(tree,node)
	y0<-min(range(obj$yy[dd]))-0.5
	y1<-max(range(obj$yy[dd]))+0.5
	polygon(c(x0,x1,x1,x0),c(y0,y0,y1,y1),col=color,
		border=0)
}

## draw tip labels as linking lines to text
## written by Liam J. Revell 2017

linklabels<-function(text,tips,link.type=c("bent","curved","straight"),
	...){
	lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(!(lastPP$direction%in%c("leftwards","rightwards")))
		stop("direction should be \"rightwards\" or \"leftwards\".")
	if(hasArg(cex)) cex<-list(...)$cex
	else cex<-1
	if(hasArg(col)) col<-list(...)$col
	else col<-"black"
	if(hasArg(lty)) lty<-list(...)$lty
	else lty<-"dashed"
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-1
	if(hasArg(link.offset)) link.offset<-list(...)$link.offset
	else link.offset<-0.1*max(lastPP$xx)
	if(hasArg(font)) font<-list(...)$font
	else font<-3
	link.type<-link.type[1]
	xpos<-lastPP$xx[tips]+strwidth("i")
	ypos<-lastPP$yy[tips]
	xmax<-rep(max(lastPP$xx),length(tips))+link.offset
	ylab<-seq(min(lastPP$yy),max(lastPP$yy),
		by=(max(lastPP$yy)-min(lastPP$yy))/(length(tips)-1))
	ylab<-ylab[rank(ypos)]
	text(xmax,ylab,gsub("_"," ",text),pos=4,font=font,cex=cex,
		offset=0)
	if(link.type=="curved"){
		for(i in 1:length(tips))
			drawCurve(c(xpos[i],xmax[i]),c(ypos[i],ylab[i]),
				scale=0.05,lty=lty,col=col,lwd=lwd)
	} else if(link.type=="bent"){
		tipmax<-max(lastPP$xx)
		for(i in 1:length(tips)){
			ff<-strwidth("W")
			segments(xpos[i],ypos[i],tipmax+link.offset/2,ypos[i],
				lty=lty,col=col,lwd=lwd)
			segments(tipmax+link.offset/2,ypos[i],tipmax+
				link.offset/2+ff,ylab[i],lty=lty,col=col,lwd=lwd)
			segments(tipmax+link.offset/2+ff,ylab[i],xmax[i],ylab[i],
				lty=lty,col=col,lwd=lwd)
		}
	} else if(link.type=="straight")
		segments(xpos,ypos,xmax,ylab,lty=lty,col=col)
}

## function forces a tree to be ultrametric using two different methods
## written by Liam J. Revell 2017

force.ultrametric<-function(tree,method=c("nnls","extend")){
	method<-method[1]
	if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
		rooted=TRUE,trace=0)
	else if(method=="extend"){
		h<-diag(vcv(tree))
		d<-max(h)-h
		ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
			y=tree$edge[,2])
		tree$edge.length[ii]<-tree$edge.length[ii]+d
	} else 
		cat("method not recognized: returning input tree\n\n")
	tree
}

## function to create curved clade labels for a fan tree
## written by Liam J. Revell 2017

arc.cladelabels<-function(tree=NULL,text,node=NULL,ln.offset=1.02,
	lab.offset=1.06,cex=1,orientation="curved",...){
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(obj$type!="fan") stop("method works only for type=\"fan\"")
	h<-max(sqrt(obj$xx^2+obj$yy^2))
	if(hasArg(mark.node)) mark.node<-list(...)$mark.node
	else mark.node<-TRUE
	if(hasArg(interactive)) interactive<-list(...)$interactive
	else {
		if(is.null(node)) interactive<-TRUE
		else interactive<-FALSE
	}
	if(interactive) node<-getnode()
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-par()$lwd
	if(hasArg(col)) col<-list(...)$col
	else col<-par()$col
	if(mark.node) points(obj$xx[node],obj$yy[node],pch=21,
		bg="red")
	if(is.null(tree)){
		tree<-list(edge=obj$edge,tip.label=1:obj$Ntip,
			Nnode=obj$Nnode)
		class(tree)<-"phylo"
	}
	d<-getDescendants(tree,node)
	d<-sort(d[d<=Ntip(tree)])
	deg<-atan(obj$yy[d]/obj$xx[d])*180/pi
	ii<-intersect(which(obj$yy[d]>=0),which(obj$xx[d]<0))
	deg[ii]<-180+deg[ii]
	ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]<0))
	deg[ii]<-180+deg[ii]
	ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]>=0))
	deg[ii]<-360+deg[ii]
	draw.arc(x=0,y=0,radius=ln.offset*h,deg1=min(deg),
		deg2=max(deg),lwd=lwd,col=col)
	if(orientation=="curved")
		arctext(text,radius=lab.offset*h,
			middle=mean(range(deg*pi/180)),cex=cex)
	else if(orientation=="horizontal"){
		x0<-lab.offset*cos(median(deg)*pi/180)*h
		y0<-lab.offset*sin(median(deg)*pi/180)*h
		text(x=x0,y=y0,label=text,
		adj=c(if(x0>=0) 0 else 1,if(y0>=0) 0 else 1),
		offset=0)
	}
}

## function to return a node index interactively from a plotted tree
## written by Liam J. Revell 2017
getnode<-function(...){
	if(hasArg(env)) env<-list(...)$env
	else env<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(hasArg(show.pt)) show.pt<-list(...)$show.pt
	else show.pt<-FALSE
	xy<-unlist(locator(n=1))
	if(show.pt) points(xy[1],xy[2])
	d<-sqrt((xy[1]-env$xx)^2+(xy[2]-env$yy)^2)
	ii<-which(d==min(d))[1]
	ii
}

## function mostly to interactively label nodes by clicking
## written by Liam J. Revell 2017
labelnodes<-function(text,node=NULL,interactive=TRUE,
	shape=c("circle","ellipse","rect"),...){
	shape<-shape[1]
	if(hasArg(circle.exp)) circle.exp<-list(...)$circle.exp
	else circle.exp<-1.3
	if(hasArg(rect.exp)) rect.exp<-list(...)$rect.exp
	else rect.exp<-1.6
	if(hasArg(cex)) cex<-list(...)$cex
	else cex<-1
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	h<-cex*strheight("A")
	w<-cex*strwidth(text)
	rad<-circle.exp*h*diff(par()$usr[1:2])/diff(par()$usr[3:4])
	if(is.null(node)){
		if(!interactive){
			cat("No nodes provided. Setting interactive mode to TRUE.\n")
			interactive<-TRUE
		}
		node<-vector(length=length(text))
	}
	for(i in 1:length(text)){
		if(interactive){
			cat(paste("Click on the node you would like to label ",
				text[i],".\n",sep=""))
			flush.console()
			ii<-getnode(env=obj)
			node[i]<-ii
		} else ii<-node[i]
		if(shape=="circle")
			draw.circle(obj$xx[ii],obj$yy[ii],rad,col="white")
		else if(shape=="ellipse")
			draw.ellipse(obj$xx[ii],obj$yy[ii],0.8*w[i],h,
				col="white")
		else if(shape=="rect")
			rect(xleft=obj$xx[ii]-0.5*rect.exp*w[i],
				ybottom=obj$yy[ii]-0.5*rect.exp*h,
				xright=obj$xx[ii]+0.5*rect.exp*w[i],
				ytop=obj$yy[ii]+0.5*rect.exp*h,col="white",
				ljoin=1)
		text(obj$xx[ii],obj$yy[ii],label=text[i],cex=cex)
	}
	invisible(node)
}

## convert object of class "birthdeath" into birth & death rates
bd<-function(x){
	if(class(x)!="birthdeath") stop("x should be an object of class 'birthdeath'")
	b<-x$para[2]/(1-x$para[1])
	d<-b-x$para[2]
	setNames(c(b,d),c("b","d"))
}

## compute AIC weights
aic.w<-function(aic){
	d.aic<-aic-min(aic)
	x<-exp(-1/2*d.aic)/sum(exp(-1/2*d.aic))
	class(x)<-"aic.w"
	x
}

print.aic.w<-function(x,...){
	if(hasArg(signif)) signif<-list(...)$signif
	else signif<-8
	print(round(unclass(x),signif))
}

## function to compute all paths towards the tips from a node
## written by  Liam J. Revell
node.paths<-function(tree,node){
	d<-Descendants(tree,node,"children")
	paths<-as.list(d)
	while(any(d>Ntip(tree))){
		jj<-1
		new.paths<-list()
		for(i in 1:length(paths)){
			if(paths[[i]][length(paths[[i]])]<=Ntip(tree)){ 
				new.paths[[jj]]<-paths[[i]]
				jj<-jj+1
			} else {
				ch<-Descendants(tree,paths[[i]][length(paths[[i]])],
					"children")
				for(j in 1:length(ch)){
					new.paths[[jj]]<-c(paths[[i]],ch[j])
					jj<-jj+1
				}
			}
		}
		paths<-new.paths
		d<-sapply(paths,function(x) x[length(x)])
	}
	paths
}

## function to compute a modification of Grafen's edge lengths
## written by Liam J. Revell 2016
modified.Grafen<-function(tree,power=2){
	max.np<-function(tree,node){
		np<-node.paths(tree,node)
		if(length(np)>0) max(sapply(np,length)) else 0
	}
	nn<-1:(Ntip(tree)+tree$Nnode)
	h<-sapply(nn,max.np,tree=tree)+1
	h<-(h/max(h))^power
	edge.length<-vector()
	for(i in 1:nrow(tree$edge)) 
		edge.length[i]<-diff(h[tree$edge[i,2:1]])
	tree$edge.length<-edge.length
	tree
}

## function to compute all rotations
## written by Liam J. Revell 2016
allRotations<-function(tree){
	if(!is.binary.tree(tree)){
		was.binary<-FALSE
		if(is.null(tree$edge.length)){ 
			tree<-compute.brlen(tree)
			had.edge.lengths<-FALSE
		} else had.edge.lengths<-TRUE
		tree<-multi2di(tree)
	} else was.binary<-TRUE
	nodes<-1:tree$Nnode+Ntip(tree)
	trees<-vector(mode="list",length=2^length(nodes))
	ii<-2
	trees[[1]]<-tree
	for(i in 1:length(nodes)){
		N<-ii-1
		for(j in 1:N){
			trees[[ii]]<-rotate(trees[[j]],nodes[i])
			ii<-ii+1
		}
	}
	trees<-lapply(trees,untangle,"read.tree")
	if(!was.binary){
		trees<-lapply(trees,di2multi)
		if(!had.edge.lengths) trees<-lapply(trees,
			function(x){
				x$edge.length<-NULL
				x
			})
	}
	class(trees)<-"multiPhylo"
	trees
}

## function to rotate a multifurcation in all possible ways
## written by Liam J. Revell 2016
rotate.multi<-function(tree,node){
	kids<-Children(tree,node)
	if(length(kids)>2){
		ii<-sapply(kids,function(x,y) which(y==x),y=tree$edge[,2])
		jj<-permn(ii)
		foo<-function(j,i,t){
			t$edge[i,]<-t$edge[j,]
			if(!is.null(t$edge.length))
				t$edge.length[i]<-t$edge.length[j]
			untangle(t,"read.tree")
		}
		obj<-lapply(jj[2:length(jj)],foo,i=ii,t=tree)
		class(obj)<-"multiPhylo"
	} else obj<-untangle(rotate(tree,node),"read.tree")
	obj
}

## wrapper for bind.tree that takes objects of class "simmap"
## written by Liam J. Revell 2016
bind.tree.simmap<-function(x,y,where="root"){
	x<-reorder(x)
	y<-reorder(y)
	rootx<-x$edge[1,1]
	rooty<-y$edge[1,1]
	xy<-read.tree(text=write.tree(bind.tree(x,y,where)))
	Mx<-rbind(matchLabels(x,xy),matchNodes(x,xy,"distances"))
	My<-rbind(matchLabels(y,xy),matchNodes(y,xy,"distances"))
	if(where!="root"&&where<=Ntip(x))
		Mx[which(is.na(Mx[,2])),2]<-findMRCA(xy,y$tip.label)
	xy$maps<-vector(mode="list",length=nrow(xy$edge))
	ix<-sapply(Mx[-which(Mx[,1]==rootx),1],
		function(x,y) which(y==x),y=x$edge[,2])
	ixy<-sapply(Mx[-which(Mx[,1]==rootx),2],
		function(x,y) which(y==x),y=xy$edge[,2])
	xy$maps[ixy]<-x$maps[ix]
	iy<-sapply(My[-which(My[,1]==rooty),1],
		function(x,y) which(y==x),y=y$edge[,2])
	ixy<-sapply(My[-which(My[,1]==rooty),2],
		function(x,y) which(y==x),y=xy$edge[,2])
	xy$maps[ixy]<-y$maps[iy]
	xy$mapped.edge<-makeMappedEdge(xy$edge,xy$maps)
	ns<-c(setNames(getStates(xy,"tips"),1:Ntip(xy)),
		getStates(xy,"nodes"))
	xy$node.states<-cbind(ns[as.character(xy$edge[,1])],
		ns[as.character(xy$edge[,2])])
	xy$states<-getStates(xy,"tips")
	attr(xy,"class")<-c("simmap",class(xy))
	xy
}

## generic function to convert an object of class "simmap" to "phylo"
## written by Liam J. Revell 2016
as.phylo.simmap<-function(x,...){
	x$maps<-NULL
	x$mapped.edge<-NULL
	if(!is.null(x$node.states)) x$node.states<-NULL
	if(!is.null(x$states)) x$states<-NULL
	if(!is.null(x$Q)) x$Q<-NULL
	if(!is.null(x$logL)) x$logL<-NULL
	if(!is.null(attr(x,"map.order"))) attr(x,"map.order")<-NULL
	class(x)<-setdiff(class(x),"simmap")
	x
}

## generic function to convert an object of class "multiSimmap" to "multiPhylo"
## written by Liam J. Revell 2016
as.multiPhylo.multiSimmap<-function(x,...){
	obj<-lapply(x,as.phylo)
	class(obj)<-setdiff(class(x),"multiSimmap")
	obj
}

## generic function to convert object of class "phylo" to "multiPhylo"
## written by Liam J. Revell 2016
as.multiPhylo.phylo<-function(x,...){
	obj<-list(x)
	class(obj)<-"multiPhylo"
	obj
}

as.multiPhylo<-function(x,...){
	if (identical(class(x),"multiPhylo")) return(x)
	UseMethod("as.multiPhylo")
}

## get mapped states
## written by Liam J. Revell 2016
mapped.states<-function(tree,...){
	if(!(inherits(tree,"simmap")||inherits(tree,"multiSimmap")))
		stop("tree should be an object of class \"simmap\" or \"multiSimmap\".")
	else {
		if(inherits(tree,"simmap")){
			if(!is.null(tree$mapped.edge)) 
				obj<-sort(colnames(tree$mapped.edge))
			else 
				obj<-sort(unique(unlist(lapply(tree$maps,function(x) names(x)))))
		} else if(inherits(tree,"multiSimmap")) {
			obj<-sapply(tree,mapped.states,...)
		}
	}
	obj
}

## match labels between trees (equivalent to matchNodes)
## written by Liam J. Revell 2016
matchLabels<-function(tr1,tr2){
	foo<-function(x,y) if(length(obj<-which(y==x))>0) obj else NA
	M<-cbind(1:Ntip(tr1),sapply(tr1$tip.label,foo,y=tr2$tip.label))
	colnames(M)<-c("tr1","tr2")
	M
}

## compute the probability of states changes along edges of the tree
## written by Liam J. Revell 2015
edgeProbs<-function(trees){
	if(!inherits(trees,"multiSimmap")) stop("trees should be an object of class \"multiSimmap\".")
	SS<-sapply(trees,getStates,"tips")
	states<-sort(unique(as.vector(SS)))
	m<-length(states)
	TT<-sapply(states,function(x,y) sapply(y,paste,x,sep="->"),
		y=states)
	nn<-c(TT[upper.tri(TT)],TT[lower.tri(TT)])
	## this function computes for a given edge
	fn<-function(edge,trees,states){
		obj<-sapply(trees,function(x,e,s) 
		if(names(x$maps[[e]])[1]==
			s[1]&&names(x$maps[[e]])[length(x$maps[[e]])]==s[2]) TRUE
		else FALSE,e=edge,s=states)
		sum(obj)/length(obj)
	}
	edge.probs<-matrix(0,nrow(trees[[1]]$edge),m,
		dimnames=list(apply(trees[[1]]$edge,1,paste,collapse=","),nn))
	k<-1
	for(i in 1:m) for(j in 1:m){
		if(i!=j){ 
			edge.probs[,k]<-sapply(1:nrow(trees[[1]]$edge),fn,
				trees=trees,states=states[c(i,j)])
			k<-k+1
		}
	}
	edge.probs<-cbind(edge.probs,1-rowSums(edge.probs))
	colnames(edge.probs)[ncol(edge.probs)]<-"no change"
	edge.probs
}

## get a position in the tree interactively
## written by Liam J. Revell 2015, 2016
get.treepos<-function(message=TRUE,...){
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(obj$type=="phylogram"&&obj$direction=="rightwards"){
		if(message){ 
			cat("Click on the tree position you want to capture...\n")
			flush.console()
		}
		if(hasArg(x)) x<-list(...)$x
		else x<-NULL
		if(hasArg(y)) y<-list(...)$y
		else y<-NULL
		if(is.null(x)||is.null(y)){
			x<-unlist(locator(1)) 	
			y<-x[2] 	
			x<-x[1]
		}
		d<-pos<-c()
		for(i in 1:nrow(obj$edge)){
			x0<-obj$xx[obj$edge[i,]]
			y0<-obj$yy[obj$edge[i,2]]
			if(x<x0[1]||x>x0[2]){
				d[i]<-min(dist(rbind(c(x,y),c(x0[1],y0))),
					dist(rbind(c(x,y),c(x0[2],y0))))
				pos[i]<-if(x>x0[2]) 0 else diff(obj$xx[obj$edge[i,]])
			} else {
				d[i]<-abs(y0-y)
				pos[i]<-obj$xx[obj$edge[i,2]]-x
			}
		}
		ii<-which(d==min(d))
		list(where=obj$edge[ii,2],pos=pos[ii])
	} else stop("Does not work for the plotted tree type.")
}

## fastDist: uses fastHeight to compute patristic distance between a pair of species
fastDist<-function(tree,sp1,sp2){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.null(tree$edge.length)) stop("tree should have edge lengths.")
	if(sp1==sp2) 0
	else fastHeight(tree,sp1,sp1)+fastHeight(tree,sp2,sp2)-
		2*fastHeight(tree,sp1,sp2)
}

# function reorders simmap tree
# written Liam Revell 2011, 2013, 2015
reorderSimmap<-function(tree,order="cladewise",index.only=FALSE,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	ii<-reorder.phylo(tree,order,index.only=TRUE,...)
	if(!index.only){
		if(inherits(ii,"phylo")) ii<-whichorder(ii$edge[,2],tree$edge[,2]) ## bug workaround
		tree$edge<-tree$edge[ii,]
		tree$edge.length<-tree$edge.length[ii]
		if(!is.null(tree$maps)){
			tree$maps<-tree$maps[ii]
			tree$mapped.edge<-tree$mapped.edge[ii,]
		}
		attr(tree,"order")<-order
		return(tree)
	} else return(ii)
}

## S3 reorder method for objects of class "simmap"
reorder.simmap<-function(x,...) reorderSimmap(x,...)

# function whichorder
# written by Liam Revell 2011, 2013, 2015
whichorder<-function(x,y) sapply(x,function(x,y) which(x==y),y=y)

# function returns random state with probability given by y
# written by Liam J. Revell 2013, 2015
rstate<-function(y){
	if(length(y)==1) return(names(y)[1])
	else {
		p<-y/sum(y)
		if(any(p<0)){ 
			warning("Some probabilities (slightly?) < 0. Setting p < 0 to zero.")
			p[p<0]<-0
		}
		return(names(which(rmultinom(1,1,p)[,1]==1)))
	}
}

## mark the changes on a plotted "simmap" object
## written by Liam J. Revell 2015
markChanges<-function(tree,colors=NULL,cex=1,lwd=2,plot=TRUE){
	states<-sort(unique(getStates(tree)))
	if(is.null(colors)) colors<-setNames(palette()[1:length(states)],
		states)
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	nc<-sapply(tree$maps,length)-1
	ii<-which(nc>0)
	nc<-nc[ii]
	xx<-yy<-vector()
	for(i in 1:length(ii)){
		for(j in 1:nc[i]){
			ss<-names(tree$maps[[ii[i]]])[j+1]
			mm<-tree$edge[ii[i],1]
			dd<-tree$edge[ii[i],2]
			x<-rep(obj$xx[mm]+cumsum(tree$maps[[ii[i]]])[j],2)
			y<-c(obj$yy[dd]-0.5*mean(strheight(LETTERS)*cex),
				obj$yy[dd]+0.5*mean(strheight(LETTERS)*cex))
			if(plot) lines(x,y,lwd=lwd,col=colors[ss],lend=2)
			xx<-c(xx,setNames(x[1],
				paste(names(tree$maps[[ii[i]]])[j:(j+1)],
				collapse="->")))
			yy<-c(yy,mean(y))
		}
	}
	XY<-cbind(xx,yy)
	colnames(XY)<-c("x","y")
	invisible(XY)
}

## function to label clades
## written by Liam J. Revell 2014, 2015
cladelabels<-function(tree=NULL,text,node,offset=NULL,wing.length=NULL,cex=1,
	orientation="vertical"){
	lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(is.null(tree)){
		wing.length<-1
		if(is.null(offset)) offset<-8
		tree<-list(edge=lastPP$edge,
			tip.label=1:lastPP$Ntip,
			Nnode=lastPP$Nnode)
		H<-matrix(lastPP$xx[tree$edge],nrow(tree$edge),2)
		tree$edge.length<-H[,2]-H[,1]
		class(tree)<-"phylo"
	}
	if(is.null(offset)) offset<-0.5
	xx<-mapply(labelSubTree,node,text,
		MoreArgs=list(tree=tree,pp=lastPP,offset=offset,wl=wing.length,cex=cex,
		orientation=orientation))
}

## internal function used by cladelabels
## written by Liam J. Revell 2014, 2015
labelSubTree<-function(tree,nn,label,pp,offset,wl,cex,orientation){
	if(is.null(wl)) wl<-1
	tree<-reorder(tree)
	tips<-getDescendants(tree,nn)
	tips<-tips[tips<=length(tree$tip.label)]
	ec<-0.7 ## expansion constant
	sw<-pp$cex*max(strwidth(tree$tip.label[tips]))
	sh<-pp$cex*max(strheight(tree$tip.label))
	cw<-mean(strwidth(LETTERS)*cex)	
	h<-max(sapply(tips,function(x,tree)
		nodeHeights(tree)[which(tree$edge[,2]==x),2],
		tree=tree))+sw+offset*cw
	y<-range(pp$yy[tips])
	lines(c(h,h),y+ec*c(-sh,sh))
	lines(c(h-wl*cw,h),
		c(y[1]-ec*sh,y[1]-ec*sh))
	lines(c(h-wl*cw,h),
		c(y[2]+ec*sh,y[2]+ec*sh))
	text(h+cw,mean(y),
		label,srt=if(orientation=="horizontal") 0 else 90,
		adj=if(orientation=="horizontal") 0 else 0.5,cex=cex)
}

## get all the extant/extinct tip names
## written by Liam J. Revell 2012, 2015

getExtant<-function(tree,tol=1e-8){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	H<-nodeHeights(tree)
	tl<-max(H)
	x<-which(H[,2]>=(tl-tol))
	y<-tree$edge[x,2]
	y<-y[y<=length(tree$tip)]
	z<-tree$tip.label[y]
	return(z)
}

getExtinct<-function(tree,tol=1e-8) setdiff(tree$tip.label,getExtant(tree,tol))

# function splits tree at split
# written by Liam Revell 2011, 2014, 2015

splitTree<-function(tree,split){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(split$node>length(tree$tip.label)){
		# first extract the clade given by shift$node
		tr2<-extract.clade(tree,node=split$node)
		tr2$root.edge<-tree$edge.length[which(tree$edge[,2]==split$node)]-split$bp
		#now remove tips in tr2 from tree
		tr1<-drop.clade(tree,tr2$tip.label)
		nn<-if(!is.null(tree$node.label)) c(tree$node.label,"NA") else "NA"
		tr1$tip.label[which(tr1$tip.label%in%nn)]<-"NA"
		tr1$edge.length[match(which(tr1$tip.label=="NA"),tr1$edge[,2])]<-split$bp
	} else {
		# first extract the clade given by shift$node
		tr2<-list(edge=matrix(c(2L,1L),1,2),tip.label=tree$tip.label[split$node],edge.length=tree$edge.length[which(tree$edge[,2]==split$node)]-split$bp,Nnode=1L)
		class(tr2)<-"phylo"
		# now remove tip in tr2 from tree
		tr1<-tree
		tr1$edge.length[match(which(tr1$tip.label==tr2$tip.label[1]),tr1$edge[,2])]<-split$bp
		tr1$tip.label[which(tr1$tip.label==tr2$tip.label[1])]<-"NA"
	}
	trees<-list(tr1,tr2)
	class(trees)<-"multiPhylo"
	trees
}


# function drops entire clade
# written by Liam Revell 2011, 2015

drop.clade<-function(tree,tip){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	nn<-if(!is.null(tree$node.label)) c(tree$node.label,"NA") else "NA"
	tree<-drop.tip(tree,tip,trim.internal=FALSE)
	while(sum(tree$tip.label%in%nn)>1)
		tree<-drop.tip(tree,tree$tip.label[tree$tip.label%in%nn],
			trim.internal=FALSE)
	tree
}


## function to re-root a phylogeny along an edge
## written by Liam J. Revell 2011-2016

reroot<-function(tree,node.number,position=NULL,interactive=FALSE,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(interactive){
		plotTree(tree,...)
		cat("Click where you would like re-root the plotted tree\n")
		flush.console()
		obj<-get.treepos(message=FALSE)
		node.number<-obj$where
		position<-tree$edge.length[which(tree$edge[,2]==node.number)]-obj$pos
	}
	if(is.null(position)) position<-tree$edge.length[which(tree$edge[,2]==node.number)]
	tt<-splitTree(tree,list(node=node.number,bp=position))
	p<-tt[[1]]
	d<-tt[[2]]
	tip<-if(length(which(p$tip.label=="NA"))>0) "NA" else p$tip.label[which(p$tip.label%in%tree$node.label)]
	p<-root(p,outgroup=tip,resolve.root=T)
	bb<-which(p$tip.label==tip)
	p$tip.label[bb]<-"NA"
	ee<-p$edge.length[which(p$edge[,2]==bb)]
	p$edge.length[which(p$edge[,2]==bb)]<-0
	cc<-p$edge[which(p$edge[,2]==bb),1]
	dd<-setdiff(p$edge[which(p$edge[,1]==cc),2],bb)
	p$edge.length[which(p$edge[,2]==dd)]<-p$edge.length[which(p$edge[,2]==dd)]+ee
	obj<-paste.tree(p,d)
	if(interactive) plotTree(obj,...)
	obj
}

## function to add an arrow pointing to a tip or node in the tree
## written by Liam J. Revell 2014, 2017

add.arrow<-function(tree=NULL,tip,...){
	lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(!is.null(tree)){
		if(inherits(tree,"contMap")) tree<-tree$tree
		else if(inherits(tree,"densityMap")) tree<-tree$tree
	}
	if(is.numeric(tip)){
		ii<-tip
		if(!is.null(tree)&&ii<=length(tree$tip.label)) tip<-tree$tip.label[ii]
		else tip<-""
	} else if(is.character(tip)&&!is.null(tree)) ii<-which(tree$tip.label==tip)
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-1
	strw<-lastPP$cex*(strwidth(tip)+offset*mean(strwidth(c(LETTERS,letters))))
	if(hasArg(arrl)) arrl<-list(...)$arrl
	else { 
		if(lastPP$type=="fan") arrl<-0.3*max(lastPP$xx)
		else if(lastPP$type=="phylogram") arrl<-0.15*max(lastPP$xx)
	}
	if(hasArg(hedl)) hedl<-list(...)$hedl
	else hedl<-arrl/3
	if(hasArg(angle)) angle<-list(...)$angle
	else angle<-45
	arra<-angle*pi/180
	asp<-if(lastPP$type=="fan") 1 else (par()$usr[4]-par()$usr[3])/(par()$usr[2]-par()$usr[1])
	if(hasArg(col)) col<-list(...)$col
	else col<-"black"
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(lastPP$type=="fan") theta<-atan2(lastPP$yy[ii],lastPP$xx[ii])
	else if(lastPP$type=="phylogram") theta<-0	
	segments(x0=lastPP$xx[ii]+cos(theta)*(strw+arrl),
		y0=lastPP$yy[ii]+sin(theta)*(strw+arrl),
		x1=lastPP$xx[ii]+cos(theta)*strw,
		y1=lastPP$yy[ii]+sin(theta)*strw,
		col=col,lwd=lwd,lend="round")
	segments(x0=lastPP$xx[ii]+cos(theta)*strw+cos(theta+arra/2)*hedl,
		y0=lastPP$yy[ii]+sin(theta)*strw+sin(theta+arra/2)*hedl*asp,
		x1=lastPP$xx[ii]+cos(theta)*strw,
		y1=lastPP$yy[ii]+sin(theta)*strw,
		col=col,lwd=lwd,lend="round")
	segments(x0=lastPP$xx[ii]+cos(theta)*strw+cos(theta-arra/2)*hedl,
		y0=lastPP$yy[ii]+sin(theta)*strw+sin(theta-arra/2)*hedl*asp,
		x1=lastPP$xx[ii]+cos(theta)*strw,
		y1=lastPP$yy[ii]+sin(theta)*strw,
		col=col,lwd=lwd,lend="round")
	invisible(list(x0=lastPP$xx[ii]+cos(theta)*(strw+arrl),
		y0=lastPP$yy[ii]+sin(theta)*(strw+arrl),
		x1=lastPP$xx[ii]+cos(theta)*strw,
		y1=lastPP$yy[ii]+sin(theta)*strw))
}

## function to ladderize phylogeny with mapped discrete character
## written by Liam J. Revell 2014, 2015

ladderize.simmap<-function(tree,right=TRUE){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	obj<-read.tree(text=write.tree(ladderize(tree,right=right)))
	rN<-Ntip(obj)+1
	T<-cbind(1:Ntip(obj),sapply(obj$tip.label,function(x,y) which(y==x),y=tree$tip.label))
	N<-matchNodes(obj,tree)
	M<-rbind(T,N[N[,1]!=rN,])
	ii<-sapply(M[,1],function(x,y) which(y==x),y=obj$edge[,2])
	jj<-sapply(M[,2],function(x,y) which(y==x),y=tree$edge[,2])
	obj$maps<-vector(length=nrow(tree$edge),mode="list")
	obj$mapped.edge<-matrix(NA,nrow(tree$edge),ncol(tree$mapped.edge),
		dimnames=list(apply(tree$edge,1,paste,collapse=","),
		colnames(tree$mapped.edge)))
	if(!is.null(tree$states)) 
		obj$states<-tree$states[sapply(obj$tip.label,function(x,y) which(y==x),y=tree$tip.label)]
	if(!is.null(tree$node.states)) obj$node.states<-matrix(NA,nrow(tree$edge),2)
	for(i in 1:length(ii)){
		obj$maps[[ii[i]]]<-tree$maps[[jj[i]]]
		obj$mapped.edge[ii[i],]<-tree$mapped.edge[jj[i],]
		if(!is.null(tree$node.states)) obj$node.states[ii[i],]<-tree$node.states[jj[i],]
	}
	obj
}

## for backward compatibility

repPhylo<-function(tree,times) rep(tree,times)

## S3 method rep for objects of class "phylo" and "multiPhylo"
## written by Liam J. Revell 2014

rep.phylo<-function(x,...){
	if(hasArg(times)) times<-list(...)$times
	else times<-(...)[[1]]
	for(i in 1:times) obj<-if(i==1) x else if(i==2) c(obj,x) else c(obj,list(x))
	class(obj)<-"multiPhylo"
	obj
}

rep.multiPhylo<-function(x,...){
	if(hasArg(times)) times<-list(...)$times
	else times<-(...)[[1]]
	for(i in 1:times) obj<-if(i==1) x else if(i>=2) c(obj,x)
	class(obj)<-"multiPhylo"
	obj
}

## function to rescale simmap style trees
## written by Liam J. Revell 2012, 2013, 2014, 2015, 2017
rescaleSimmap<-function(tree,...){
	if(inherits(tree,"multiPhylo")){
		cls<-class(tree)
		trees<-unclass(tree)
		trees<-lapply(trees,rescaleSimmap,...)
		class(trees)<-cls
		return(trees)
	} else if(inherits(tree,"phylo")){
		if(hasArg(lambda)) lambda<-list(...)$lambda
		else lambda<-1
		if(hasArg(totalDepth)) depth<-totalDepth<-list(...)$totalDepth
		else if(hasArg(depth)) depth<-totalDepth<-list(...)$depth
		else depth<-totalDepth<-max(nodeHeights(tree))
		if(lambda!=1){
			e<-lambdaTree(tree,lambda)$edge.length/tree$edge.length
			tree$edge.length<-tree$edge.length*e
			tree$maps<-mapply(function(x,y) x*y,tree$maps,e)
			tree$mapped.edge<-tree$mapped.edge*matrix(rep(e,ncol(tree$mapped.edge)),length(e),ncol(tree$mapped.edge))
		}
		if(depth!=max(nodeHeights(tree))){
			h<-max(nodeHeights(tree))
			s<-depth/h
			tree$edge.length<-tree$edge.length*s
			tree$maps<-lapply(tree$maps,"*",s)
			tree$mapped.edge<-tree$mapped.edge*s
		}
		return(tree)
	} else message("tree should be an object of class \"phylo\" or \"multiPhylo\"")
}

## function to drop one or more tips from a tree but retain all ancestral nodes as singletons
## written by Liam J. Revell 2014, 2015
drop.tip.singleton<-function(tree,tip){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	N<-Ntip(tree)
	m<-length(tip)
	ii<-sapply(tip,function(x,y) which(y==x),y=tree$tip.label)
	tree$tip.label<-tree$tip.label[-ii]
	ii<-sapply(ii,function(x,y) which(y==x),y=tree$edge[,2])
	tree$edge<-tree$edge[-ii,]
	tree$edge.length<-tree$edge.length[-ii]
	tree$edge[tree$edge<=N]<-as.integer(rank(tree$edge[tree$edge<=N]))
	tree$edge[tree$edge>N]<-tree$edge[tree$edge>N]-m
	N<-N-m
	if(any(sapply(tree$edge[tree$edge[,2]>N,2],"%in%",tree$edge[,1])==FALSE)) internal<-TRUE
	else internal<-FALSE
	while(internal){
		ii<-which(sapply(tree$edge[,2],"%in%",c(1:N,tree$edge[,1]))==FALSE)[1]
		nn<-tree$edge[ii,2]
		tree$edge<-tree$edge[-ii,]
		tree$edge.length<-tree$edge.length[-ii]
		tree$edge[tree$edge>nn]<-tree$edge[tree$edge>nn]-1
		tree$Nnode<-tree$Nnode-length(ii)
		if(any(sapply(tree$edge[tree$edge[,2]>N,2],"%in%",tree$edge[,1])==FALSE)) internal<-TRUE
		else internal<-FALSE
	}
	tree
}

## S3 print method for object of class 'describe.simmap'
## written by Liam J. Revell 2014, 2015
print.describe.simmap<-function(x,...){
	if(inherits(x$tree,"multiPhylo")){
		cat(paste(length(x$tree),"trees with a mapped discrete character with states:\n",paste(colnames(x$ace),collapse=", "),"\n\n"))
		cat(paste("trees have",colMeans(x$count)["N"],"changes between states on average\n\n"))
		cat(paste("changes are of the following types:\n"))
		aa<-t(as.matrix(colMeans(x$count)[2:ncol(x$count)]))
		rownames(aa)<-"x->y"
		print(aa)
		cat(paste("\nmean total time spent in each state is:\n"))
		print(matrix(c(colMeans(x$times),colMeans(x$times[,1:ncol(x$times)]/x$times[,ncol(x$times)])),2,ncol(x$times),byrow=TRUE,
			dimnames=list(c("raw","prop"),c(colnames(x$times)))))
		cat("\n")
	} else if(inherits(x$tree,"phylo")){
		cat(paste("1 tree with a mapped discrete character with states:\n",paste(colnames(x$Tr),collapse=", "),"\n\n"))
		cat(paste("tree has",x$N,"changes between states\n\n"))
		cat(paste("changes are of the following types:\n"))
		print(x$Tr)
		cat(paste("\nmean total time spent in each state is:\n"))
		print(x$times)
		cat("\n")
	}
}

## S3 plot method for object of class 'describe.simmap'
## written by Liam J. Revell 2014, 2015
plot.describe.simmap<-function(x,...){
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(cex)) cex<-list(...)$cex
	else cex<-c(0.6,0.4)
	if(length(cex)==1) cex<-rep(cex,2)
	if(hasArg(type)) type<-list(...)$type
	else type<-"phylogram"
	if(inherits(x$tree,"multiPhylo")){
		states<-colnames(x$ace)
		if(hasArg(colors)) colors<-list(...)$colors
		else colors<-setNames(palette()[1:length(states)],states)
		plotTree(if(is.null(x$ref.tree)) x$tree[[1]] else x$ref.tree,lwd=lwd,
			offset=cex[2],...)
		nodelabels(pie=x$ace,piecol=colors[colnames(x$ace)],cex=cex[1])
		if(!is.null(x$tips)) tips<-x$tips else tips<-to.matrix(getStates(x$tree[[1]],"tips"),
			seq=states) 
		tiplabels(pie=tips[if(is.null(x$ref.tree)) x$tree[[1]]$tip.label else 
			x$ref.tree$tip.label,],piecol=colors[colnames(tips)],cex=cex[2])
	} else if(inherits(x$tree,"phylo")){
		states<-colnames(x$Tr)
		if(hasArg(colors)) colors<-list(...)$colors
		else colors<-setNames(palette()[1:length(states)],states)
		plotSimmap(x$tree,lwd=lwd,colors=colors,type=type)
	}
}


# function to summarize the results of stochastic mapping
# written by Liam J. Revell 2013, 2014, 2015
describe.simmap<-function(tree,...){
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-FALSE
	if(hasArg(check.equal)) check.equal<-list(...)$check.equal
	else check.equal<-FALSE
	if(hasArg(message)) message<-list(...)$message
	else message<-FALSE
	if(hasArg(ref.tree)) ref.tree<-list(...)$ref.tree
	else ref.tree<-NULL
	if(inherits(tree,"multiPhylo")){
		if(check.equal){
			TT<-sapply(tree,function(x,y) sapply(y,all.equal.phylo,x),y=tree)
			check<-all(TT)
			if(!check) cat("Note: Some trees are not equal.\nA \"reference\" tree will be computed if none was provided.\n\n")
		} else check<-TRUE
		YY<-getStates(tree)
		states<-sort(unique(as.vector(YY)))
		if(is.null(ref.tree)&&check) ZZ<-t(apply(YY,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=states,Nsim=length(tree)))
		else {
			if(is.null(ref.tree)){
				cat("No reference tree provided & some trees are unequal.\nComputing majority-rule consensus tree.\n")
				ref.tree<-consensus(tree,p=0.5)
			}
			YYp<-matrix(NA,ref.tree$Nnode,length(tree),dimnames=list(1:ref.tree$Nnode+Ntip(ref.tree),NULL))
			for(i in 1:length(tree)){
				M<-matchNodes(ref.tree,tree[[i]])
				jj<-sapply(M[,2],function(x,y) if(x%in%y) which(as.numeric(y)==x) else NA,y=as.numeric(rownames(YY)))
				YYp[,i]<-YY[jj,i]
			}
			ZZ<-t(apply(YYp,1,function(x,levels) summary(factor(x[!is.na(x)],levels))/sum(!is.na(x)),levels=states))
		}
		XX<-countSimmap(tree,states,FALSE)
		XX<-XX[,-(which(as.vector(diag(-1,length(states)))==-1)+1)]
		AA<-t(sapply(unclass(tree),function(x) c(colSums(x$mapped.edge),sum(x$edge.length))))
		colnames(AA)[ncol(AA)]<-"total"
		BB<-getStates(tree,type="tips")
		CC<-t(apply(BB,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=states,Nsim=length(tree)))
		x<-list(count=XX,times=AA,ace=ZZ,tips=CC,tree=tree,ref.tree=if(!is.null(ref.tree)) ref.tree else NULL)
		class(x)<-"describe.simmap"
	} else if(inherits(tree,"phylo")){
		XX<-countSimmap(tree,message=FALSE)
		YY<-getStates(tree)
		states<-sort(unique(YY))
		AA<-setNames(c(colSums(tree$mapped.edge),sum(tree$edge.length)),c(colnames(tree$mapped.edge),"total"))
		AA<-rbind(AA,AA/AA[length(AA)]); rownames(AA)<-c("raw","prop")
		x<-list(N=XX$N,Tr=XX$Tr,times=AA,states=YY,tree=tree)
		class(x)<-"describe.simmap"
	}
	if(message) print(x)
	if(plot) plot(x)
	x
}

## function finds the height of a given node
## written by Liam Revell 2014, 2015, 2016
nodeheight<-function(tree,node,...){
	if(hasArg(root.edge)) root.edge<-list(...)$root.edge
	else root.edge<-FALSE
	if(root.edge) ROOT<-if(!is.null(tree$root.edge)) tree$root.edge else 0
	else ROOT<-0 
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(node==(length(tree$tip.label)+1)) h<-0
	else {
		a<-setdiff(c(getAncestors(tree,node),node),length(tree$tip.label)+1)
		h<-sum(tree$edge.length[sapply(a,function(x,e) which(e==x),e=tree$edge[,2])])
	}
	h+ROOT
}

# fast pairwise MRCA function
# written by Liam Revell 2012, 2014, 2015
fastMRCA<-function(tree,sp1,sp2){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	x<-match(sp1,tree$tip.label)
	y<-match(sp2,tree$tip.label)
	a<-c(x,getAncestors(tree,x))
	b<-c(y,getAncestors(tree,y))
	z<-a%in%b
	return(a[min(which(z))])
}

## function to find the height of the MRCA of sp1 & sp2
## written by Liam J. Revell 2014, 2015
fastHeight<-function(tree,sp1,sp2){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.null(tree$edge.length)) stop("tree should have edge lengths.")
	sp1<-which(tree$tip.label==sp1)
	sp2<-which(tree$tip.label==sp2)
	a1<-c(sp1,getAncestors(tree,sp1))
	a2<-c(sp2,getAncestors(tree,sp2))
	a12<-intersect(a1,a2)
	if(length(a12)>1){ 
		a12<-a12[2:length(a12)-1]
		h<-sapply(a12,function(i,tree) tree$edge.length[which(tree$edge[,2]==i)],tree=tree)
		return(sum(h))
	} else return(0)
}

## function gets ancestor node numbers, to be used internally by 
## written by Liam J. Revell 2014
getAncestors<-function(tree,node,type=c("all","parent")){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	type<-type[1]
	if(type=="all"){
		aa<-vector()
		rt<-length(tree$tip.label)+1
		currnode<-node
		while(currnode!=rt){
			currnode<-getAncestors(tree,currnode,"parent")
			aa<-c(aa,currnode)
		}
		return(aa)
	} else if(type=="parent"){
		aa<-tree$edge[which(tree$edge[,2]==node),1]
		return(aa)
	} else stop("do not recognize type")
}

## function for midpoint rooting
## written by Liam J. Revell 2014
midpoint.root<-function(tree){
	D<-cophenetic(tree)
	dd<-max(D)
	ii<-which(D==dd)[1]
	ii<-c(ceiling(ii/nrow(D)),ii%%nrow(D))
	if(ii[2]==0) ii[2]<-nrow(D)
	spp<-rownames(D)[ii]
	nn<-which(tree$tip.label==spp[2])
	tree<-reroot(tree,nn,tree$edge.length[which(tree$edge[,2]==nn)])
	aa<-getAncestors(tree,which(tree$tip.label==spp[1]))
	D<-c(0,dist.nodes(tree)[which(tree$tip.label==spp[1]),aa])
	names(D)[1]<-which(tree$tip.label==spp[1])
	i<-0
	while(D[i+1]<(dd/2)) i<-i+1
	tree<-reroot(tree,as.numeric(names(D)[i]),D[i+1]-dd/2)
	tree
}

# function computes phylogenetic variance-covariance matrix, including for internal nodes
# written by Liam J. Revell 2011, 2013, 2014, 2015
vcvPhylo<-function(tree,anc.nodes=TRUE,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	# get & set optional arguments
	if(hasArg(internal)) internal<-list(...)$internal
	else internal<-anc.nodes
	if(internal!=anc.nodes){ 
		message(paste("arguments \"internal\" and \"anc.nodes\" are synonyms; setting internal =",anc.nodes))
		internal<-anc.nodes
	}
	if(hasArg(model)) model<-list(...)$model
	else model<-"BM"

	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-12
	if(model=="OU"){
		if(hasArg(alpha)) alpha<-list(...)$alpha
		else alpha<-0
	}
	if(model=="OU"&&alpha<tol) model<-"BM"
	if(model=="lambda"){
		if(hasArg(lambda)){ 
			lambda<-list(...)$lambda
			tree<-lambdaTree(tree,lambda)
		} else model<-"BM"
		model<-"BM"
	}	
	# done settings
	n<-length(tree$tip.label)
	h<-nodeHeights(tree)[order(tree$edge[,2]),2]
	h<-c(h[1:n],0,h[(n+1):length(h)])
	M<-mrca(tree,full=internal)[c(1:n,internal*(n+2:tree$Nnode)),c(1:n,internal*(n+2:tree$Nnode))]
	C<-matrix(h[M],nrow(M),ncol(M))
	if(internal) rownames(C)<-colnames(C)<-c(tree$tip.label,n+2:tree$Nnode)
	else rownames(C)<-colnames(C)<-tree$tip.label
	if(model=="OU"){
		D<-dist.nodes(tree)
		rownames(D)[1:n]<-colnames(D)[1:n]<-tree$tip.label
		D<-D[rownames(C),colnames(C)]
		# not sure 
		C<-(1-exp(-2*alpha*C))*exp(-alpha*D)/(2*alpha) # Hansen (2007)
		# C<-(1-exp(-2*alpha*C))*exp(-2*alpha*(1-C))/(2*alpha) # Butler & King (2004)
	}
	return(C)
}

## simplified lambdaTree to be used internally by vcvPhylo
## written by Liam J. Revell 2014
lambdaTree<-function(tree,lambda){
	ii<-which(tree$edge[,2]>length(tree$tip.label))
	H1<-nodeHeights(tree)
	tree$edge.length[ii]<-lambda*tree$edge.length[ii]
	H2<-nodeHeights(tree)
	tree$edge.length[-ii]<-tree$edge.length[-ii]+     H1[-ii,2]-H2[-ii,2]
	tree
}

## di2multi method for tree with mapped state
## written by Liam J. Revell 2013, 2015, 2016
di2multi.simmap<-function(phy,...){
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-08
	if(!inherits(phy,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.null(phy$maps)){
		cat("Warning: tree does not contain mapped state. Using di2multi.\n")
		return(di2multi(phy,tol))
	}
	N<-length(phy$tip.label)
	n<-length(intersect(which(phy$edge.length<tol),which(phy$edge[,2]>N)))
	if(n==0) return(phy)
	edge<-phy$edge
	edge[edge>N]<--edge[edge>N]+N
	edge.length<-phy$edge.length
	maps<-phy$maps
	Nnode<-phy$Nnode
	for(i in 1:n){
		ii<-intersect(which(edge.length<tol),which(edge[,2]<0))[1]
		node<-edge[ii,2]
		edge[edge==node]<-edge[ii,1]
		jj<-which(apply(edge,1,function(x) x[1]==x[2]))[1]
		edge<-edge[-jj,]
		edge.length<-edge.length[-jj]
		maps<-maps[-jj]	
		Nnode<-Nnode-1
	}
	nn<-sort(unique(edge[edge<0]),decreasing=TRUE)
	mm<-1:Nnode+N
	for(i in 1:length(edge)) if(edge[i]%in%nn) edge[i]<-mm[which(nn==edge[i])]
	mapped.edge<-makeMappedEdge(edge,maps)
	tt<-list(edge=edge,Nnode=Nnode,tip.label=phy$tip.label,edge.length=edge.length,
		maps=maps,mapped.edge=mapped.edge)
	class(tt)<-"phylo"
	if(!is.null(attr(phy,"order"))) attr(tt,"order")<-attr(phy,"order")
	if(!is.null(phy$node.states)) tt$node.states<-getStates(tt,"nodes")
	if(!is.null(phy$states)) tt$states<-getStates(tt,"tips")
	return(tt)
}

# returns the heights of each node
# written by Liam J. Revell 2011, 2012, 2013, 2015, 2016
nodeHeights<-function(tree,...){
	if(hasArg(root.edge)) root.edge<-list(...)$root.edge
	else root.edge<-FALSE
	if(root.edge) ROOT<-if(!is.null(tree$root.edge)) tree$root.edge else 0
	else ROOT<-0 
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(attr(tree,"order")!="cladewise"||is.null(attr(tree,"order"))) t<-reorder(tree)
	else t<-tree
	root<-length(t$tip.label)+1
	X<-matrix(NA,nrow(t$edge),2)
	for(i in 1:nrow(t$edge)){
		if(t$edge[i,1]==root){
			X[i,1]<-0.0
			X[i,2]<-t$edge.length[i]
		} else {
			X[i,1]<-X[match(t$edge[i,1],t$edge[,2]),2]
			X[i,2]<-X[i,1]+t$edge.length[i]
		}
	}
	if(attr(tree,"order")!="cladewise"||is.null(attr(tree,"order")))
		o<-apply(matrix(tree$edge[,2]),1,function(x,y) which(x==y),y=t$edge[,2])
	else o<-1:nrow(t$edge)
	return(X[o,]+ROOT)
}

## function drops all the leaves from the tree & collapses singleton nodes
## written by Liam J. Revell 2013, 2014, 2015
drop.leaves<-function(tree,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	## optional arguments
	if(hasArg(keep.tip.labels)) keep.tip.labels<-list(...)$keep.tip.labels
	else keep.tip.labels<-FALSE
	## end optional arguments
	n<-length(tree$tip.label)
	edge<-tree$edge
	edge[edge>n]<--edge[edge>n]+n
	ii<-which(edge[,2]>0)
	edge<-edge[-ii,]
	if(!is.null(tree$edge.length)){
		edge.length<-tree$edge.length
		edge.length<-edge.length[-ii]
	}
	zz<-sapply(edge[,2],function(x,y) !(x%in%y),y=edge[,1])
	if(is.null(tree$node.label)) tree$node.label<-1:tree$Nnode+n
	nn<-matrix(tree$node.label[-edge],nrow(edge),ncol(edge))
	tip.label<-nn[zz,2]
	node.label<-c(nn[1,1],nn[!zz,2])
	edge[zz,2]<-1:sum(zz)
	Nnode<-length(unique(edge[edge<0]))
	rr<-cbind(sort(unique(edge[edge<0]),decreasing=TRUE),1:Nnode+sum(zz))
	for(i in 1:nrow(rr)) edge[edge==rr[i,1]]<-rr[i,2]
	tt<-list(edge=edge,Nnode=Nnode,tip.label=tip.label,edge.length=edge.length,node.label=node.label)
	class(tt)<-"phylo"
	tt<-collapse.singles(tt)
	if(keep.tip.labels){
		for(i in 1:length(tt$tip.label)){
			yy<-getDescendants(tree,node=which(tree$node.label==tt$tip.label[i])+n)
			tt$tip.label[i]<-paste(tree$tip.label[yy[yy<=n]],collapse=",")
		}
	}
	return(tt)
}

# function rounds the branch lengths of the tree & applies rounding to simmap tree
# written by Liam J. Revell 2012, 2013, 2015
roundBranches<-function(tree,digits=0){
	if(inherits(tree,"multiPhylo")){
		trees<-lapply(tree,roundBranches,digits=digits)
		class(trees)<-"multiPhylo"
		return(trees)
	} else if(inherits(tree,"phylo")) {
		tree$edge.length<-round(tree$edge.length,digits)
		if(!is.null(tree$maps)){
			for(i in 1:nrow(tree$edge)){
				temp<-tree$maps[[i]]/sum(tree$maps[[i]])
				tree$maps[[i]]<-temp*tree$edge.length[i]
			}
		}
		if(!is.null(tree$mapped.edge)){
			a<-vector()
			for(i in 1:nrow(tree$edge)) a<-c(a,names(tree$maps[[i]]))
			a<-unique(a)
			tree$mapped.edge<-matrix(data=0,length(tree$edge.length),length(a),dimnames=list(apply(tree$edge,1,function(x) paste(x,collapse=",")),state=a))
			for(i in 1:length(tree$maps)) for(j in 1:length(tree$maps[[i]])) tree$mapped.edge[i,names(tree$maps[[i]])[j]]<-tree$mapped.edge[i,names(tree$maps[[i]])[j]]+tree$maps[[i]][j]
		}
		return(tree)
	} else stop("tree should be an object of class \"phylo\" or \"multiPhylo\".")
}

# function to merge mapped states
# written by Liam J. Revell 2013, 2015
mergeMappedStates<-function(tree,old.states,new.state){
	if(inherits(tree,"multiPhylo")){
		tree<-unclass(tree)
		tree<-lapply(tree,mergeMappedStates,old.states=old.states,new.state=new.state)
		class(tree)<-"multiPhylo"
	} else if(inherits(tree,"phylo")) {
		maps<-tree$maps
		rr<-function(map,oo,nn){ 
			for(i in 1:length(map)) if(names(map)[i]%in%oo) names(map)[i]<-nn
			map
		}
		mm<-function(map){
			if(length(map)>1){
				new.map<-vector()
				j<-1
				new.map[j]<-map[1]
				names(new.map)[j]<-names(map)[1]
				for(i in 2:length(map)){
					if(names(map)[i]==names(map)[i-1]){ 
						new.map[j]<-map[i]+new.map[j]
						names(new.map)[j]<-names(map)[i]
					} else {
						j<-j+1
						new.map[j]<-map[i]
						names(new.map)[j]<-names(map)[i]
					}

				}
				map<-new.map
			}
			map
		}
		maps<-lapply(maps,rr,oo=old.states,nn=new.state)
		if(length(old.states)>1){ 
			maps<-lapply(maps,mm)
			mapped.edge<-tree$mapped.edge
			mapped.edge<-cbind(rowSums(mapped.edge[,colnames(mapped.edge)%in%old.states]),
				mapped.edge[,setdiff(colnames(mapped.edge),old.states)])
			colnames(mapped.edge)<-c(new.state,setdiff(colnames(tree$mapped.edge),old.states))
		} else {
			mapped.edge<-tree$mapped.edge
			colnames(mapped.edge)[which(colnames(mapped.edge)==old.states)]<-new.state
		}
		tree$maps<-maps
		tree$mapped.edge<-mapped.edge
	} else stop("tree should be an object of class \"phylo\" or \"multiPhylo\".")
	return(tree)
}

# function rotates a node or multiple nodes
# written by Liam J. Revell 2013, 2015
rotateNodes<-function(tree,nodes,polytom=c(1,2),...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	n<-length(tree$tip.label)
	if(nodes[1]=="all") nodes<-1:tree$Nnode+n
	for(i in 1:length(nodes)) tree<-rotate(tree,nodes[i],polytom)
	if(hasArg(reversible)) reversible<-list(...)$reversible
	else reversible<-TRUE
	if(reversible){ 
		ii<-which(tree$edge[,2]<=n)
		jj<-tree$edge[ii,2]
		tree$edge[ii,2]<-1:n
		tree$tip.label<-tree$tip.label[jj]
	}
	return(tree)
}

# function simulates random sampling from xbar, xvar, with sample sizes n
# written by Liam J. Revell 2012
sampleFrom<-function(xbar=0,xvar=1,n=1,randn=NULL,type="norm"){
	if(length(xvar)==1&&length(xbar)!=length(xvar)) xvar<-rep(xvar,length(xbar))
	if(!is.null(randn))
		for(i in 1:length(xbar)) n[i]<-floor(runif(n=1,min=randn[1],max=(randn[2]+1)))
	x<-vector()

	for(i in 1:length(xbar)){
		y<-rnorm(n=n[i],mean=xbar[i],sd=sqrt(xvar[i]))
   		names(y)<-rep(names(xbar)[i],length(y))
   		x<-c(x,y)
	}
	return(x)
}

# function adds a new tip to the tree
# written by Liam J. Revell 2012, 2013, 2014, 2015
bind.tip<-function(tree,tip.label,edge.length=NULL,where=NULL,position=0,interactive=FALSE,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	use.edge.length<-if(is.null(tree$edge.length)) FALSE else TRUE
	if(use.edge.length==FALSE) tree<-compute.brlen(tree)
	if(interactive==TRUE){
		plotTree(tree,...)
		cat(paste("Click where you would like to bind the tip \"",tip.label,"\"\n",sep=""))
		flush.console()
		obj<-get.treepos(message=FALSE)
		where<-obj$where
		position<-obj$pos
	} else if(is.null(where)) where<-length(tree$tip.label)+1
	if(where<=length(tree$tip.label)&&position==0){
		pp<-1e-12
		if(tree$edge.length[which(tree$edge[,2]==where)]<=1e-12){
			tree$edge.length[which(tree$edge[,2]==where)]<-2e-12
			ff<-TRUE
		} else ff<-FALSE
	} else pp<-position
	if(is.null(edge.length)&&is.ultrametric(tree)){
		H<-nodeHeights(tree)
		if(where==(length(tree$tip.label)+1)) edge.length<-max(H)
		else edge.length<-max(H)-H[tree$edge[,2]==where,2]+position
	}
	tip<-list(edge=matrix(c(2,1),1,2),
		tip.label=tip.label,
		edge.length=edge.length,
		Nnode=1)
		class(tip)<-"phylo"
	obj<-bind.tree(tree,tip,where=where,position=pp)
	if(where<=length(tree$tip.label)&&position==0){
		nn<-obj$edge[which(obj$edge[,2]==which(obj$tip.label==tip$tip.label)),1]
		obj$edge.length[which(obj$edge[,2]==nn)]<-obj$edge.length[which(obj$edge[,2]==nn)]+1e-12
		obj$edge.length[which(obj$edge[,2]==which(obj$tip.label==tip$tip.label))]<-0
		obj$edge.length[which(obj$edge[,2]==which(obj$tip.label==tree$tip.label[where]))]<-0
	}
	root.time<-if(!is.null(obj$root.time)) obj$root.time else NULL
	obj<-untangle(obj,"read.tree")
	if(!is.null(root.time)) obj$root.time<-root.time
	if(interactive) plotTree(obj,...)
	if(!use.edge.length) obj$edge.length<-NULL
	obj
}

# function collapses the subtree descended from node to a star tree
# written by Liam J. Revell 2013, 2015
collapse.to.star<-function(tree,node){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	nel<-if(is.null(tree$edge.length)) TRUE else FALSE
	if(nel) tree$edge.length<-rep(1,nrow(tree$edge))
	tt<-splitTree(tree,split=list(node=node,bp=tree$edge.length[which(tree$edge[,2]==node)]))
	ss<-starTree(species=tt[[2]]$tip.label,branch.lengths=diag(vcv(tt[[2]])))
	ss$root.edge<-0
	tree<-paste.tree(tt[[1]],ss)
	if(nel) tree$edge.length<-NULL 
	tree
}

## function returns the MRCA, or its height above the root, for a set of taxa (in tips)
## written by Liam Revell 2012, 2013, 2015, 2016
findMRCA<-function(tree,tips=NULL,type=c("node","height")){
	type<-type[1]
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.null(tips)){ 
		X<-mrca(tree)
		if(type=="height"){
			H<-nodeHeights(tree)
			X<-apply(X,c(1,2),function(x,y,z) y[which(z==x)[1]],y=H,z=tree$edge)
		}
		return(X)
    } else {
		node<-getMRCA(tree,tips)
		if (type == "node") return(node)
		else if(type=="height") return(nodeheight(tree,node))
	}
}

# function works like extract.clade in ape but will preserve a discrete character mapping
# written by Liam J. Revell 2013
extract.clade.simmap<-function(tree,node){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	x<-getDescendants(tree,node)
	x<-x[x<=length(tree$tip.label)]
	drop.tip.simmap(tree,tree$tip.label[-x])
}

# function gets all subtrees that cannot be further subdivided into two clades of >= clade.size
# written by Liam J. Revell 2013, 2015
getCladesofSize<-function(tree,clade.size=2){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	n<-length(tree$tip.label)
	nn<-1:(tree$Nnode+n)
	ndesc<-function(tree,node){
		x<-getDescendants(tree,node)
		sum(x<=length(tree$tip.label))
	}
	dd<-setNames(sapply(nn,ndesc,tree=tree),nn)
	aa<-n+1 # root
	nodes<-vector()
	while(length(aa)){
		bb<-lapply(aa,function(x,tree) tree$edge[which(tree$edge[,1]==x),2],tree=tree)
		cc<-lapply(bb,function(x) dd[as.character(x)])
		gg<-sapply(cc,function(x,cs) any(x<cs),cs=clade.size)
		nodes<-c(nodes,aa[gg])
		aa<-unlist(bb[!gg])
	}
	trees<-lapply(nodes,extract.clade,phy=tree)
	class(trees)<-"multiPhylo"
	return(trees)
}

# function to get states at internal nodes from simmap style trees
# written by Liam J. Revell 2013, 2014, 2015
getStates<-function(tree,type=c("nodes","tips")){
	type<-type[1]
	if(inherits(tree,"multiPhylo")){
		tree<-unclass(tree)
		obj<-lapply(tree,getStates,type=type)
		nn<-names(obj[[1]])
		y<-sapply(obj,function(x,n) x[n],n=nn)
	} else if(inherits(tree,"phylo")){ 
		if(type=="nodes"){
			y<-setNames(sapply(tree$maps,function(x) names(x)[1]),tree$edge[,1])
			y<-y[as.character(length(tree$tip.label)+1:tree$Nnode)]
		} else if(type=="tips"){
			y<-setNames(sapply(tree$maps,function(x) names(x)[length(x)]),tree$edge[,2])
			y<-setNames(y[as.character(1:length(tree$tip.label))],tree$tip.label)
		}
	} else stop("tree should be an object of class \"phylo\" or \"multiPhylo\".")
	return(y)
}

# function counts transitions from a mapped history
# written by Liam J. Revell 2013, 2015
countSimmap<-function(tree,states=NULL,message=TRUE){
	if(inherits(tree,"multiPhylo")){
		ff<-function(zz){
 			XX<-countSimmap(zz,states,message)
			setNames(c(XX$N,as.vector(t(XX$Tr))),c("N",
			sapply(rownames(XX$Tr),paste,colnames(XX$Tr),sep=",")))
		}
		tree<-unclass(tree)
		XX<-t(sapply(tree,ff))	
		if(!message) return(XX)
		else return(list(Tr=XX,message=
			c("Column N is the total number of character changes on the tree",
			"Other columns give transitions x,y from x->y")))
	} else if(inherits(tree,"phylo")) {
		n<-sum(sapply(tree$maps,length))-nrow(tree$edge)
		if(is.null(states)) states<-colnames(tree$mapped.edge)
		m<-length(states)	
		TT<-matrix(NA,m,m,dimnames=list(states,states))
		gg<-function(map,a,b){
			if(length(map)==1) zz<-0
			else {
				zz<-0; i<-2
				while(i<=length(map)){
					if(names(map)[i]==b&&names(map)[i-1]==a) zz<-zz+1
				i<-i+1
				}
			} 
			return(zz)
		}
		for(i in 1:m) for(j in 1:m)
			if(i==j) TT[i,j]<-0
			else TT[i,j]<-sum(sapply(tree$maps,gg,a=states[i],b=states[j]))
		if(!message) return(list(N=n,Tr=TT))
		else return(list(N=n,Tr=TT,message=c(
			"N is the total number of character changes on the tree",
			"Tr gives the number of transitions from row state->column state")))
	}
}

# function to match nodes between trees
# written by Liam J. Revell 2012, 2013, 2015
matchNodes<-function(tr1,tr2,method=c("descendants","distances"),...){
	if(!inherits(tr1,"phylo")||!inherits(tr1,"phylo")) stop("tr1 & tr2 should both be objects of class \"phylo\".")
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	method<-method[1]
	method<-matchType(method,c("descendants","distances"))
	if(method=="descendants"){
		desc.tr1<-lapply(1:tr1$Nnode+length(tr1$tip),function(x) extract.clade(tr1,x)$tip.label)
		names(desc.tr1)<-1:tr1$Nnode+length(tr1$tip)
		desc.tr2<-lapply(1:tr2$Nnode+length(tr2$tip),function(x) extract.clade(tr2,x)$tip.label)
		names(desc.tr2)<-1:tr2$Nnode+length(tr2$tip)
		Nodes<-matrix(NA,length(desc.tr1),2,dimnames=list(NULL,c("tr1","tr2")))
		for(i in 1:length(desc.tr1)){
			Nodes[i,1]<-as.numeric(names(desc.tr1)[i])
			for(j in 1:length(desc.tr2))
				if(all(desc.tr1[[i]]%in%desc.tr2[[j]])&&all(desc.tr2[[j]]%in%desc.tr1[[i]]))
					Nodes[i,2]<-as.numeric(names(desc.tr2)[j])
		}
	} else if(method=="distances"){
		if(hasArg(tol)) tol<-list(...)$tol
		else tol<-1e-6
		if(hasArg(corr)) corr<-list(...)$corr
		else corr<-FALSE
		if(corr) tr1$edge.length<-tr1$edge.length/max(nodeHeights(tr1))
		if(corr) tr2$edge.length<-tr2$edge.length/max(nodeHeights(tr2))
		D1<-dist.nodes(tr1)[1:length(tr1$tip),1:tr1$Nnode+length(tr1$tip)]
		D2<-dist.nodes(tr2)[1:length(tr2$tip),1:tr2$Nnode+length(tr2$tip)]
		rownames(D1)<-tr1$tip.label
		rownames(D2)<-tr2$tip.label
		common.tips<-intersect(tr1$tip.label,tr2$tip.label)
		D1<-D1[common.tips,]
		D2<-D2[common.tips,]
		Nodes<-matrix(NA,tr1$Nnode,2,dimnames=list(NULL,c("tr1","tr2")))
		for(i in 1:tr1$Nnode){
			if(corr) z<-apply(D2,2,function(X,y) cor(X,y),y=D1[,i])
			else z<-apply(D2,2,function(X,y) 1-sum(abs(X-y)),y=D1[,i])
			Nodes[i,1]<-as.numeric(colnames(D1)[i])
			if(any(z>=(1-tol))){
				a<-as.numeric(names(which(z>=(1-tol))))
				if(length(a)==1) Nodes[i,2]<-a
				else {
					Nodes[i,2]<-a[1]
					if(!quiet) warning("polytomy detected; some node matches may be arbitrary")
				}
			}
		}
	}
	return(Nodes)
}

# function applies the branch lengths of a reference tree to a second tree, including mappings

# written by Liam J. Revell 2012, 2015
applyBranchLengths<-function(tree,edge.length){
	if(inherits(tree,"multiPhylo")){
		trees<-lapply(tree,applyBranchLengths,edge.length=edge.length)
		class(trees)<-"multiPhylo"
		return(trees)
	} else if(inherits(tree,"phylo")) {
		tree$edge.length<-edge.length
		if(!is.null(tree$maps)){
			for(i in 1:nrow(tree$edge)){
				temp<-tree$maps[[i]]/sum(tree$maps[[i]])
				tree$maps[[i]]<-temp*tree$edge.length[i]
			}
		}
		if(!is.null(tree$mapped.edge)){
			a<-vector()
			for(i in 1:nrow(tree$edge)) a<-c(a,names(tree$maps[[i]]))
			a<-unique(a)
			tree$mapped.edge<-matrix(data=0,length(tree$edge.length),length(a),dimnames=list(apply(tree$edge,1,function(x) paste(x,collapse=",")),state=a))
			for(i in 1:length(tree$maps)) for(j in 1:length(tree$maps[[i]])) tree$mapped.edge[i,names(tree$maps[[i]])[j]]<-tree$mapped.edge[i,names(tree$maps[[i]])[j]]+tree$maps[[i]][j]
		}
		return(tree)
	}
}

# function to compute phylogenetic VCV using joint Pagel's lambda
# written by Liam Revell 2011

phyl.vcv<-function(X,C,lambda){
	C<-lambda.transform(lambda,C)
	invC<-solve(C)
	a<-matrix(colSums(invC%*%X)/sum(invC),ncol(X),1)
	A<-matrix(rep(a,nrow(X)),nrow(X),ncol(X),byrow=T)
	V<-t(X-A)%*%invC%*%(X-A)/(nrow(C)-1)
	return(list(C=C,R=V,alpha=a))
}

# lambda transformation of C
# written by Liam Revell 2011
lambda.transform<-function(lambda,C){
	if(lambda==1) return(C)
	else {
		V<-diag(diag(C))
		C<-C-V
		C.lambda<-(V+lambda*C)
		return(C.lambda)
	}
}

# likelihood function for joint estimation of lambda for multiple traits
# written by Liam Revell 2011/2012
likMlambda<-function(lambda,X,C){
	# compute R, conditioned on lambda
	temp<-phyl.vcv(X,C,lambda);
	C<-temp$C; R<-temp$R; a<-temp$alpha
	# prep
	n<-nrow(X); m<-ncol(X); D<-matrix(0,n*m,m)
	for(i in 1:(n*m)) for(j in 1:m) if((j-1)*n<i&&i<=j*n) D[i,j]=1.0
	y<-as.matrix(as.vector(X))
	# compute the log-likelihood
	kronRC<-kronecker(R,C)
	logL<--t(y-D%*%a)%*%solve(kronRC,y-D%*%a)/2-n*m*log(2*pi)/2-determinant(kronRC,logarithm=TRUE)$modulus/2
	return(logL)
}

# function matches data to tree
# written by Liam J. Revell 2011, 2015
matchDatatoTree<-function(tree,x,name){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.matrix(x)) x<-x[,1]
	if(is.null(names(x))){
		if(length(x)==length(tree$tip.label)){
			print(paste(name,"has no names; assuming x is in the same order as tree$tip.label"))
			names(x)<-tree$tip.label
		} else

			stop(paste(name,"has no names and is a different length than tree$tip.label"))

	}
	if(any(is.na(match(names(x),tree$tip.label)))){
		print(paste("some species in",name,"are missing from tree, dropping missing taxa from",name))
		x<-x[intersect(tree$tip.label,names(x))]
	}
	return(x)
}

# function matches tree to data
# written by Liam J. Revell 2011, 2015
matchTreetoData<-function(tree,x,name){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(any(is.na(match(tree$tip.label,names(x))))){
		print(paste("some species in tree are missing from",name,", dropping missing taxa from the tree"))
		tree<-drop.tip(tree,setdiff(tree$tip.label,names(x)))
	}
	if(any(is.na(x))){
		print(paste("some data in",name,"given as 'NA', dropping corresponding species from tree"))
		tree<-drop.tip(tree,names(which(is.na(x))))
	}
	return(tree)
}

# function finds the maximum value of Pagel's lambda
# written by Liam J. Revell 2011
maxLambda<-function(tree){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.ultrametric(tree)){
		H<-nodeHeights(tree)
		return(max(H[,2])/max(H[,1]))
	} else return(1)
}

# function reorders the columns of mapped.edge from a set of simmap trees

# written by Liam J. Revell 2013, 2015
orderMappedEdge<-function(trees,ordering=NULL){
	if(!inherits(trees,"phylo")&&!inherits(trees,"multiPhylo")) 
		stop("trees should be an object of class \"phylo\" or \"multiPhylo\".")
	f1<-function(tree,ordering){
		mapped.edge<-matrix(0,nrow(tree$mapped.edge),length(ordering),
			dimnames=list(rownames(tree$mapped.edge),ordering))
		mapped.edge[,colnames(tree$mapped.edge)]<-tree$mapped.edge
		tree$mapped.edge<-mapped.edge
		return(tree)
	}
	f2<-function(tree) colnames(tree$mapped.edge)
	if(inherits(trees,"phylo")) states<-colnames(trees$mapped.edge)
	else if(inherits(trees,"multiPhylo")) states<-unique(as.vector(sapply(trees,f2)))
	if(length(ordering)>1) if(length(intersect(states,ordering))<length(states)){
		warning("not all states represented in input ordering; setting to default")
		ordering<-NULL
	}
	if(is.null(ordering)) ordering<-"alphabetical"
	if(length(ordering)==1){
		ordering<-matchType(ordering,c("alphabetical","numerical"))
		if(ordering=="alphabetical") ordering<-sort(states)
		else if(ordering=="numerical") ordering<-as.character(sort(as.numeric(states)))
	}
	if(inherits(trees,"phylo")) trees<-f1(trees,ordering)
	else { 
		trees<-lapply(trees,f1,ordering=ordering)
		class(trees)<-"multiPhylo"
	}
	return(trees)
}

# function gets sister node numbers or names
# written by Liam J. Revell 2013, 2015
getSisters<-function(tree,node,mode=c("number","label")){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	mode<-mode[1]
	if(is.character(node)) node<-match(node,c(tree$tip.label,tree$node.label))
	sisters<-tree$edge[which(tree$edge[,1]==tree$edge[which(tree$edge[,2]==node),1]),2]
	sisters<-setdiff(sisters,node)
	if(mode=="number") return(sisters)
	else if(mode=="label"){
		result<-list()
		n<-length(tree$tip.label)
		if(is.null(tree$node.label)&&any(sisters>n)) result$nodes<-sisters[which(sisters>n)] 
		else if(any(sisters>n)) result$nodes<-tree$node.label[sisters[which(sisters>n)]-n]
		if(any(sisters<=n)) result$tips<-tree$tip.label[sisters[which(sisters<=n)]]
		return(result)
	}
}

# gets descendant node numbers
# written by Liam Revell 2012, 2013, 2014
getDescendants<-function(tree,node,curr=NULL){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.null(curr)) curr<-vector()
	daughters<-tree$edge[which(tree$edge[,1]==node),2]
	curr<-c(curr,daughters)
	if(length(curr)==0&&node<=length(tree$tip.label)) curr<-node
	w<-which(daughters>length(tree$tip.label))
	if(length(w)>0) for(i in 1:length(w)) 
		curr<-getDescendants(tree,daughters[w[i]],curr)
	return(curr)
}

# function computes vcv for each state, and stores in array
# written by Liam J. Revell 2011, 2012, 2016
multiC<-function(tree,internal=FALSE){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(!inherits(tree,"simmap")) stop("tree should be an object of class \"simmap\".")
	m<-ncol(tree$mapped.edge)
	# compute separate C for each state
	mC<-list()
	for(i in 1:m){
		mtree<-list(edge=tree$edge,
			Nnode=tree$Nnode,
			tip.label=tree$tip.label,
			edge.length=tree$mapped.edge[,i])
		class(mtree)<-"phylo"
		mC[[i]]<-if(internal) vcvPhylo(mtree,internal=TRUE) else vcv.phylo(mtree)
	}
	names(mC)<-colnames(tree$mapped.edge)
	mC
}

# function pastes subtree onto tip
# written by Liam Revell 2011, 2015
paste.tree<-function(tr1,tr2){
	if(!inherits(tr1,"phylo")||!inherits(tr2,"phylo")) stop("tr1 & tr2 should be objects of class \"phylo\".")
	if(length(tr2$tip)>1){ 
		temp<-tr2$root.edge; tr2$root.edge<-NULL
		tr1$edge.length[match(which(tr1$tip.label=="NA"),tr1$edge[,2])]<-tr1$edge.length[match(which(tr1$tip.label=="NA"),tr1$edge[,2])]+temp
	}
	tr.bound<-bind.tree(tr1,tr2,where=which(tr1$tip.label=="NA"))	
	return(tr.bound)
}

# match type
# written by Liam J. Revell 2012
matchType<-function(type,types){
	for(i in 1:length(types))
		if(all(strsplit(type,split="")[[1]]==strsplit(types[i],split="")[[1]][1:length(strsplit(type,split="")[[1]])]))
			type=types[i]
	return(type)
}

# wraps around MatrixExp
# written by Liam Revell 2011
expm<-function(Y){
	Z<-MatrixExp(Y); dimnames(Z)<-dimnames(Y)
	return(Z)
}
	
# function 'untangles' (or attempts to untangle) a tree with crossing branches
# written by Liam J. Revell 2013, 2015
untangle<-function(tree,method=c("reorder","read.tree")){
	if(inherits(tree,"multiPhylo")){
		tree<-lapply(tree,untangle,method=method)
		class(tree)<-"multiPhylo"
	} else {
		if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
		obj<-attributes(tree)
		method<-method[1]
		if(method=="reorder") tree<-reorder(reorder(tree,"pruningwise"))
		else if(method=="read.tree"){
			if(inherits(tree,"simmap")) tree<-read.simmap(text=write.simmap(tree))
			else tree<-if(Ntip(tree)>1) read.tree(text=write.tree(tree)) else read.newick(text=write.tree(tree))
		}
		ii<-!names(obj)%in%names(attributes(tree))
		attributes(tree)<-c(attributes(tree),obj[ii])
	}
	tree
}

