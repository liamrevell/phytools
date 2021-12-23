## creates an object of class "cophylo"
## written by Liam J. Revell 2015, 2016, 2017, 2019

cophylo<-function(tr1,tr2,assoc=NULL,rotate=TRUE,...){
	if(!inherits(tr1,"phylo")||!inherits(tr2,"phylo")) 
		stop("tr1 & tr2 should be objects of class \"phylo\".")
	## check optional arguments
	if(hasArg(methods)) methods<-list(...)$methods
	else methods<-"pre"
	if("exhaustive"%in%methods) methods<-"exhaustive"
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
		if("exhaustive"%in%methods){
			tt1<-allRotations(tr1)
			tt2<-allRotations(tr2)
			M1<-M2<-matrix(NA,length(tt1),length(tt2))
			for(i in 1:length(tt1)){
				for(j in 1:length(tt2)){
					x<-setNames(sapply(assoc[,2],match,table=tt2[[j]]$tip.label),assoc[,1])
					y<-setNames(sapply(assoc[,1],match,table=tt1[[i]]$tip.label),assoc[,2])
					M1[i,j]<-attr(tipRotate(tt1[[i]],x*Ntip(tr1)/Ntip(tr2),methods="just.compute"),"minRotate")
					M2[i,j]<-attr(tipRotate(tt2[[j]],y*Ntip(tr2)/Ntip(tr1),methods="just.compute"),"minRotate")
				}
			}
			MM<-M1+M2
			ij<-which(MM==min(MM),arr.ind=TRUE)
			obj<-list()
			for(i in 1:nrow(ij)){
				tr1<-tt1[[ij[i,1]]]
				attr(tr1,"minRotate")<-M1[ij[i,1],ij[i,2]]
				tr2<-tt2[[ij[i,2]]]
				attr(tr2,"minRotate")<-M2[ij[i,1],ij[i,2]]
				tt<-list(tr1,tr2)
				class(tt)<-"multiPhylo"
				obj[[i]]<-list(trees=tt,assoc=assoc)
				class(obj[[i]])<-"cophylo"
			}
			if(length(obj)>1) class(obj)<-"multiCophylo"
			else obj<-obj[[1]]
		} else if ("all"%in%methods){
			tt1<-allRotations(tr1)
			tt2<-allRotations(tr2)
			obj<-vector(mode="list",length=length(tt1)*length(tt2))
			ij<-1
			for(i in 1:length(tt1)){
				for(j in 1:length(tt2)){
					x<-setNames(sapply(assoc[,2],match,table=tt2[[j]]$tip.label),assoc[,1])
					y<-setNames(sapply(assoc[,1],match,table=tt1[[i]]$tip.label),assoc[,2])
					obj[[ij]]<-list(trees=c(
						tipRotate(tt1[[i]],x*Ntip(tr1)/Ntip(tr2),methods="just.compute"),
						tipRotate(tt2[[j]],y*Ntip(tr2)/Ntip(tr1),methods="just.compute")),
						assoc=assoc)
					class(obj[[ij]])<-"cophylo"
					ij<-ij+1
				}
			}
			class(obj)<-"multiCophylo"
		} else {
			x<-setNames(sapply(assoc[,2],match,table=tr2$tip.label),assoc[,1])
			tr1<-tipRotate(tr1,x*Ntip(tr1)/Ntip(tr2),right=tr2,assoc=assoc,...)
			best.tr1<-Inf
			x<-setNames(sapply(assoc[,1],match,table=tr1$tip.label),assoc[,2])
			tr2<-tipRotate(tr2,x*Ntip(tr2)/Ntip(tr1),left=tr1,assoc=assoc,...)
			best.tr2<-Inf
			while((best.tr2-attr(tr2,"minRotate"))>0||(best.tr1-attr(tr1,"minRotate"))>0){
				best.tr1<-attr(tr1,"minRotate")
				best.tr2<-attr(tr2,"minRotate")
				x<-setNames(sapply(assoc[,2],match,table=tr2$tip.label),assoc[,1])
				tr1<-tipRotate(tr1,x*Ntip(tr1)/Ntip(tr2),right=tr2,assoc=assoc,...)
				x<-setNames(sapply(assoc[,1],match,table=tr1$tip.label),assoc[,2])
				tr2<-tipRotate(tr2,x*Ntip(tr2)/Ntip(tr1),left=tr1,assoc=assoc,...)
			}
			tt<-list(tr1,tr2)
			class(tt)<-"multiPhylo"
			obj<-list(trees=tt,assoc=assoc)
			class(obj)<-"cophylo"
		}
		cat("Done.\n")
	} else {
		tt<-list(tr1,tr2)
		class(tt)<-"multiPhylo"
		obj<-list(trees=tt,assoc=assoc)
		class(obj)<-"cophylo"
	}
	obj
}

## called internally by plot.cophylo to plot a phylogram
## written by Liam J. Revell

phylogram<-function(tree,part=1,direction="rightwards",fsize=1,ftype="i",lwd=1,...){
	if(hasArg(pts)) pts<-list(...)$pts
	else pts<-TRUE
	if(hasArg(edge.col)) edge.col<-list(...)$edge.col
	else edge.col<-rep("black",nrow(tree$edge))
	if(hasArg(tip.lwd)) tip.lwd<-list(...)$tip.lwd
	else tip.lwd<-1
	if(hasArg(tip.lty)) tip.lty<-list(...)$tip.lty
	else tip.lty<-"dotted"
	if(hasArg(tip.len)) tip.len<-list(...)$tip.len
	else tip.len<-0.1
	if(pts==TRUE&&tip.len==0) tip.len<-0.1
	d<-if(direction=="rightwards") 1 else -1
	## sub "_" for " "
	tree$tip.label<-gsub("_"," ",tree$tip.label)
	## check if edge lenths
	if(is.null(tree$edge.length)) tree<-compute.brlen(tree)
	## rescale tree so it fits in one half of the plot
	## with enough space for labels
	if(ftype=="off") fsize<-0
	n<-Ntip(tree)
	sh<-fsize*strwidth(tree$tip.label)
	H<-nodeHeights(tree)
	th<-sapply(1:n,function(i,x,e) x[which(e==i)],x=H[,2],
		e=tree$edge[,2])+tip.len*max(H)
	tree$edge.length<-tree$edge.length/max(th/(part-sh))
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
	for(i in 1:nrow(X)) lines(d*X[i,],rep(y[cw$edge[i,2]],2),lwd=lwd,lend=2,
		col=edge.col[i])
	## plot vertical relationships
	for(i in 1:tree$Nnode+n){
		ee<-which(cw$edge[,1]==i)
		p<-if(i%in%cw$edge[,2]) which(cw$edge[,2]==i) else NULL
		if(!is.null(p)){
			xx<-c(X[ee,1],X[p,2])
			yy<-sort(c(y[cw$edge[ee,2]],y[cw$edge[p,2]]))
		} else {
			xx<-c(X[ee,1],X[ee[1],1])
			yy<-sort(c(y[cw$edge[ee,2]],mean(y[cw$edge[ee,2]])))
		}
		segments(x0=d*xx[1:(length(xx)-1)],y0=yy[1:(length(yy)-1)],
			x1=d*xx[2:length(xx)],y1=yy[2:length(yy)],lwd=lwd,lend=2,col=edge.col[ee])
	}
	h<-part-0.5-tip.len*(max(X)-min(X))-fsize*strwidth(tree$tip.label)
	## plot links to tips
	for(i in 1:n){ 
		lines(d*c(X[which(cw$edge[,2]==i),2],h[i]+tip.len*(max(X)-min(X))),rep(y[i],2),
			lwd=tip.lwd,lty=tip.lty)
		if(pts) points(d*X[which(cw$edge[,2]==i),2],y[i],pch=16,cex=pts*0.7*sqrt(lwd))
	}	
	## plot tip labels
	font<-which(c("off","reg","b","i","bi")==ftype)-1
	if(font>0){
		for(i in 1:n) TEXTBOX(d*(h[i]+fsize*strwidth(tree$tip.label[i])+
			tip.len*(max(X)-min(X))),y[i],
			tree$tip.label[i], pos=if(d<0) 4 else 2,offset=0,
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
	invisible(d*max(h+fsize*strwidth(tree$tip.label)+tip.len*(max(X)-min(X))))
}

## internally used function

TEXTBOX<-function(x,y,label,pos,offset,cex,font){
	rect(x,y-0.5*strheight(label,cex=cex,font=font),x+if(pos==4) strwidth(label,
		cex=cex,font=font) else -strwidth(label,cex=cex,font=font),
		y+0.5*strheight(label,cex=cex,font=font),border=FALSE,
		col=if(par()$bg%in%c("white","transparent")) "white" else par()$bg)
	text(x=x,y=y,label=label,pos=pos,offset=offset,cex=cex,font=font)
}
	

## plot links between tip taxa according to assoc
## written by Liam J. Revell 2015, 2016, 2019

makelinks<-function(obj,x,link.type="curved",link.lwd=1,link.col="black",
	link.lty="dashed"){
	if(length(link.lwd)==1) link.lwd<-rep(link.lwd,nrow(obj$assoc))
	if(length(link.col)==1) link.col<-rep(link.col,nrow(obj$assoc))
	if(length(link.lty)==1) link.lty<-rep(link.lty,nrow(obj$assoc))
	for(i in 1:nrow(obj$assoc)){
		ii<-which(obj$trees[[1]]$tip.label==obj$assoc[i,1])
		jj<-which(obj$trees[[2]]$tip.label==obj$assoc[i,2])
		for(j in 1:length(ii)) for(k in 1:length(jj)){
			y<-c((ii[j]-1)/(Ntip(obj$trees[[1]])-1),
				(jj[k]-1)/(Ntip(obj$trees[[2]])-1))
			if(link.type=="straight") lines(x,y,lty=link.lty[i],
				lwd=link.lwd[i],col=link.col[i])
			else if(link.type=="curved") drawCurve(x,y,lty=link.lty[i],
				lwd=link.lwd[i],col=link.col[i])
		}
	}
}

## plot method for class "multiCophylo"

plot.multiCophylo<-function(x,...){
	par(ask=TRUE)
	for(i in 1:length(x)) plot.cophylo(x[[i]],...)
}

## plot an object of class "cophylo"
## written by Liam J. Revell 2015, 2016, 2017

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
	if(hasArg(link.type)) link.type<-list(...)$link.type
	else link.type<-"straight"
	if(hasArg(link.lwd)) link.lwd<-list(...)$link.lwd
	else link.lwd<-1
	if(hasArg(link.col)) link.col<-list(...)$link.col
	else link.col<-"black"
	if(hasArg(link.lty))  link.lty<-list(...)$link.lty
	else link.lty<-"dashed"
	if(hasArg(edge.col)) edge.col<-list(...)$edge.col
	else edge.col<-list(
		left=rep("black",nrow(x$trees[[1]]$edge)),
		right=rep("black",nrow(x$trees[[2]]$edge)))
	obj<-list(...)
	if(is.null(obj$part)) obj$part<-0.4
	par(mar=mar)
	plot.window(xlim=xlim,ylim=ylim)
	leftArgs<-rightArgs<-obj
	leftArgs$edge.col<-edge.col$left
	rightArgs$edge.col<-edge.col$right
	if(!is.null(obj$fsize)){
		if(length(obj$fsize)>1){
			leftArgs$fsize<-obj$fsize[1]
			rightArgs$fsize<-obj$fsize[2]
			sb.fsize<- if(length(obj$fsize)>2) obj$fsize[3] else 1
		} else sb.fsize<-1
	} else sb.fsize<-1
	x1<-do.call("phylogram",c(list(tree=x$trees[[1]]),leftArgs))
	left<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	x2<-do.call("phylogram",c(list(tree=x$trees[[2]],direction="leftwards"),
		rightArgs))
	right<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(!is.null(x$assoc)) makelinks(x,c(x1,x2),link.type,link.lwd,link.col,
		link.lty)
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

## print method for "multiCophylo" object

print.multiCophylo<-function(x, ...)
	cat("Object of class \"multiCophylo\" containg",length(x),"objects of class \"cophylo\".\n\n")

## written by Liam J. Revell 2015, 2016

tipRotate<-function(tree,x,...){
	if(hasArg(fn)) fn<-list(...)$fn
	else fn<-function(x) x^2
	if(hasArg(methods)) methods<-list(...)$methods
	else methods<-"pre"
	if("exhaustive"%in%methods) methods<-"exhaustive"
	if(hasArg(print)) print<-list(...)$print
	else print<-FALSE
	if(hasArg(max.exhaustive)) max.exhaustive<-list(...)$max.exhaustive
	else max.exhaustive<-20
	if(hasArg(rotate.multi)) rotate.multi<-list(...)$rotate.multi
	else rotate.multi<-FALSE
	if(rotate.multi) rotate.multi<-!is.binary(tree)
	if(hasArg(anim.cophylo)) anim.cophylo<-list(...)$anim.cophylo
	else anim.cophylo<-FALSE
	if(anim.cophylo){
		if(hasArg(left)) left<-list(...)$left
		else left<-NULL
		if(hasArg(right)) right<-list(...)$right
		else right<-NULL
		if(hasArg(assoc)) assoc<-list(...)$assoc
		else assoc<-NULL
		if(is.null(left)&&is.null(right)) anim.cophylo<-FALSE
		if(hasArg(only.accepted)) only.accepted<-list(...)$only.accepted
		else only.accepted<-TRUE
	}
	tree<-reorder(tree)
	nn<-1:tree$Nnode+length(tree$tip.label)
	if("just.compute"%in%methods){
		foo<-function(phy,x) sum(fn(x-setNames(1:length(phy$tip.label),phy$tip.label)[names(x)]))
		oo<-pp<-foo(tree,x)
	}
	if("exhaustive"%in%methods){
		if(Ntip(tree)>max.exhaustive){
			cat(paste("\nmethods=\"exhaustive\" not permitted for more than",
				max.exhaustive,"tips.\n",
				"If you are sure you want to run an exhaustive search for a tree of this size\n",
				"increasing argument max.exhaustive & re-run.\n",
				"Setting methods to \"pre\".\n\n"))
			methods<-"pre"
		} else {
			cat("Running exhaustive search. May be slow!\n")
			oo<-Inf
			tt<-allRotations(tree)
			foo<-function(phy,x) sum(fn(x-setNames(1:length(phy$tip.label),phy$tip.label)[names(x)]))
			pp<-sapply(tt,foo,x=x)
			ii<-which(pp==min(pp))
			ii<-if(length(ii)>1) sample(ii,1) else ii
			tt<-tt[[ii]]
			pp<-pp[ii]
		}
		if(print) message(paste("objective:",pp))
		tree<-tt
	}
	ANIM.COPHYLO<-function(tree){
		dev.hold()
		if(is.null(left)) plot(cophylo(tree,right,assoc=assoc,rotate=FALSE),...)
		else if(is.null(right)) plot(cophylo(left,tree,assoc=assoc,rotate=FALSE),...)
		nodelabels.cophylo(node=i+Ntip(tree),pie=1,col="red",cex=0.4,
			which=if(is.null(left)) "left" else "right")
		dev.flush()
	}
	if("pre"%in%methods){
		for(i in 1:tree$Nnode){
			if(anim.cophylo) ANIM.COPHYLO(tree)
			tt<-if(rotate.multi) rotate.multi(tree,nn[i]) else untangle(rotate(tree,nn[i]),"read.tree")
			oo<-sum(fn(x-setNames(1:length(tree$tip.label),tree$tip.label)[names(x)]))
			if(inherits(tt,"phylo")) pp<-sum(fn(x-setNames(1:length(tt$tip.label),tt$tip.label)[names(x)]))
			if(anim.cophylo&&!only.accepted) ANIM.COPHYLO(tt)
			else if(inherits(tt,"multiPhylo")){
				foo<-function(phy,x) sum(fn(x-setNames(1:length(phy$tip.label),phy$tip.label)[names(x)]))
				pp<-sapply(tt,foo,x=x)
				if(anim.cophylo&&!only.accepted) nulo<-lapply(tt,ANIM.COPHYLO)
				ii<-which(pp==min(pp))
				ii<-if(length(ii)>1) sample(ii,1) else ii
				tt<-tt[[ii]]
				pp<-pp[ii]
			}
			if(oo>pp) tree<-tt
			if(print) message(paste("objective:",min(oo,pp)))
		}
	}
	if("post"%in%methods){
		for(i in tree$Nnode:1){
			if(anim.cophylo) ANIM.COPHYLO(tree)
			tt<-if(rotate.multi) rotate.multi(tree,nn[i]) else untangle(rotate(tree,nn[i]),"read.tree")
			oo<-sum(fn(x-setNames(1:length(tree$tip.label),tree$tip.label)[names(x)]))
			if(inherits(tt,"phylo")) pp<-sum(fn(x-setNames(1:length(tt$tip.label),tt$tip.label)[names(x)]))
			if(anim.cophylo&&!only.accepted) ANIM.COPHYLO(tt)
			else if(inherits(tt,"multiPhylo")){
				foo<-function(phy,x) sum(fn(x-setNames(1:length(phy$tip.label),phy$tip.label)[names(x)]))
				pp<-sapply(tt,foo,x=x)
				if(anim.cophylo&&!only.accepted) nulo<-lapply(tt,ANIM.COPHYLO)
				ii<-which(pp==min(pp))
				ii<-if(length(ii)>1) sample(ii,1) else ii
				tt<-tt[[ii]]
				pp<-pp[ii]
			}
			if(oo>pp) tree<-tt
			if(print) message(paste("objective:",min(oo,pp)))
		}
	}
	attr(tree,"minRotate")<-min(oo,pp)
	if(anim.cophylo) ANIM.COPHYLO(tree)
	tree
}

## multi2di for "multiPhylo" object

MULTI2DI<-function(x,...){
	obj<-lapply(x,multi2di,...)
	class(obj)<-"multiPhylo"
	obj
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

## function to draw sigmoidal links
## modified from https://stackoverflow.com/questions/32046889/connecting-two-points-with-curved-lines-s-ish-curve-in-r

drawCurve<-function(x,y,scale=0.01,...){
	x1<-x[1]
	x2<-x[2]
	y1<-y[1]
	y2<-y[2]
	curve(plogis(x,scale=scale,location=(x1+x2)/2)*(y2-y1)+y1, 
		x1,x2,add=TRUE,...)
}

## S3 summary method
## written by Liam J. Revell 2016

summary.cophylo<-function(object,...){
	cat("\nCo-phylogenetic (\"cophylo\") object:",deparse(substitute(object)),
		"\n\n")
	cat(paste("Tree 1 (left tree) is an object of class \"phylo\" containing",
		Ntip(object$trees[[1]]),"species.\n\n"))
	cat(paste("Tree 2 (right tree) is an object of class \"phylo\" containing",
		Ntip(object$trees[[2]]),"species.\n\n"))
	cat("Association (assoc) table as follows:\n\n")
	maxl<-max(sapply(strsplit(object$assoc[,1],""),length))
	cat(paste("\tleft:",paste(rep(" ",max(0,maxl-5)),collapse=""),
		"\t----\tright:\n",sep=""))
	nulo<-apply(object$assoc,1,function(x,maxl) cat(paste("\t",x[1],
		paste(rep(" ",maxl-length(strsplit(x[1],split="")[[1]])),
		collapse=""),"\t----\t",x[2],"\n",sep="")),maxl=maxl)
	cat("\n")
}

