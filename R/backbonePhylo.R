## function to manipulate & plot objects of class "backbonePhylo"
## written by Liam J. Revell 2013, 2016

backbone.toTrans<-function(obj){
	tip.label<-obj$tip.label
	clade.label<-sapply(obj$tip.clade,function(x) x$label)
	N<-sapply(obj$tip.clade,function(x) x$N)
	depth<-sapply(obj$tip.clade,function(x) x$depth)
	data.frame(tip.label,clade.label,N,depth)
}

## convert "phylo" to "backbonePhylo"
phylo.toBackbone<-function(x,trans=NULL,...){
	if(!inherits(x,"phylo")) 
		stop("tree should be an object of class \"phylo\".")
	if(hasArg(interactive)) interactive<-list(...)$interactive
	else interactive<-is.null(trans)
	if(hasArg(height)) height<-list(...)$height
	else height<-"mean"
	if(interactive){
		options(locatorBell=FALSE)
		if(is.null(x$node.label)) x$node.label<-as.character(1:x$Nnode+Ntip(x))
		if(inherits(x,"backbonePhylo")) plot(x,...) else plotTree(x,...)
		cat("   Select clade to collapse.\n")
		flush.console()
		check<-textbox(x=c(par()$usr[1],par()$usr[1]+
			0.1*(par()$usr[2]-par()$usr[1])),
			y=par()$usr[4],c("click to stop"),justify="c")
		xy<-unlist(locator(1))
		while(!(xy[1]>par()$usr[1]&&xy[1]<par()$usr[1]+0.1*(par()$usr[2]-
			par()$usr[1])&&xy[2]>check&&xy[2]<par()$usr[4])){
			obj<-get.treepos(message=FALSE,x=xy[1],y=xy[2])
			cat("   Enter name of clade (or press ENTER). > ")
			flush.console()
			clab<-readLines(n=1)
			if(obj$where<=Ntip(x)){
				if(clab=="") clab<-x$tip.label[obj$where]
				tlab<-x$tip.label[obj$where]
				depth<-obj$pos
				N<-if(inherits(x,"backbonePhylo")) 
					x$tip.clade[[obj$where]]$N else 1
				if(!inherits(x,"backbonePhylo")){
					trans<-data.frame(tip.label=as.character(tlab),
						clade.label=as.character(clab),N=N,depth=depth)
					x<-phylo.toBackbone(x,trans)
				} else {
					x$tip.clade[[obj$where]]<-list(label=as.character(clab),
						N=N,depth=depth)
					x$tip.label[obj$where]<-clab
				}
			} else {
				if(clab=="") clab<-x$node.label[obj$where-Ntip(x)]
				tlab<-x$node.label[obj$where-Ntip(x)]
				split<-list(node=obj$where,
					bp=x$edge.length[which(x$edge[,2]==obj$where)]-obj$pos)
				tip.clade<-if(!is.null(x$tip.clade)) x$tip.clade else NULL
				aclass<-class(x)
				tmp<-splitTree(x,split)
				tip<-which(tmp[[1]]$tip.label=="NA")
				if(height=="mean") depth<-obj$pos+
					mean(sapply(1:Ntip(tmp[[2]]),nodeheight,tree=tmp[[2]]))
				else if(height=="max") depth<-obj$pos+
					max(sapply(1:Ntip(tmp[[2]]),nodeheight,tree=tmp[[2]]))
				else if(height=="min") depth<-obj$pos+
					min(sapply(1:Ntip(tmp[[2]]),nodeheight,tree=tmp[[2]]))
				tmp[[1]]$edge.length[which(tmp[[1]]$edge[,2]==tip)]<-
					tmp[[1]]$edge.length[which(tmp[[1]]$edge[,2]==tip)]+depth
				if(inherits(x,"backbonePhylo")){
					dd<-getDescendants(x,obj$where)
					N<-sum(sapply(x$tip.clade[dd[dd<=Ntip(x)]],function(x) x$N))
				} else {
					dd<-getDescendants(x,obj$where)
					N<-sum(dd<=Ntip(x))
				}
				x<-tmp[[1]]
				x$tip.label[tip]<-tlab
				trans<-data.frame(tip.label=as.character(tlab),
					clade.label=as.character(clab),N=N,depth=depth)
				x<-phylo.toBackbone(x,trans)
			}
			if(inherits(x,"backbonePhylo")) plot(x,...) else plotTree(x,...)
			cat("   Select clade to collapse (or STOP).\n")
			flush.console()
			check<-textbox(x=c(par()$usr[1],par()$usr[1]+
				0.1*(par()$usr[2]-par()$usr[1])),
				y=par()$usr[4],c("click to stop"),justify="c")
			xy<-unlist(locator(1))
		}
	} else {
		if(!inherits(x,"backbonePhylo")){
			x$tip.clade<-list()
			for(i in 1:length(x$tip.label)){
				if(x$tip.label[i]%in%trans$tip.label){
					ii<-which(trans$tip.label==x$tip.label[i])
					x$tip.clade[[i]]<-list()
					x$tip.clade[[i]]$label<-as.character(trans$clade.label[ii])
					x$tip.clade[[i]]$N<-trans$N[ii]
					x$tip.clade[[i]]$depth<-trans$depth[ii]
				} else {
					x$tip.clade[[i]]<-list()
					x$tip.clade[[i]]$label<-as.character(x$tip.label[i])
					x$tip.clade[[i]]$N<-1
					x$tip.clade[[i]]$depth<-x$edge.length[which(x$edge[,2]==i)]
				}
			}
			x$tip.label<-sapply(x$tip.clade,function(y) as.character(y$label))
			class(x)<-c("backbonePhylo","phylo")
		} else {
			for(i in 1:nrow(trans)){
				ii<-which(x$tip.label==trans$tip.label[i])
				x$tip.label[ii]<-as.character(trans$clade.label[i])
				x$tip.clade[[length(x$tip.clade)+1]]<-
					list(label=as.character(trans$clade.label[i]),
					N=trans$N[i],depth=trans$depth[i])
				clade.label<-sapply(x$tip.clade,function(y) y$label)
				ii<-sapply(clade.label,function(x,y) if(x%in%y) which(y==x) else -1,
					y=x$tip.label)
				x$tip.clade<-x$tip.clade[ii>0]
				x$tip.clade<-x$tip.clade[ii[ii>0]]
			}
		}
	}
	x
}		

## convert to object of class "phylo"
backbone.toPhylo<-function(x){
	if(!inherits(x,"backbonePhylo")) 
		stop("x not an object of class \"backbonePhylo\"")
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
	cat(paste("\nBackbone phylogenetic tree with",length(x$tip.clade),
		"subtrees and",x$Nnode,"resolved internal nodes.\n"))
	n<-min(length(x$tip.clade),5)
	cat("\nLabels: ")
	cat(paste(sapply(x$tip.clade[1:n],function(y) y$label),collapse=", "))
	cat(", ...\nDiversities: ")
	cat(paste(sapply(x$tip.clade[1:n],function(y) y$N),collapse=", "))
	cat(", ...\n\n")
}

## scale N
scaleN<-function(x,k){
	for(i in 1:length(x$tip.clade)) if(x$tip.clade[[i]]$N>1) 
		x$tip.clade[[i]]$N<-x$tip.clade[[i]]$N*k
	x
}

## plot backbone phylogeny with triangles
plot.backbonePhylo<-function(x,...){
	if(!inherits(x,"backbonePhylo")) 
		stop("x not an object of class \"backbonePhylo\"")
	if(hasArg(vscale)) vscale<-list(...)$vscale
	else vscale<-1
	x<-scaleN(x,vscale)
	tt<-backbone.toPhylo(x)
	n<-sum(sapply(x$tip.clade,function(x) x$N))
	cw<-reorder.backbonePhylo(x,"cladewise")
	y<-vector(length=length(cw$tip.clade)+cw$Nnode)
	ii<-order(cw$edge[,2][cw$edge[,2]<=Ntip(cw)])
	z<-c(0,cumsum(sapply(cw$tip.clade[order(ii)],function(x) x$N)))
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
	plot.new()
	par(mar=rep(0.1,4))
	pp<-par("pin")[1]
	sw<-par("cex")*(max(strwidth(sapply(cw$tip.clade,function(x) x$label),
		units="inches")))+1.37*par("cex")*strwidth("W",units="inches")
	alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=X,sw=sw,pp=pp,
		interval=c(0,1e6))$minimum
	xlim<-c(min(X),max(X)+sw/alp)
	plot.window(xlim=xlim,ylim=c(0,n))
	# plot horizontal edges
	for(i in 1:nrow(X)){
		if(cw$edge[i,2]>length(cw$tip.clade)) lines(X[i,],rep(y[cw$edge[i,2]],2),
			lwd=2,lend=2)
		else lines(c(X[i,1],X[i,2]-cw$tip.clade[[cw$edge[i,2]]]$depth),
			rep(y[cw$edge[i,2]],2),lwd=2,lend=2)
	}
	# plot vertical relationships
	for(i in 1:x$Nnode+length(x$tip.clade)) lines(X[which(cw$edge[,1]==i),1],
		range(y[cw$edge[which(cw$edge[,1]==i),2]]),lwd=2,lend=2)
	for(i in 1:length(x$tip.clade)){
		xx<-c(X[which(cw$edge[,2]==i),2]-cw$tip.clade[[i]]$depth,
			X[which(cw$edge[,2]==i),2],X[which(cw$edge[,2]==i),2])
		yy<-c(y[cw$edge[which(cw$edge[,2]==i),2]],
			y[cw$edge[which(cw$edge[,2]==i),2]]+cw$tip.clade[[i]]$N/2-
			0.5,y[cw$edge[which(cw$edge[,2]==i),2]]-
			cw$tip.clade[[i]]$N/2+0.5)
		if(yy[2]<yy[3]) yy[2]<-yy[3]<-yy[1]
		polygon(x=xx,y=yy,col="grey",lwd=2)
	}
	for(i in 1:length(cw$tip.clade)) 
		text(X[which(cw$edge[,2]==i),2],y[i],cw$tip.clade[[i]]$label,pos=4,
			offset=0.1)
	PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
		show.tip.label=TRUE,show.node.label=FALSE,
		font=1,cex=1,adj=0,srt=0,no.margin=FALSE,label.offset=0.1,
		x.lim=par()$usr[1:2],y.lim=par()$usr[3:4],
		direction="rightwards",tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
		edge=cw$edge,xx=sapply(1:(Ntip(cw)+cw$Nnode),
		function(x,y,z) y[match(x,z)],y=X,z=cw$edge),yy=y)
	assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
}



	

