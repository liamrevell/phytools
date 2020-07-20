## function to plot probability or trait value by branch
## written by Liam J. Revell 2013, 2014, 2016, 2020

plotBranchbyTrait<-function(tree,x,mode=c("edges","tips","nodes"),palette="rainbow",
	legend=TRUE,xlims=NULL,...){
	mode<-mode[1]
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(mode=="tips"){
		x<-c(x[tree$tip.label],fastAnc(tree,x))
		names(x)[1:length(tree$tip.label)]<-1:length(tree$tip.label)
		XX<-matrix(x[tree$edge],nrow(tree$edge),2)
		x<-rowMeans(XX)
	} else if(mode=="nodes"){
		XX<-matrix(x[tree$edge],nrow(tree$edge),2)
		x<-rowMeans(XX)
	}
	# begin optional arguments
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-6
	if(hasArg(prompt)) prompt<-list(...)$prompt
	else prompt<-FALSE
	if(hasArg(type)) type<-list(...)$type
	else type<-"phylogram"
	if(hasArg(show.tip.label)) show.tip.label<-list(...)$show.tip.label
	else show.tip.label<-TRUE
	if(hasArg(show.node.label)) show.node.label<-list(...)$show.node.label
	else show.node.label<-FALSE
	if(hasArg(edge.width)) edge.width<-list(...)$edge.width
	else edge.width<-4
	if(hasArg(edge.lty)) edge.lty<-list(...)$edge.lty
	else edge.lty<-1
	if(hasArg(font)) font<-list(...)$font
	else font<-3
	if(hasArg(cex)) cex<-list(...)$cex
	else cex<-par("cex")
	if(length(cex)==1) cex<-rep(cex,2)
	if(hasArg(adj)) adj<-list(...)$adj
	else adj<-NULL
	if(hasArg(srt)) srt<-list(...)$srt
	else srt<-0
	if(hasArg(no.margin)) no.margin<-list(...)$no.margin
	else no.margin<-TRUE
	if(hasArg(root.edge)) root.edge<-list(...)$root.edge
	else root.edge<-FALSE
	if(hasArg(label.offset)) label.offset<-list(...)$label.offset
	else label.offset<-0.01*max(nodeHeights(tree))
	if(hasArg(underscore)) underscore<-list(...)$underscore
	else underscore<-FALSE
	if(hasArg(x.lim)) x.lim<-list(...)$x.lim
	else x.lim<-NULL
	if(hasArg(y.lim)) y.lim<-list(...)$y.lim
	else y.lim<-if(legend&&!prompt&&type%in%c("phylogram","cladogram")) 
		c(1-0.06*length(tree$tip.label),length(tree$tip.label)) else NULL
	if(hasArg(direction)) direction<-list(...)$direction
	else direction<-"rightwards"
	if(hasArg(lab4ut)) lab4ut<-list(...)$lab4ut
	else lab4ut<-NULL
	if(hasArg(tip.color)) tip.color<-list(...)$tip.color
	else tip.color<-"black"
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	if(hasArg(rotate.tree)) rotate.tree<-list(...)$rotate.tree
	else rotate.tree<-0
	if(hasArg(open.angle)) open.angle<-list(...)$open.angle
	else open.angle<-0
	# end optional arguments
	if(is.function(palette)) cols<-palette(n=1000)
	else {
		if(palette=="heat.colors") cols<-heat.colors(n=1000)
		if(palette=="gray") cols<-gray(1000:1/1000)
		if(palette=="rainbow")	cols<-rainbow(1000,start=0.7,end=0) # blue->red
	}
	if(is.null(xlims)) xlims<-range(x)+c(-tol,tol)
	breaks<-0:1000/1000*(xlims[2]-xlims[1])+xlims[1]
	whichColor<-function(p,cols,breaks){
		i<-1
		while(p>=breaks[i]&&p>breaks[i+1]) i<-i+1
		cols[i]
	}
	colors<-sapply(x,whichColor,cols=cols,breaks=breaks)
	par(lend=2)
	# now plot
	xx<-plot.phylo(tree,type=type,show.tip.label=show.tip.label,show.node.label=show.node.label,
		edge.color=colors,edge.width=edge.width,edge.lty=edge.lty,font=font,cex=cex[1],adj=adj,
		srt=srt,no.margin=no.margin,root.edge=root.edge,label.offset=label.offset,
		underscore=underscore,x.lim=x.lim,y.lim=y.lim,direction=direction,lab4ut=lab4ut,
		tip.color=tip.color,plot=plot,rotate.tree=rotate.tree,open.angle=open.angle,lend=2,
		new=FALSE)
	if(legend==TRUE&&is.logical(legend)) legend<-round(0.3*max(nodeHeights(tree)),2)
	if(legend){
		if(hasArg(title)) title<-list(...)$title
		else title<-"trait value"
		if(hasArg(digits)) digits<-list(...)$digits
		else digits<-1
		if(prompt) add.color.bar(legend,cols,title,xlims,digits,prompt=TRUE,fsize=cex[2])
		else add.color.bar(legend,cols,title,xlims,digits,prompt=FALSE,
			x=par()$usr[1]+0.05*(par()$usr[2]-par()$usr[1]),
			y=par()$usr[3]+0.05*(par()$usr[4]-par()$usr[3]),fsize=cex[2])
	}
	invisible(xx)
}

# function to add color bar
# written by Liam J. Revell 2013, 2015, 2016

add.color.bar<-function(leg,cols,title=NULL,lims=c(0,1),digits=1,prompt=TRUE,lwd=4,outline=TRUE,...){
	if(prompt){
		cat("Click where you want to draw the bar\n")
		flush.console()
		x<-unlist(locator(1))
		y<-x[2]
		x<-x[1]
	} else {
		if(hasArg(x)) x<-list(...)$x
		else x<-0
		if(hasArg(y)) y<-list(...)$y
		else y<-0
	}
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-1.0
	if(hasArg(subtitle)) subtitle<-list(...)$subtitle
	else subtitle<-NULL
	if(hasArg(direction)) direction<-list(...)$direction
	else direction<-"rightwards"
	if(direction%in%c("rightwards","leftwards")){
		X<-x+cbind(0:(length(cols)-1)/length(cols),1:length(cols)/length(cols))*(leg)
		if(direction=="leftwards"){ 
			X<-X[nrow(X):1,]
			if(!is.null(lims)) lims<-lims[2:1]
		}
		Y<-cbind(rep(y,length(cols)),rep(y,length(cols)))
	} else if(direction%in%c("upwards","downwards")){
		Y<-y+cbind(0:(length(cols)-1)/length(cols),1:length(cols)/length(cols))*(leg)
		if(direction=="downwards"){ 
			X<-X[nrow(X):1,]
			if(!is.null(lims)) lims<-lims[2:1]
		}
		X<-cbind(rep(x,length(cols)),rep(x,length(cols)))
	}
	if(outline) lines(c(X[1,1],X[nrow(X),2]),c(Y[1,1],Y[nrow(Y),2]),lwd=lwd+2,lend=2) 
	for(i in 1:length(cols)) lines(X[i,],Y[i,],col=cols[i],lwd=lwd,lend=2)
	if(direction%in%c("rightwards","leftwards")){
		if(!is.null(lims)) text(x=x,y=y,
			round(lims[1],digits),pos=3,cex=fsize)
		if(!is.null(lims)) text(x=x+leg,y=y,
			round(lims[2],digits),pos=3,cex=fsize)
		if(is.null(title)) title<-"P(state=1)"
		text(x=(2*x+leg)/2,y=y,title,pos=3,cex=fsize)
		if(is.null(subtitle)) 
			text(x=(2*x+leg)/2,y=y,paste("length=",round(leg,3),sep=""),pos=1,cex=fsize)
		else text(x=(2*x+leg)/2,y=y,subtitle,pos=1,cex=fsize)
	} else if(direction%in%c("upwards","downwards")){
		if(!is.null(lims)) text(x=x,y=y-0.02*diff(par()$usr[3:4]),round(lims[1],digits),
			pos=1,cex=fsize)
		if(!is.null(lims)) text(x=x,y=y+leg+0.02*diff(par()$usr[3:4]),
			round(lims[2],digits),
			pos=3,cex=fsize)
		if(is.null(title)) title<-"P(state=1)"
		text(x=x-0.04*diff(par()$usr[1:2]),y=(2*y+leg)/2,title,
			pos=3,cex=fsize,srt=90)
		if(is.null(subtitle)) 
			text(x=x+0.04*diff(par()$usr[1:2]),y=(2*y+leg)/2,
				paste("length=",round(leg,3),sep=""),pos=1,
				srt=90,cex=fsize)
		else text(x=x+0.04*diff(par()$usr[1:2]),y=(2*y+leg)/2,
			subtitle,pos=1,cex=fsize,srt=90)
	}
}

## function to create "edge.widthMap" object
edge.widthMap<-function(tree,x,...){
	if(!inherits(tree,"phylo")) 
		stop("tree should be an object of class \"phylo\".")
	tree<-as.phylo(tree)
	a<-fastAnc(tree,x)
	node.values<-c(x[tree$tip.label],a)
	edge.values<-apply(tree$edge,1,function(e,nv)
		mean(nv[e]),nv=node.values)
	edge.widths<-edge.values
	object<-list(tree=tree,edge.widths=edge.widths,
		node.values=node.values)
	class(object)<-"edge.widthMap"
	object
}

## print method
print.edge.widthMap<-function(x,...){
	cat("Object of class \"edge.widthMap\" containing:\n")
	cat(paste("(1) Phylogenetic tree with",Ntip(x$tree),
		"tips and",x$tree$Nnode,"internal nodes.\n"))
	cat("(2) Vector of node values for a mapped quantitative\n")
	cat("    trait.\n\n")
}

## plot method
plot.edge.widthMap<-function(x,max.width=0.9,legend="trait value",...){
	if(hasArg(min.width)) min.width<-list(...)$min.width
	else min.width<-0
	if(hasArg(vertical.as.polygon)) 
		vertical.as.polygon<-list(...)$vertical.as.polygon
	else vertical.as.polygon<-TRUE
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	h<-max(nodeHeights(x$tree))
	node.values<-x$node.values-min(x$node.values)
	node.values<-node.values*((max.width-min.width)/
		max(node.values))+min.width
	args.list<-list(...)
	args.list$tree<-x$tree
	args.list$type<-"phylogram"
	if(!is.null(args.list$direction)){
		if(!args.list$direction%in%c("leftwards","rightwards"))
			args.list$direction<-"rightwards"
	} else args.list$direction<-"rightwards"
	if(is.null(args.list$ylim)) 
		args.list$ylim<-c(1,Ntip(x$tree)+Ntip(x$tree)/25)
	if(is.null(args.list$ftype)) args.list$ftype<-"i"
	if(is.null(args.list$fsize)) 
		args.list$fsize<-36*par()$pin[2]/par()$pin[1]/
			Ntip(x$tree)
	if(is.null(args.list$color)){ 
		args.list$color<-"transparent"
		color<-"gray62"
	} else {
		color<-args.list$color
		args.list$color<-"transparent"
	}
	if(is.null(args.list$border)){
		border<-color
	} else {
		border<-args.list$border
		args.list$border<-NULL
	}
	do.call(plotTree,args.list)
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	asp<-par()$pin[2]/par()$pin[1]/
		(diff(par()$usr[4:3])/diff(par()$usr[2:1]))
	for(i in 1:nrow(x$tree$edge)){
		if(vertical.as.polygon){
			xx<-obj$xx[c(x$tree$edge[i,1],
				x$tree$edge[i,1:2],
				x$tree$edge[i,2:1],
				x$tree$edge[i,1])]+c(0,
				asp*node.values[x$tree$edge[i,1]]/2,
				0,0,asp*node.values[x$tree$edge[i,1]]/2,0)
			yy<-rep(obj$yy[x$tree$edge[i,2]],6)+
				c(node.values[x$tree$edge[i,1]],
				node.values[x$tree$edge[i,1:2]],
				-node.values[x$tree$edge[i,2:1]],
				-node.values[x$tree$edge[i,1]])/2
		}
		else {
			xx<-obj$xx[c(x$tree$edge[i,1:2],
				x$tree$edge[i,2:1])]
			yy<-rep(obj$yy[x$tree$edge[i,2]],4)+
				c(node.values[x$tree$edge[i,1:2]],
				-node.values[x$tree$edge[i,2:1]])/2
		}
		polygon(x=crop.to.h(xx,h),y=yy,
			border=border,col=color,lwd=lwd)	
	}
	for(i in 1:x$tree$Nnode+Ntip(x$tree)){
		nn<-x$tree$edge[which(x$tree$edge[,1]==i),2]
		yy<-range(obj$yy[nn])
		if(vertical.as.polygon){
			xx<-rep(obj$xx[i],4)+
				asp*c(-node.values[i]/2,node.values[i]/2,
				node.values[i]/2,-node.values[i]/2)
			polygon(x=crop.to.h(xx,h),
				y=c(rep(yy[1],2),rep(yy[2],2))+
				c(-rep(node.values[i],2),
				rep(node.values[i],2))/2,
				border=border,col=color,lwd=lwd)
		} else {
			lines(rep(obj$xx[i],2),yy+c(-node.values[i],
				node.values[i])/2,lend=2,col=border,
				lwd=lwd)
		}
	}
	if(border!=color&&vertical.as.polygon){
		for(i in 1:nrow(x$tree$edge)){
			xx<-obj$xx[c(x$tree$edge[i,1],
				x$tree$edge[i,1:2],
				x$tree$edge[i,2:1],
				x$tree$edge[i,1])]+c(0,
				asp*node.values[x$tree$edge[i,1]]/2,
				0,0,asp*node.values[x$tree$edge[i,1]]/2,0)
			yy<-rep(obj$yy[x$tree$edge[i,2]],6)+
				c(node.values[x$tree$edge[i,1]],
				node.values[x$tree$edge[i,1:2]],
				-node.values[x$tree$edge[i,2:1]],
				-node.values[x$tree$edge[i,1]])/2
			polygon(x=crop.to.h(xx,h),y=yy,
				border=FALSE,col=color)
		}
		for(i in 1:x$tree$Nnode+Ntip(x$tree)){
			nn<-x$tree$edge[which(x$tree$edge[,1]==i),2]
			yy<-range(obj$yy[nn])
			xx<-rep(obj$xx[i],4)+
				asp*c(-node.values[i]/2,node.values[i]/2,
				node.values[i]/2,-node.values[i]/2)
			polygon(x=crop.to.h(xx,h),
				y=c(rep(yy[1],2),rep(yy[2],2))+
				c(-rep(node.values[i],2),
				rep(node.values[i],2))/2,
				border=FALSE,col=color)
		}
	}
	leg.length<-0.4*h
	x.adj<-if(obj$direction=="rightwards") 0 else obj$x.lim[2]-leg.length
	polygon(x=c(0,0,leg.length,leg.length)+x.adj,
		y=Ntip(x$tree)+Ntip(x$tree)/25+
		c(-min.width/2,min.width/2,max(node.values)/2,
		-max(node.values)/2),
		border=border,col=color,lwd=lwd)
	if(border!=color)
		polygon(x=c(0,0,leg.length,leg.length)+x.adj,
			y=Ntip(x$tree)+Ntip(x$tree)/25+
			c(-min.width/2,min.width/2,max(node.values)/2,
			-max(node.values)/2),
			border=FALSE,col=color)
	text(0+x.adj,Ntip(x$tree)+Ntip(x$tree)/25-0.2*max.width,
		round(min(x$node.values),2),pos=1,
		cex=0.8)
	text(leg.length+x.adj,Ntip(x$tree)+Ntip(x$tree)/25-0.2*max.width,
		round(max(x$node.values),2),pos=1,cex=0.8)
	text(leg.length/2+x.adj,Ntip(x$tree)+Ntip(x$tree)/25-0.2*max.width,
		legend,pos=1,cex=0.8)
	
}

crop.to.h<-function(x,h) sapply(x,function(x,h) if(x<=h) x else h,h=h)