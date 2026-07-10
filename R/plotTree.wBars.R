## plotTree.boxplot
## written by Liam J. Revell 2016, 2021, 2022, 2024

plotTree.boxplot<-function(tree,x,args.plotTree=list(),
	args.boxplot=list(),...){
	tree<-untangle(tree,"read.tree")
	cw<-reorder(tree)
	if(inherits(x,"formula")){ 
		obj<-x
	} else {
		if(!is.list(x)){
			obj<-setNames(
				lapply(cw$tip.label,function(x,y) y[which(names(y)==x)],
					y=x),cw$tip.label)
		} else obj<-x
	}
	if(inherits(obj,"formula")) 
		args.boxplot$formula<-obj else args.boxplot$x<-obj
	args.boxplot$horizontal<-TRUE
	args.boxplot$axes<-FALSE
	args.boxplot$names.arg<-""
	args.boxplot$xlim<-c(1,Ntip(cw))
	if(is.null(args.boxplot$space)) args.boxplot$space<-0.7
	if(is.null(args.boxplot$mar)) 
		args.boxplot$mar<-c(5.1,0,2.1,1.1)
	else args.boxplot$mar[2]<-0.1
	args.plotTree$tree<-cw
	if(is.null(args.plotTree$mar)) 
		args.plotTree$mar<-c(5.1,1.1,2.1,0)
	else {
		args.plotTree$mar[4]<-0
	}
	if(args.plotTree$mar[1]!=args.boxplot$mar[1])
		args.plotTree$mar[1]<-args.boxplot$mar[1]
	if(args.plotTree$mar[3]!=args.boxplot$mar[3])
		args.plotTree$mar[3]<-args.boxplot$mar[3]
	if(is.null(args.plotTree$ftype)) args.plotTree$ftype<-"i"
	if(is.null(args.plotTree$lwd)) args.plotTree$lwd<-1
	par(mfrow=c(1,2))
	ii<-which(names(args.boxplot)%in%c("formula","x"))
	args.boxplot<-c(args.boxplot[ii],args.boxplot[-ii])
	args.boxplot$plot<-FALSE
	obj<-do.call(boxplot,args.boxplot)
	N<-ncol(obj$stats)
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-c(0.5,N+0.5)
	args.boxplot$xlim<-ylim
	args.boxplot$plot<-TRUE
	args.plotTree$tips<-setNames(1:Ntip(cw),obj$names)
	args.plotTree$ylim<-ylim
	do.call(plotTree,args.plotTree)
	par(mar=args.boxplot$mar)
	ii<-which(names(args.boxplot)%in%c("formula","x"))
	args.boxplot<-c(args.boxplot[ii],args.boxplot[-ii])
	obj<-do.call(boxplot,args.boxplot)
	axis(1)
	if(!is.null(args.boxplot$xlab)) title(xlab=args.boxplot$xlab)
	else title(xlab="x")
	invisible(obj)
}


## plotTree.barplot
## written by Liam J. Revell 2016, 2017, 2018, 2021, 2024

plotTree.barplot<-function(tree,x,args.plotTree=list(),
	args.barplot=list(),...){
	tree$tip.label<-gsub(" ","_",tree$tip.label)
	tree<-untangle(tree,"read.tree")
	if(hasArg(add)) add<-list(...)$add
	else add<-FALSE
	if(hasArg(args.axis)) args.axis<-list(...)$args.axis
	else args.axis<-list()
	args.axis$side<-1
	cw<-reorder(tree)
	if(is.vector(x)) names(x)<-gsub(" ","_",names(x))
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)){
		x<-x[cw$tip.label,]
		x<-t(x)
		colnames(x)<-gsub(" ","_",colnames(x))
	}
	args.barplot$height<-if(is.matrix(x)) x[,cw$tip.label] else x[cw$tip.label]
	args.barplot$plot<-FALSE
	args.barplot$horiz<-TRUE
	args.barplot$axes<-FALSE
	args.barplot$names.arg<-rep("",Ntip(cw))
	if(is.null(args.barplot$beside)) args.barplot$beside<-FALSE
	if(is.null(args.barplot$space)) 
		args.barplot$space<-if(args.barplot$beside) c(0,1) else 0.7
	if(is.null(args.barplot$mar)) 
		args.barplot$mar<-c(5.1,0,2.1,1.1)
	else args.barplot$mar[2]<-0.1
	obj<-as.matrix(do.call(barplot,args.barplot))
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-c(min(obj)-mean(args.barplot$space),
		max(obj)+mean(args.barplot$space))
	if(dim(obj)[2]==1) obj<-t(obj)
	args.plotTree$tips<-setNames(colMeans(obj),cw$tip.label)
	args.barplot$plot<-TRUE
	args.barplot$ylim<-ylim
	args.plotTree$ylim<-ylim
	args.plotTree$tree<-cw
	if(is.null(args.plotTree$mar)) 
		args.plotTree$mar<-c(5.1,1.1,2.1,0)
	else {
		args.plotTree$mar[4]<-0.1
	}
	if(args.plotTree$mar[1]!=args.barplot$mar[1])
		args.plotTree$mar[1]<-args.barplot$mar[1]
	if(args.plotTree$mar[3]!=args.barplot$mar[3])
		args.plotTree$mar[3]<-args.barplot$mar[3]
	if(is.null(args.plotTree$ftype)) args.plotTree$ftype<-"i"
	if(is.null(args.plotTree$lwd)) args.plotTree$lwd<-1
	if(!add) par(mfrow=c(1,2))
	do.call(plotTree,args.plotTree)
	if(!is.null(args.plotTree$plot)&&args.plotTree$plot==FALSE) par(new=TRUE)
	par(mar=args.barplot$mar)
	obj<-do.call(barplot,args.barplot)
	if(!is.null(args.barplot$xlab)) args.axis$xlab<-args.barplot$xlab
	else args.axis$xlab<-"x"
	do.call(axis,args.axis)	
	invisible(obj)
}

## function to plot bars at the tips of a plotted tree
## written by Liam J. Revell 2014, 2015, 2018, 2019

plotTree.wBars<-function(tree,x,scale=NULL,width=NULL,type="phylogram",
	method="plotTree",tip.labels=FALSE,col="grey",border="black",...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.null(scale)){ 
		scale<-0.3*max(nodeHeights(tree))/diff(range(x))
	}
	if(is.matrix(x)){
		x.neg<-apply(x,1,function(x) sum(x[x<0]))
		x.pos<-apply(x,1,function(x) sum(x[x>0]))
	} else {
		d<-scale*(max(x)-min(0,min(x)))
		H<-nodeHeights(tree)
		if(tip.labels==FALSE){
			lims<-c(-max(H)-d,max(H)+d)
			sw<-0
		} else {
			if(hasArg(fsize)) fsize<-list(...)$fsize
			else fsize<-1
			if(type=="phylogram"){
				pp<-par("pin")[1]
				sw<-fsize*(max(strwidth(tree$tip.label,
					units="inches")))+1.37*fsize*strwidth("W",
					units="inches")
				alp<-optimize(function(a,H,sw,pp,d) 
					(a*1.04*(max(H)+d)+sw-pp)^2,H=H,sw=sw,
					pp=pp,d=d,interval=c(0,1e6))$minimum
				lims<-c(min(H),max(H)+d+sw/alp)
			} else if(type=="fan"){
				pp<-par("pin")[1]
				sw<-fsize*(max(strwidth(tree$tip.label,
					units="inches")))+1.37*fsize*strwidth("W",
					units="inches")
				alp<-optimize(function(a,H,sw,pp,d) 
					(a*2*1.04*(max(H)+d)+2*sw-pp)^2,H=H,sw=sw,pp=pp,
						d=d,interval=c(0,1e6))$minimum
				lims<-c(-max(H)-d-sw/alp,max(H)+d+sw/alp)
			}	
		}
		if(hasArg(lims)) lims<-list(...)$lims
		um<-tree
		if(!is.ultrametric(um)){
			tip.h<-sapply(1:Ntip(tree),nodeheight,tree=tree)
			for(i in 1:Ntip(tree)){	
				ii<-which(um$edge[,2]==i)
				um$edge.length[ii]<-um$edge.length[ii]+(max(tip.h)-tip.h[i])
			}
		}
		if(type=="phylogram"){
			fg<-par()$fg
			if(!is.ultrametric(tree)){
				plotTree(um,ftype=if(tip.labels) "i" else "off",
					xlim=c(0,lims[2]),lwd=1,color="transparent",...)
				for(i in 1:Ntip(tree)) lines(c(max(tip.h),
					tip.h[i]),rep(i,2),lty="dotted")
				add<-TRUE
				par(fg="transparent")
			} else add=FALSE
			if(method=="plotTree") capture.output(plotTree(tree,
				ftype=if(tip.labels) "i" else "off",xlim=c(0,lims[2]),
				add=add,...))
			else if(method=="plotSimmap") capture.output(plotSimmap(tree,
				ftype=if(tip.labels) "i" else "off",xlim=c(0,lims[2]),add=add,...))
			par(fg=fg)
		} else if(type=="fan"){
			fg<-par()$fg
			if(!is.ultrametric(tree)){
				plotTree(um,type="fan",ftype=if(tip.labels) "i" else "off",xlim=lims,ylim=lims,
					lwd=1,color="transparent",...)
				um<-get("last_plot.phylo",envir=.PlotPhyloEnv)
				par(fg="transparent")
				plotTree(tree,type="fan",ftype=if(tip.labels) "i" else "off",xlim=lims,
					ylim=lims,lwd=1,color="transparent",add=TRUE,...)
				tt<-get("last_plot.phylo",envir=.PlotPhyloEnv)
				par(fg="black",lty="solid")
				for(i in 1:Ntip(tree)) lines(c(um$xx[i],tt$xx[i]),c(um$yy[i],tt$yy[i]),lty="dotted")
				par(fg="transparent")
				add<-TRUE
			} else add<-FALSE
			if(method=="plotTree") capture.output(plotTree(tree,type="fan",
				ftype=if(tip.labels) "i" else "off",xlim=lims,ylim=lims,add=add,...))
			else if(method=="plotSimmap") capture.output(plotSimmap(tree,
				type="fan",ftype=if(tip.labels) "i" else "off",xlim=lims,
				ylim=lims,add=add,...))
			par(fg=fg)
		}
		obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		x<-x[tree$tip.label]*scale
		if(is.null(width))
			width<-if(type=="fan") (par()$usr[4]-par()$usr[3])/(max(c(max(x)/max(nodeHeights(tree)),1))*length(tree$tip.label)) 
				else if(type=="phylogram") (par()$usr[4]-par()$usr[3])/(2*length(tree$tip.label))
		w<-width
		if(length(col)<Ntip(tree)) col<-rep(col,ceiling(Ntip(tree)/length(col)))[1:Ntip(tree)]
		if(is.null(names(col))) names(col)<-tree$tip.label
		col<-col[tree$tip.label]
		if(length(border)<Ntip(tree)) border<-rep(border,ceiling(Ntip(tree)/length(border)))[1:Ntip(tree)]
		if(is.null(names(border))) names(border)<-tree$tip.label
		border<-border[tree$tip.label]
		if(type=="phylogram"){
			if(hasArg(direction)) direction<-list(...)$direction
			else direction<-"rightwards"
			sw<-if(tip.labels) fsize*(max(strwidth(tree$tip.label)))+fsize*strwidth("1") else strwidth("l")
			for(i in 1:length(x)){
				dx<-max(obj$xx)
				dy<-obj$yy[i]
				x1<-x2<-dx+sw
				x3<-x4<-x1+x[i]
				y1<-y4<-dy-w/2
				y2<-y3<-dy+w/2
				polygon(c(x1,x2,x3,x4)-min(0,min(x)),
					c(y1,y2,y3,y4),col=col[i],border=border[i])
			}
		} else if(type=="fan"){
			h<-max(nodeHeights(tree))
			sw<-if(tip.labels) fsize*(max(strwidth(tree$tip.label)))+fsize*strwidth("1") else strwidth("l")
			for(i in 1:length(x)){
				theta<-atan(obj$yy[i]/obj$xx[i])
				s<-if(obj$xx[i]>0) 1 else -1
				dx<-s*h*cos(theta)+s*cos(theta)*sw
				dy<-s*h*sin(theta)+s*sin(theta)*sw
				x1<-dx+(w/2)*cos(pi/2-theta)-s*min(0,min(x))*cos(theta)
				y1<-dy-(w/2)*sin(pi/2-theta)-s*min(0,min(x))*sin(theta)
				x2<-dx-(w/2)*cos(pi/2-theta)-s*min(0,min(x))*cos(theta)
				y2<-dy+(w/2)*sin(pi/2-theta)-s*min(0,min(x))*sin(theta)
				x3<-s*x[i]*cos(theta)+x2
				y3<-s*x[i]*sin(theta)+y2
				x4<-s*x[i]*cos(theta)+x1
				y4<-s*x[i]*sin(theta)+y1
				polygon(c(x1,x2,x3,x4),c(y1,y2,y3,y4),col=col[i],
					border=border[i])
			}	
		}
	}
	invisible(obj)
}

