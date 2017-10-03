## plot tree with error bars around divergence times at nodes
## written by Liam J. Revell 2017

plotTree.errorbars<-function(tree,CI,...){
	args<-list(...)
	if(!is.null(args$gridlines)){ 
		gridlines<-args$gridlines
		args$gridlines<-NULL
	} else gridlines<-TRUE
	if(is.null(args$mar)) args$mar<-c(4.1,1.1,1.1,1.1)
	if(is.null(args$ftype)) args$ftype<-"i"
	fsize<-if(!is.null(args$fsize)) args$fsize else 1
	if(is.null(args$direction)) args$direction<-"leftwards"
	if(!is.null(args$bar.width)){
		bar.width<-args$bar.width
		args$bar.width<-NULL
	} else bar.width<-11
	if(!is.null(args$cex)){
		cex<-args$cex
		args$cex<-NULL
	} else cex<-1.2
	if(!is.null(args$bar.col)){
		bar.col<-args$bar.col
		args$bar.col<-NULL
	} else bar.col<-"blue"
	par(mar=args$mar)
	plot.new()		
	th<-max(nodeHeights(tree))
	h<-max(th,max(CI))
	if(is.null(args$xlim)){
		m<-min(min(nodeHeights(tree)),min(CI))
		d<-diff(c(m,h))
		pp<-par("pin")[1]
		sw<-fsize*(max(strwidth(tree$tip.label,units="inches")))+
			1.37*fsize*strwidth("W",units="inches")
		alp<-optimize(function(a,d,sw,pp) (a*1.04*d+sw-pp)^2,
			d=d,sw=sw,pp=pp,
			interval=c(0,1e6))$minimum
		args$xlim<-if(args$direction=="leftwards") c(h,m-sw/alp) else 
			c(m,h+sw/alp)
	}
	if(is.null(args$at)) at<-seq(0,h,by=h/5)
	else {
		at<-args$at
		args$at<-NULL
	}
	args$tree<-tree
	args$add<-TRUE
	do.call(plotTree,args=args)
	if(gridlines) abline(v=at,lty="dashed",
		col=make.transparent("grey",0.5))
	axis(1,at=at,labels=signif(at,3))
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	for(i in 1:tree$Nnode+Ntip(tree))
		lines(x=c(CI[i-Ntip(tree),1],CI[i-Ntip(tree),2]),
			y=rep(obj$yy[i],2),lwd=bar.width,lend=0,
			col=make.transparent(bar.col,0.4))
	points(obj$xx[1:tree$Nnode+Ntip(tree)],
		obj$yy[1:tree$Nnode+Ntip(tree)],pch=19,col=bar.col,
		cex=cex)
}

