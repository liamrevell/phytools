plotTree.lollipop<-function(tree,x,
	args.plotTree=list(),args.lollipop=list(),...){
	if(!inherits(x,c("matrix","data.frame"))) x<-as.matrix(x)
	h<-max(nodeHeights(tree))
	if(hasArg(panel_height)) panel_height<-list(...)$panel_height
	else panel_height<-1.0
	panel_height<-panel_height*h
	args.plotTree$tree<-tree
	args.plotTree$direction<-"upwards"
	if(is.null(args.plotTree$mar)) 
		args.plotTree$mar<-c(0.1,5.1,0.1,0.1)
	if(is.null(args.plotTree$ylim)) 
		args.plotTree$ylim<-c(0,h+ncol(x)*panel_height)
	if(is.null(args.plotTree$ftype)) 
		args.plotTree$ftype<-"off"
	if(is.null(args.plotTree$lwd)) args.plotTree$lwd<-1
	do.call(plotTree,args.plotTree)
	pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(pp$font){
		dx<-abs(diff(pp$x.lim))
		pdin<-par()$din[2]
		sh<-(pp$cex*strwidth(paste(" ",tree$tip.label,sep=""))+
			0.3*pp$cex*strwidth("W"))*(par()$din[1]/par()$din[2])*
			(diff(par()$usr[3:4])/diff(par()$usr[1:2]))
		new_h<-h+max(sh)
		panel_height<-(h-new_h+ncol(x)*panel_height)/ncol(x)
		h<-new_h
	}
	if(hasArg(ylab)) ylab<-list(...)$ylab
	else ylab<-if(!is.null(colnames(x))) 
		colnames(x) else rep("",ncol(x))
	for(i in ncol(x):1){
		d<-max(c(diff(range(x[,i])),max(x[,i])))
		y<-setNames(x[,i]/d*0.8*panel_height,rownames(x))
		lower<-h+(i-1)*panel_height+panel_height*0.05
		upper<-h+(i-1)*panel_height+panel_height*0.95
		polygon(c(0,max(pp$xx)+1,max(pp$xx)+1,0),
			c(lower,lower,upper,upper),
			border=FALSE,col="#F2F2F2")
		hh<-lower-min(c(0,min(y)))+0.05*panel_height
		lines(range(pp$xx),rep(hh,2),col="black",lty="dotted")
		segments(x0=pp$xx[1:Ntip(tree)],y0=rep(hh,Ntip(tree)),
			x1=pp$xx[1:Ntip(tree)],y1=y[tree$tip.label]+hh)
		labs<-pretty(c(min(c(0,min(x[,i]))),x[,i]),n=4)
		labs[!(labs>max(x[,i]))]->labs
		labs[!(labs<min(c(0,min(x[,i]))))]->labs
		axis(2,at=hh+max(y)/max(x[,i])*labs,
			labels=labs,las=1,cex.axis=0.6)
		args.lollipop$bg<-setNames(
			hcl.colors(n=100)[ceiling(99*((y-
				min(y))/diff(range(y))))+1],
			names(y))
		args.lollipop$bg<-args.lollipop$bg[tree$tip.label]
		if(is.null(args.lollipop$pch)) args.lollipop$pch<-21
		if(is.null(args.lollipop$cex)) args.lollipop$cex<-1.2
		args.lollipop$x<-pp$xx[1:Ntip(tree)]
		args.lollipop$y<-y[tree$tip.label]+hh
		do.call(points,args.lollipop)
		mtext(ylab[i],2,line=3,at=mean(hh+max(y)/
			max(x[,i])*labs),cex=0.8)
	}
}
