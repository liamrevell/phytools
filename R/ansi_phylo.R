## function for ANSI graphic style phylogenetic tree

Asp<-function(){
	w<-par("pin")[1]/diff(par("usr")[1:2])
	h<-par("pin")[2]/diff(par("usr")[3:4])
	w/h
}

ansi_phylo<-function(tree,vertical=c("|","-"),...){
	vertical<-vertical[1]
	if(hasArg(horizontal)) horizontal<-list(...)$horizontal
	else horizontal<-"-"
	if(hasArg(x_spacer)) x_spacer<-list(...)$x_spacer
	else x_spacer<-1
	if(hasArg(y_spacer)) y_spacer<-list(...)$y_spacer
	else y_spacer<-1.4
	args<-list(...)
	args$direction<-"rightwards"
	args$type<-"phylogram"
	args$plot<-FALSE
	args$tree<-tree
	do.call(plotTree,args)
	pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	family<-par()$family
	par(family="mono")
	w<-x_spacer*strwidth(horizontal)
	if(vertical!="-") h<-y_spacer*strheight(vertical)
	else h<-x_spacer*strwidth("-")*Asp()
	ee<-pp$edge
	for(i in 1:nrow(pp$edge)){
		d<-diff(pp$xx[ee[i,]])
		n<-floor(d/w)
		x<-mean(pp$xx[ee[i,]])
		y<-pp$yy[ee[i,2]]
		text(x,y,paste(rep(horizontal,n),collapse=""))
		if(ee[i,2]>Ntip(tree)) {
			dd<-ee[which(ee[,1]==ee[i,2]),2]
			d<-diff(pp$yy[dd])
			n<-floor(d/h)
			if(n>0){
				if(vertical!="-"){
					y<-seq(h/2,by=h,length.out=n)
					y<-(y-mean(y))+mean(pp$yy[dd])
					x<-rep(pp$xx[ee[i,2]],length(y))
					text(x,y,vertical)
				} else {
					y<-mean(pp$yy[dd])
					x<-pp$xx[ee[i,2]]
					text(x,y,paste(rep("-",n),collapse=""),
						srt=90)
				}
			}
		}
	}
	root<-Ntip(tree)+1
	dd<-ee[which(ee[,1]==root),2]
	d<-diff(pp$yy[dd])
	n<-floor(d/h)
	if(n>0){
		if(vertical!="-"){
			y<-seq(h/2,by=h,length.out=n)
			y<-(y-mean(y))+mean(pp$yy[dd])
			x<-rep(pp$xx[root],length(y))
			text(x,y,vertical)
		} else {
			y<-mean(pp$yy[dd])
			x<-mean(pp$xx[root])
			text(x,y,paste(rep("-",n),collapse=""),srt=90)
		}
	}
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-1
	for(i in 1:Ntip(tree)) text(pp$xx[i],pp$yy[i],
		tree$tip.label[i],pos=4,cex=fsize)
	par(family=family)
}
