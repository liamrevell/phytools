## function to plot a grid of discrete character data next to the tips of a tree
## written by Liam J. Revell 2018, 2020, 2021

plotTree.datamatrix<-function(tree,X,...){
	N<-Ntip(tree)
	ss<-lapply(X,function(x) levels(x))
	k<-sapply(ss,length)
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-40*par()$pin[2]/par()$pin[1]/Ntip(tree)
	if(hasArg(xexp)) xexp<-list(...)$xexp
	else xexp<-1.3
	if(hasArg(yexp)) yexp<-list(...)$yexp
	else yexp<-1.05
	if(hasArg(colors)) colors<-list(...)$colors
	else {
		chk<-.check.pkg("RColorBrewer")
		if(!chk) brewer.pal<-function(...) NULL
		else {
			if(hasArg(palettes)) palettes<-list(...)$palettes
			else {
				palettes<-c("Accent","Dark2","Paired","Pastel1","Pastel2",
					"Set1","Set2","Set3")
			}
			while(length(palettes)<length(k)) palettes<-c(palettes,palettes)
			BREWER.PAL<-function(k,pal){
				ii<-if(k==2) c(1,3) else 1:k
				brewer.pal(max(k,3),pal)[ii]
			}
			colors<-mapply(setNames,mapply(BREWER.PAL,k,
				palettes[1:length(ss)],SIMPLIFY=FALSE),ss,
				SIMPLIFY=FALSE)
        	}
    	}
	if(!is.list(colors)) colors<-list(colors)
	if(!is.null(colnames(X))) names(colors)<-colnames(X) else 
		names(colors)<-rep("",length(colors))
	if(hasArg(sep)) sep<-list(...)$sep
	else sep<-0.5
	if(hasArg(srt)) srt<-list(...)$srt
	else srt<-60
	if(hasArg(space)) space<-list(...)$space
	else space<-0
	if(hasArg(header)) header<-list(...)$header
	else header<-TRUE
	cw<-reorder(tree,"cladewise")
	X<-X[cw$tip.label,,drop=FALSE]
	plotTree(cw,plot=FALSE,fsize=fsize)
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	plotTree(cw,lwd=1,ylim=c(0,obj$y.lim[2]*yexp),
		xlim=c(0,obj$x.lim[2]*xexp),fsize=fsize,
		ftype="off",add=TRUE)
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	h<-max(obj$xx)
	for(i in 1:Ntip(cw)){
		lines(c(obj$xx[i],h),rep(obj$yy[i],2),lty="dotted")
		text(h,obj$yy[i],sub("_"," ",cw$tip.label[i]),cex=fsize,pos=4,
			font=3,offset=0.1)
	}
	s<-max(fsize*strwidth(cw$tip.label))
	start.x<-1.05*h+s
	half<-0.5*(1-space)
	for(i in 1:ncol(X)){
    	if(header) text(start.x,max(obj$yy)+1,colnames(X)[i],pos=4,srt=srt,
			cex=fsize,offset=0)
		for(j in 1:nrow(X)){
			xy<-c(start.x,obj$yy[j])
			y<-c(xy[2]-half,xy[2]+half,xy[2]+half,xy[2]-half)
			asp<-(par()$usr[2]-par()$usr[1])/(par()$usr[4]-par()$usr[3])*
				par()$pin[2]/par()$pin[1]
			x<-c(xy[1]-half*asp,xy[1]-half*asp,xy[1]+half*asp,xy[1]+half*asp)
			polygon(x,y,col=colors[[i]][as.character(X[[i]][j])])
 	   	}
		start.x<-start.x+(1+sep)*asp
	}
	obj<-list(fsize=fsize,
		colors=colors,
		end.x=start.x)
	invisible(obj)
}