## function to compare two time-trees
## written by Liam J. Revell 2017

compare.chronograms<-function(t1,t2,...){
	if(hasArg(colors)) colors<-list(...)$colors
	else colors<-sapply(c("blue","red"),make.transparent,alpha=0.4)
	if(hasArg(arr.colors)) arr.colors<-list(...)$arr.colors
	else arr.colors<-sapply(c("blue","red"),make.transparent,alpha=0.7)
	h1<-sapply(1:Ntip(t1),nodeheight,tree=t1)
	h2<-sapply(1:Ntip(t2),nodeheight,tree=t2)
	plotTree(if(max(h1)>max(h2)) t1 else t2,plot=FALSE,
		mar=c(4.1,1.1,1.1,1.1),direction="leftwards")
	xlim<-get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim[2:1]
	par(fg="transparent")
	plotTree(t1,color=colors[1],mar=c(4.1,1.1,1.1,1.1),
		xlim=xlim,direction="leftwards",lwd=3)
	T1<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	par(fg="black")
	axis(1)
	par(fg="transparent")
	plotTree(t2,color=colors[2],mar=c(4.1,1.1,1.1,1.1),
		xlim=xlim,add=TRUE,direction="leftwards",ftype="off",lwd=3)
	T2<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	par(fg="black")
	for(i in 1:t1$Nnode+Ntip(t1)){
		arrows(T1$xx[i],T1$yy[i],T2$xx[i],T2$yy[i],lwd=2,
			col=if(T1$xx[i]>T2$xx[i]) arr.colors[2] else arr.colors[1],
			length=0.1)
	}
	h<-mapply(function(x,y) if(x<y) x else y,x=T1$xx[1:Ntip(t1)],
		y=T2$xx[1:Ntip(t2)])
	text(rep(min(h),Ntip(t1)),T1$yy[1:Ntip(t1)],
		labels=t1$tip.label,font=3,pos=4,offset=0.1*max(c(h1,h2)))
	for(i in 1:Ntip(t1)) lines(c(h[i]+
		if(h[i]>min(h)) 0.005*diff(xlim) else 0,min(h)),
		rep(T1$yy[i],2),lty="dotted")	
}

