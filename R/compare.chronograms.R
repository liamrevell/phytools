## function to compare two time-trees
## written by Liam J. Revell 2017

compare.chronograms<-function(t1,t2,...){
	if(hasArg(colors)) colors<-list(...)$colors
	else colors<-sapply(c("blue","red"),make.transparent,alpha=0.4)
	if(hasArg(arr.colors)) arr.colors<-list(...)$arr.colors
	else arr.colors<-sapply(c("blue","red"),make.transparent,alpha=0.7)
	h1<-max(nodeHeights(t1))
	h2<-max(nodeHeights(t2))
	plotTree(if(h1>h2) t1 else t2,plot=FALSE,mar=c(4.1,1.1,1.1,1.1))
	xlim<-get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim[2:1]
	plotTree(t1,color=colors[1],mar=c(4.1,1.1,1.1,1.1),
		xlim=xlim,direction="leftwards",lwd=3)
	T1<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	axis(1)
	plotTree(t2,color=colors[2],mar=c(4.1,1.1,1.1,1.1),
		xlim=xlim,add=TRUE,direction="leftwards",ftype="off",lwd=3)
	T2<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	for(i in 1:t1$Nnode+Ntip(t1)){
		arrows(T1$xx[i],T1$yy[i],T2$xx[i],T2$yy[i],lwd=2,
			col=if(T1$xx[i]>T2$xx[i]) arr.colors[2] else arr.colors[1],
			length=0.1)
	}
}