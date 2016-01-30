## function for phylogenetic heat map
## written by Liam J. Revell 2016

phylo.heatmap<-function(tree,X,fsize=1,colors=NULL,standardize=FALSE,...){
	if(length(fsize)!=3) fsize<-rep(fsize[1],3)
	if(hasArg(legend)) legend<-list(...)$legend
	else legend<-TRUE
	if(hasArg(labels)) labels<-list(...)$labels
	else labels<-TRUE
	if(is.null(colnames(X))) colnames(X)<-paste("var",1:ncol(X),sep="")
	if(standardize){
		sd<-apply(X,2,function(x) sqrt(var(x)))
		X<-(X-matrix(rep(1,Ntip(tree)),Ntip(tree),1)%*%colMeans(X))/
			(matrix(rep(1,Ntip(tree)),Ntip(tree),1)%*%sd)
	}
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-c(-0.5,2)
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-if(legend) c(if(standardize) -0.15 else -0.1,
		if(labels) 1.1 else 1) else c(0,if(labels) 1.1 else 1)
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-rep(1.1,4)
	if(is.null(colors)) colors<-heat.colors(n=20)[20:1]
	cw<-reorder(tree,"cladewise")
	plot.new()
	par(mar=mar)
	plot.window(xlim=xlim,ylim=ylim)
	h<-phylogram(tree,fsize=fsize[1])
	START<-h+1/2*(2-h)/(ncol(X)-1)+0.5*strwidth("W")*fsize[1]
	END<-2-1/2*(2-START)/(ncol(X)-1)
	image(x=seq(START,END,by=(END-START)/(ncol(X)-1)),
		z=t(X[cw$tip.label,]),add=TRUE,
		col=colors,...)
	add.color.bar(leg=END-START,cols=colors,lims=range(X),
		title=if(standardize) "standardized value" else "value",
		subtitle=if(standardize) "SD units" else "",prompt=FALSE,x=START,
		y=-1/(2*(Ntip(tree)-1))-2*fsize[3]*strheight("W"),
		digits=if(max(abs(X))<1) round(log10(1/max(abs(X))))+1 else 2,
		fsize=fsize[3])
	text(x=seq(START,END,by=(END-START)/(ncol(X)-1)),
		y=rep(1+1/(2*(Ntip(tree)-1))+0.4*fsize[2]*strwidth("I"),ncol(X)),
		colnames(X),srt=70,adj=c(0,0.5),cex=fsize[2])
}
