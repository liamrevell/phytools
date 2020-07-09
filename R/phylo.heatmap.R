## function for phylogenetic heat map
## written by Liam J. Revell 2016, 2017, 2020

phylo.heatmap<-function(tree,X,fsize=1,colors=NULL,standardize=FALSE,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(length(fsize)!=3) fsize<-rep(fsize,3)
	if(hasArg(legend)) legend<-list(...)$legend
	else legend<-TRUE
	if(hasArg(labels)) labels<-list(...)$labels
	else labels<-TRUE
	if(hasArg(split)) split<-list(...)$split
	else split<-c(0.5,0.5)
	split<-split/sum(split)
	if(is.null(colnames(X))) colnames(X)<-paste("var",1:ncol(X),sep="")
	if(standardize){
		sd<-apply(X,2,function(x) sqrt(var(x,na.rm=TRUE)))
		X<-(X-matrix(rep(1,Ntip(tree)),Ntip(tree),1)%*%colMeans(X,na.rm=TRUE))/
			(matrix(rep(1,Ntip(tree)),Ntip(tree),1)%*%sd)
	}
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-c(-0.5,(2-0.5)*split[2]/split[1]+0.5)
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-if(legend) c(if(standardize) -0.15 else -0.1,
		if(labels) 1.1 else 1) else c(0,if(labels) 1.1 else 1)
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-rep(1.1,4)
	if(is.null(colors)) colors<-heat.colors(n=20)[20:1]
	if(hasArg(grid)) add.grid <- list(...)$grid
	else add.grid <- FALSE
	cw<-untangle(as.phylo(tree),"read.tree")
	plot.new()
	par(mar=mar)
	plot.window(xlim=xlim,ylim=ylim)
	h<-phylogram(cw,fsize=fsize[1],...)
	START<-h+1/2*((2-0.5)*split[2]/split[1]+0.5-h)/(ncol(X)-1)+
		0.5*strwidth("W")*fsize[1]
	END<-(2-0.5)*split[2]/split[1]+0.5-1/2*((2-0.5)*split[2]/split[1]+
		0.5-START)/(ncol(X)-1)
	X<-X[cw$tip.label,]
	image(x=seq(START,END,by=(END-START)/(ncol(X)-1)),
		z=t(X[cw$tip.label,]),add=TRUE,
		col=colors,...)
	if(add.grid){
	    dx <- (END - START)/(ncol(X) - 1)
	    x <- seq(START - dx/2, END + dx/2, by = dx) 
	    nTips <- length(tree$tip.label)
	    y <- c(-1/(2*(nTips-1)), seq(0, 1, length = nTips) + 1/(2*(nTips-1)) )
	    segments(x, y[1], x, y[length(y)])
	    segments(x[1], y, x[length(x)], y)
	}
	if(legend){
		if(hasArg(leg.title)) leg.title<-list(...)$leg.title
		else leg.title<-if(standardize) "standardized value" else "value"
		if(hasArg(leg.subtitle)) leg.subtitle<-list(...)$leg.subtitle
		else leg.subtitle<-if(standardize) "SD units" else ""
		add.color.bar(leg=END-START,cols=colors,lims=range(X,na.rm=TRUE),
			title=leg.title,subtitle=leg.subtitle,prompt=FALSE,x=START,
			y=-1/(2*(Ntip(cw)-1))-3*fsize[3]*strheight("W"),
			digits=if(max(abs(X),na.rm=TRUE)<1) round(log10(1/max(abs(X),na.rm=TRUE)))+1 
			else 2,fsize=fsize[3])
	}
	if(labels) text(x=seq(START,END,by=(END-START)/(ncol(X)-1)),
		y=rep(1+1/(2*(Ntip(cw)-1))+0.4*fsize[2]*strwidth("I"),ncol(X)),
		colnames(X),srt=70,adj=c(0,0.5),cex=fsize[2])
	if(any(is.na(X))){
		ii<-which(is.na(X),arr.ind=TRUE)
		x.na<-seq(START,END,by=(END-START)/(ncol(X)-1))[ii[,2]]
		y.na<-seq(0,1,by=1/(nrow(X)-1))[ii[,1]]
		for(i in 1:length(x.na)){
			xx<-x.na[i]+c(1/2,-1/2)*(END-START)/(ncol(X)-1)
			yy<-y.na[i]+c(-1/2,1/2)*1/(nrow(X)-1)
			lines(xx,yy)
		}
	}
}
