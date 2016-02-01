## function to plot a tree with dots/circles for a plotted phenotype
## written by Liam J. Revell 2016

dotTree<-function(tree,x,legend=TRUE,...){
	## plot tree
	tree<-reorder(tree,"cladewise")
	plotTree(tree,offset=1.7,ylim=c(-1/25*Ntip(tree),
		Ntip(tree)),...)
	## get last phylo plot parameters
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	x.tip<-obj$xx[1:obj$Ntip]
	y.tip<-obj$yy[1:obj$Ntip]
	## in case any x<0
	min.x<-min(x)
	max.x<-max(x)
	if(any(x<0)) x<-x-min(x)
	## plot points
	x<-x[tree$tip.label]
	draw.circle(x.tip+1.2*strwidth("W"),y.tip,nv=200,
		radius=(0.8*x/max(x)+0.1)/2*diff(par()$usr[1:2])/
		diff(par()$usr[3:4]),col="blue")
	## add legend
	if(legend) dot.legend(x=par()$usr[1]+0.1*max(nodeHeights(tree)),
		y=0.1*(1+par()$usr[3]),min.x,max.x,length=5,...)
}

## dot legend
## written by Liam J. Revell 2016

dot.legend<-function(x,y,min,max,length=5,prompt=FALSE,...){
	if(hasArg(cex)) cex<-list(...)$cex
	else cex<-1
	if(prompt){
		obj<-locator(1)
		x<-obj$x
		y<-obj$y
	}
	text(x,y-0.5,round(min,2),pos=1,cex=cex)
	s<-(0.8*max(min,0)/min(max,max+min)+0.1)/2*diff(par()$usr[1:2])/
		diff(par()$usr[3:4])
	e<-(0.8/2*diff(par()$usr[1:2])+0.1)/diff(par()$usr[3:4])
	rr<-seq(s,e,length.out=length)
	temp<-c(0,cumsum(1.1*rep(2*max(rr),length-1)))
	draw.circle(x+temp,rep(y,length),nv=200,radius=rr,col="blue")
	text(max(x+temp),y-0.5,round(max,2),pos=1,cex=cex)
	y1<-0.1/25*Ntip(tree)
	lines(c(x,max(x+temp)),rep(y-0.5-y1,2))
	lines(c(x,x),y-c(y1+0.5,2*y1+0.5))
	lines(c(max(x+temp),max(x+temp)),y-c(y1+0.5,2*y1+0.5))	
}