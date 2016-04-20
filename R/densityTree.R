## function to make a color (e.g., "blue") transparent with alpha level alpha
make.transparent<-function(color,alpha){
	RGB<-col2rgb(color)[,1]/255
	rgb(RGB[1],RGB[2],RGB[3],alpha)
}

## function to plot a posterior density of trees (e.g., densiTree in phangorn)
## written by Liam J. Revell 2016
densityTree<-function(trees,colors="blue",alpha=NULL,method="plotTree",
	fix.depth=FALSE,use.edge.length=TRUE,compute.consensus=FALSE,...){
	N<-length(trees)
	if(any(sapply(trees,function(x) is.null(x$edge.length)))) 
		use.edge.length<-FALSE
	if(!use.edge.length) trees<-lapply(trees,compute.brlen)
	if(!fix.depth){
		h<-sapply(trees,function(x) max(nodeHeights(x)))
		ii<-order(h,decreasing=TRUE)
		trees<-trees[ii]
		h<-h[ii]
	}
	tips<-setNames(1:Ntip(trees[[1]]), 
		if(compute.consensus) untangle(consensus(trees),
			"read.tree")$tip.label 
		else trees[[1]]$tip.label)
	if(is.null(alpha)) alpha<-1/N
	colors<-setNames(sapply(colors,make.transparent,alpha),names(colors))
	if(method=="plotSimmap") foo<-plotSimmap else foo<-plotTree
	foo(trees[[1]],color=colors,tips=tips,...)
	xlim<-get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim
	xlim[1]<-xlim[1]+0.03703704*diff(xlim)
	xlim[2]<-xlim[2]-0.03703704*diff(xlim)
	par(fg="transparent")
	for(i in 2:length(trees))
		foo(trees[[i]],
			tips=tips,
			color=colors,add=TRUE,
			xlim=if(fix.depth) NULL else xlim-(h[1]-h[i]),...)
	par(fg="black")
}