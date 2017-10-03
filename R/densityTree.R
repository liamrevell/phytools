## function to make a color (e.g., "blue") transparent with alpha level alpha
make.transparent<-function(color,alpha){
	RGB<-col2rgb(color)[,1]/255
	rgb(RGB[1],RGB[2],RGB[3],alpha)
}

## function to be used internally
rescaleTree<-function(tree,scale){
	tree$edge.length<-tree$edge.length/max(nodeHeights(tree))*scale
	tree
}
## function to plot a posterior density of trees (e.g., densiTree in phangorn)
## written by Liam J. Revell 2016, 2017
densityTree<-function(trees,colors="blue",alpha=NULL,method="plotTree",
	fix.depth=FALSE,use.edge.length=TRUE,compute.consensus=TRUE,
	use.gradient=FALSE,show.axis=TRUE,...){
	N<-length(trees)
	if(any(sapply(trees,function(x) is.null(x$edge.length)))) 
		use.edge.length<-FALSE
	if(!use.edge.length){ 
		trees<-lapply(trees,compute.brlen)
		class(trees)<-"multiPhylo"
	}
	h<-sapply(trees,function(x) max(nodeHeights(x)))
	if(fix.depth){
		if(method=="plotTree"){
			trees<-lapply(trees,rescaleTree,mean(h))
			class(trees)<-"multiPhylo"
		} else if(method=="plotSimmap"){ 
			trees<-rescaleSimmap(trees,depth=mean(h))
		}
		h<-sapply(trees,function(x) max(nodeHeights(x)))
	}
	tips<-setNames(1:Ntip(trees[[1]]), 
		if(compute.consensus) untangle(consensus(trees),
			"read.tree")$tip.label 
		else trees[[1]]$tip.label)
	if(is.null(alpha)) alpha<-max(c(1/N,0.01))
	args<-list(...)
	args$direction<-"leftwards"
	args$tips<-tips
	args$add<-FALSE
	if(is.null(args$nodes)) args$nodes<-"inner"
	if(is.null(args$mar)) args$mar<-if(show.axis) c(4.1,1.1,1.1,1.1) else rep(1.1,4)
	if(is.null(args$ftype)) args$ftype<-"i"
	if(!use.gradient){
		plotTree(trees[[which(h==max(h))[1]]],direction="leftwards",mar=args$mar,
			plot=FALSE)
		args$xlim<-get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim[2:1]
		if(method=="plotTree"){
			args$color<-make.transparent(colors[1],alpha)
			for(i in 1:length(trees)){
				args$tree<-trees[[i]]
				do.call(plotTree,args)
				if(i==1){ 
					if(show.axis) axis(1)
					args$ftype<-"off"
					args$add<-TRUE
				}
			}
		} else if(method=="plotSimmap"){
			states<-sort(unique(as.vector(mapped.states(trees))))
			if(length(colors)!=length(states)){
				colors<-setNames(c("grey",palette()[2:length(states)]),
					states)
			}
			colors<-sapply(colors,make.transparent,alpha=alpha)
			args$colors<-colors
			for(i in 1:length(trees)){
				args$tree<-trees[[i]]
				do.call(plotSimmap,args)
				if(i==1){ 
					if(show.axis) axis(1)
					args$ftype<-"off"
					args$add<-TRUE
				}
			}
		}
	} else if(use.gradient){
		rf<-multiRF(trees,quiet=TRUE)
		mds<-cmdscale(rf,k=1)[,1]
		trees<-trees[order(mds)]
		h<-h[order(mds)]
		args$ylim<-c(0,Ntip(trees[[1]])+1)
		plotTree(trees[[which(h==max(h))[1]]],direction="leftwards",mar=args$mar,
			ylim=args$ylim,plot=FALSE)
		args$xlim<-get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim[2:1]
		colors<-sapply(rainbow(n=length(trees)),make.transparent,alpha=alpha)
		ftype<-args$ftype
		for(i in 1:length(trees)){
			y.shift<-(i-median(1:length(trees)))/length(trees)/2
			args$tree<-trees[[i]]
			args$tips<-tips+y.shift
			args$color<-colors[i]
			args$ftype<-if(i==floor(median(1:length(trees)))) ftype else "off"
			do.call(plotTree,args)
			if(i==1){ 
				if(show.axis) axis(1)
				args$ftype<-"off"
				args$add<-TRUE
			}
		}
	}
}
