# function slices a tree at slice and returns all subtrees
# it uses extract.clade(), if trivial==FALSE subtrees with length than 2 taxa are ignored
# written by Liam Revell 2011, 2012, 2015

treeSlice<-function(tree,slice=NULL,trivial=FALSE,prompt=FALSE,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(prompt){
		plotTree(tree,mar=c(3.1,rep(0.1,3)),...)
		axis(1)
		cat("Click at the tree height where cutting is desired...\n")
		flush.console()
		xy<-unlist(locator(1))
		slice<-xy[1]
		cat(paste("Slice height is ",signif(slice,6),". Thank you!\n",sep=""))
		flush.console()
		lines(rep(slice,2),par()$usr[3:4],lty="dashed")
		obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		X<-cbind(obj$xx[obj$edge[,1]],obj$xx[obj$edge[,2]])
		y<-obj$yy[obj$edge[,2]]
		if(trivial){ 
			for(i in 1:nrow(X)) if(X[i,1]<slice&&X[i,2]>slice) 
				points(slice,y[i],pch=19)
		} else {
			for(i in 1:nrow(X)) 
				if(X[i,1]<slice&&X[i,2]>slice&&obj$edge[i,2]>Ntip(tree))
					points(slice,y[i],pch=19)
		}
	}
	tree<-reorder(tree) # reorder cladewise
	H<-nodeHeights(tree)
	edges<-which(H[,2]>slice&H[,1]<slice)
	nodes<-tree$edge[edges,2]
	if(!trivial) nodes<-nodes[nodes>length(tree$tip)]
	trees<-list()
	class(trees)<-"multiPhylo"
	for(i in 1:length(nodes)){
		if(nodes[i]>Ntip(tree)){ 
			trees[[i]]<-extract.clade(tree,node=nodes[i])
			trees[[i]]$root.edge<-H[which(tree$edge[,2]==nodes[i]),2]-slice
		} else {
			z<-list(edge=matrix(c(2,1),1,2),
				edge.length=H[which(tree$edge[,2]==nodes[i]),2]-slice,
				tip.label=tree$tip.label[nodes[i]],Nnode=1L)
			class(z)<-"phylo"
			trees[[i]]<-z
		}
	}
	return(trees)
}