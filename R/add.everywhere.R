# function takes a tree and adds a tip in all possible places
# returns the set of unrooted trees without branch lengths
# written by Liam J. Revell 2011, 2013, 2015

add.everywhere<-function(tree,tip.name){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	tree<-unroot(tree) # unroot tree
	tree$edge.length<-rep(1,nrow(tree$edge)) # set all edge lengths to 1.0
	# create new tip
	new.tip<-list(edge=matrix(c(2L,1L),1,2),tip.label=tip.name,edge.length=1,Nnode=1L)
	class(new.tip)<-"phylo"
	# add the new tip to all edges of the tree
	trees<-list()
	class(trees)<-"multiPhylo"
	for(i in 1:nrow(tree$edge)){
		trees[[i]]<-bind.tree(tree,new.tip,where=tree$edge[i,2],position=0.5)
		trees[[i]]$edge.length<-NULL
	}
	return(trees)
}

