# function adds a set of tips to random positions in the tree
# written by Liam J. Revell 2013, 2015

add.random<-function(tree,n=NULL,tips=NULL,edge.length=NULL,order=c("random","input")){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	# internal function
	randomPosn<-function(tree){
		# sum edges cumulatively
		cum.edge<-cumsum(tree$edge.length)
		index<-tree$edge[,2]
		# pick random position
		pos<-runif(1)*sum(tree$edge.length)
		edge<-1; while(pos>cum.edge[edge]) edge<-edge+1
		return(list(node=index[edge],posn=cum.edge[edge]-pos))
	}
	# check if tree is ultrametric (required)
	if(is.ultrametric(tree)) um<-TRUE
	else um<-FALSE
	# set n and/or name tips (if not provided)
	if(is.null(tips)){
		if(is.null(n)) n<-1
		tips<-paste("t",length(tree$tip)+1:n,sep="")
	} else n<-length(tips)
	if(is.null(edge.length)) if(!um) edge.length<-runif(n=n,min=min(tree$edge.length),max=max(tree$edge.length))
	# set order
	if(order[1]=="random"){
		o<-sample(1:n)
		tips<-tips[o]
		if(!is.null(edge.length)) edge.length<-edge.length[o]
	}
	# add tips
	for(i in 1:n){
		where<-randomPosn(tree)
		if(is.null(edge.length)) tree<-bind.tip(tree,tips[i],where=where$node,position=where$posn)
		else tree<-bind.tip(tree,tips[i],where=where$node,position=where$posn,edge.length=edge.length[i])
	}
	return(tree)
}

