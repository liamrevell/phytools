## compute consensus edge lengths from a set of trees given (or not) a consensus topology
## written by Liam J. Revell 2016

consensus.edges<-function(trees,method=c("mean.edge","least.squares"),...){
	if(hasArg(consensus.tree)) consensus.tree<-list(...)$consensus.tree
	else consensus.tree<-consensus(trees,p=0.5)
	tree<-consensus.tree ## get rid of this cumbersome name
	if(hasArg(if.absent)) if.absent<-list(...)$if.absent
	else if.absent<-"zero"
	N<-length(trees)
	if(method[1]=="mean.edge"){
		M<-lapply(trees,function(x,y) rbind(matchLabels(y,x),matchNodes(y,x)),y=tree)
		nodes<-M[[1]][,1]
		edge.length<-vector(mode="numeric",length=length(nodes))
		for(i in 1:length(nodes)){
			ii<-which(tree$edge[,2]==nodes[i])
			n.absent<-0
			for(j in 1:N){
				edge.length[ii]<-edge.length[ii] +
					if(!is.na(M[[j]][i,2])) trees[[j]]$edge.length[which(trees[[j]]$edge[,2]==M[[j]][i,2])]/N
					else 0
				if(is.na(M[[j]][i,2])) n.absent<-n.absent+1
			}
			if(if.absent=="ignore") edge.length[ii]<-edge.length[ii]*N/(N-n.absent)
		}
		tree$edge.length<-edge.length
	} else if(method[1]=="least.squares"){
		D<-Reduce('+',lapply(trees,function(x,t) cophenetic(x)[t,t],t=tree$tip.label))/N
		tree<-nnls.tree(D,tree=tree,rooted=all(sapply(trees,is.ultrametric)))
	}
	tree
}