# function performs least-squares phylogeny inference by nni
# written by Liam J. Revell 2011, 2013, 2015

optim.phylo.ls<-function(D,stree=NULL,set.neg.to.zero=TRUE,fixed=FALSE,tol=1e-10,collapse=TRUE){

	# change D to a matrix (if actually an object of class "dist")
	if(class(D)=="dist") D<-as.matrix(D)

	# compute the number of species
	n<-nrow(D)

	if(is.null(stree))
		stree<-rtree(n=n,tip.label=rownames(D),br=NULL,rooted=F) # random starting tree
	else if(!inherits(stree,"phylo")){
		cat("starting tree must be an object of class \"phylo.\" using random starting tree.\n")
		stree<-rtree(n=n,tip.label=rownames(D),br=NULL,rooted=F) # random starting tree
	}
	if(!is.binary.tree(stree)) stree<-multi2di(stree)
	if(is.rooted(stree)) stree<-unroot(stree)

	# get ls branch lengths for stree
	best.tree<-ls.tree(stree,D)
	Q<-attr(best.tree,"Q-score")
	bestQ<-0 # to start the loop

	# for search
	Nnni<-0

	# loop while Q is not improved by nni
	while(bestQ-Q<tol&&fixed==FALSE){

		nni.trees<-nni(best.tree)
		nniQ<-vector()

		bestQ<-Inf
		for(i in 1:length(nni.trees)){

			# compute least squares branch lengths and Q
			nni.trees[[i]]<-ls.tree(nni.trees[[i]],D)

			# compute Q
			nniQ[i]<-attr(nni.trees[[i]],"Q-score")

			# is this the best one so far?
			if(nniQ[i]<bestQ){ 
				bestQ<-nniQ[i]
				ind<-i
			}
		}

		# set new best tree
		if(bestQ<Q){
			best.tree<-nni.trees[[ind]]
			Nnni<-Nnni+1
			Q<-attr(best.tree,"Q-score")
			cat(paste(Nnni,"set(s) of nearest neighbor interchanges. best Q so far =",round(Q,10),"\n",collapse=""))
		} else bestQ<-Inf
	}

	cat(paste("best Q score of",round(Q,10),"found after",Nnni,"nearest neighbor interchange(s).\n",collapse=""))

	if(set.neg.to.zero) best.tree$edge.length[best.tree$edge.length<0]<-0

	attr(best.tree,"Q-score")<-Q
	if(collapse) best.tree<-di2multi(best.tree)
	return(best.tree)

}

# function computes the ls branch lengths and Q score for a tree
# written by Liam J. Revell 2011

ls.tree<-function(tree,D){

	# compute design matrix for tree i
	X<-phyloDesign(tree)

	# sort and columnarize D
	D<-D[tree$tip.label,tree$tip.label]
	colD<-D[lower.tri(D)]

	# compute the least squares branches conditioned on tree i
	v<-solve(t(X)%*%X)%*%t(X)%*%colD

	# give the tree its estimated branch lengths
	tree$edge.length<-v
		
	# compute the distances for this tree
	d<-X%*%v
	
	# compute Q
	Q<-sum((colD-d)^2)

	# assign attribute to tree
	attr(tree,"Q-score")<-Q

	return(tree)

}

# function computes design matrix for least squares given a topology
# written by Liam J. Revell 2011

phyloDesign<-function(tree){

	n<-length(tree$tip)

	cnames<-vector(); for(i in 1:nrow(tree$edge)) cnames[i]<-paste(as.character(tree$edge[i,]),collapse=",")
	k<-1; rnames<-vector(); for(i in 1:n) for(j in 1:n) if(j>i) { rnames[k]<-paste(c(i,j),collapse=","); k<-k+1 }

	X<-matrix(0,n*(n-1)/2,nrow(tree$edge),dimnames=list(rnames,cnames))

	anc.nodes<-compute.ancestor.nodes(tree)

	anc.nodes<-anc.nodes[,c(n+1:tree$Nnode,1:n)]

	for(i in 1:(n-1)){
		for(j in (i+1):n){
			nodes.i<-names(anc.nodes[i,anc.nodes[i,]==1])
			for(k in 1:length(nodes.i)) 
				X[paste(c(i,j),collapse=","),match(as.numeric(nodes.i[k]),tree$edge[,2])]<-X[paste(c(i,j),collapse=","),match(as.numeric(nodes.i[k]),tree$edge[,2])]+1
			nodes.j<-names(anc.nodes[j,anc.nodes[j,]==1])
			for(k in 1:length(nodes.j)) 
				X[paste(c(i,j),collapse=","),match(as.numeric(nodes.j[k]),tree$edge[,2])]<-X[paste(c(i,j),collapse=","),match(as.numeric(nodes.j[k]),tree$edge[,2])]+1
		}
	}
	X[X>1]<-0

	return(X)
}

# function computes the ancestor node numbers for each tip number
# written by Liam J. Revell 2011

compute.ancestor.nodes<-function(tree){
	n<-length(tree$tip)
	m<-tree$Nnode
	X<-matrix(0,n,n+m,dimnames=list(1:n,1:(n+m)))
	for(i in 1:n){
		currnode<-i
		while(currnode!=(n+1)){
			X[i,currnode]<-1
			currnode<-tree$edge[match(currnode,tree$edge[,2]),1]
		}
		X[i,currnode]<-1
	}
	return(X)
}
