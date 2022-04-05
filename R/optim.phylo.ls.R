## function performs least-squares phylogeny inference by nni
## written by Liam J. Revell 2011, 2013, 2015, 2019, 2022

optim.phylo.ls<-function(D,stree=NULL,set.neg.to.zero=TRUE,fixed=FALSE,tol=1e-10,collapse=TRUE){
	# change D to a matrix (if actually an object of class "dist")
	if(inherits(D,"dist")) D<-as.matrix(D)
	# compute the number of species
	n<-nrow(D)
	if(is.null(stree))
		stree<-rtree(n=n,tip.label=rownames(D),br=NULL,rooted=F) # random starting tree
	else if(!inherits(stree,"phylo")){
		cat("starting tree must be an object of class \"phylo.\" using random starting tree.\n")
		stree<-rtree(n=n,tip.label=rownames(D),br=NULL,rooted=F) # random starting tree
	}
	if(!is.binary(stree)) stree<-multi2di(stree)
	if(is.rooted(stree)) stree<-unroot(stree)
	# get ls branch lengths for stree
	best.tree<-ls.tree(stree,D)
	Q<-attr(best.tree,"Q-score")
	bestQ<-0 # to start the loop
	# for search
	Nnni<-0
	# loop while Q is not improved by nni
	while(bestQ-Q<tol&&fixed==FALSE){
		nni.trees<-lapply(nni(best.tree),ls.tree,D=D)
		nniQ<-sapply(nni.trees,function(x) attr(x,"Q-score"))
		ii<-which(nniQ==min(nniQ))
		bestQ<-nniQ[ii]
		if(bestQ<Q){
			best.tree<-nni.trees[[ii]]
			Nnni<-Nnni+1
			Q<-attr(best.tree,"Q-score")
			cat(paste(Nnni,"set(s) of nearest neighbor interchanges. best Q so far =",round(Q,10),"\n",collapse=""))
			flush.console()
		} else bestQ<-Inf
	}
	cat(paste("best Q score of",round(Q,10),"found after",Nnni,"nearest neighbor interchange(s).\n",collapse=""))
	if(set.neg.to.zero) best.tree$edge.length[best.tree$edge.length<0]<-0
	attr(best.tree,"Q-score")<-Q
	if(collapse) best.tree<-di2multi(best.tree)
	best.tree
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
	tree
}

# function computes design matrix for least squares given a topology
# written by Liam J. Revell 2011, totally re-written 2015
phyloDesign<-function(tree){
	N<-Ntip(tree)
	A<-lapply(1:N,function(n,t) c(getAncestors(t,n),n),t=tree)
	X<-matrix(0,N*(N-1)/2,nrow(tree$edge))
	colnames(X)<-apply(tree$edge,1,paste,collapse=",")
	rn<-sapply(1:N,function(x,y) sapply(y,paste,x=x,sep=","),y=1:N)
	rownames(X)<-rn[upper.tri(rn)]
	ii<-1
	for(i in 1:(N-1)) for(j in (i+1):N){ 
		e<-c(setdiff(A[[i]],A[[j]]),setdiff(A[[j]],A[[i]]))
		e<-sapply(e,function(x,y) which(y==x),y=tree$edge[,2])
		X[ii,e]<-1
		ii<-ii+1	
	}
	X
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
	X
}
