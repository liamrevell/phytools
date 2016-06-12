ls.consensus<-function(trees,start=NULL,tol=1e-12,quiet=FALSE,...){
	D<-Reduce("+",lapply(trees,function(x,t) cophenetic(x)[t,t], 
		t=trees[[1]]$tip.label))/length(trees)
	if(is.null(start)) start<-NJ(D)
	if(hasArg(ultrametric)) ultrametric<-list(...)$ultrametric ## should the consensus tree be ultrametric
	else ultrametric<-all(sapply(trees,is.ultrametric))
	if(ultrametric&&!is.rooted(start)) start<-midpoint(start)
	curr<-nnls.tree(D,tree=start,rooted=ultrametric,trace=0)
	Q<-Inf
	Qp<-attr(curr,"RSS")
	if(is.null(Qp)) Qp<-rss(D,curr)
	ct<-0
	while((Q-Qp)>tol){
		Q<-Qp
		NNIs<-.uncompressTipLabel(nni(curr))
		curr<-list(curr)
		class(curr)<-"multiPhylo"
		NNIs<-c(NNIs,curr)
		NNIs<-lapply(NNIs,nnls.tree,dm=D,rooted=ultrametric,trace=0)
		qs<-sapply(NNIs,rss,D=D)
		ii<-which(qs==min(qs))[1]
		if(!quiet) message(paste("Best Q =",qs[ii]))
		Qp<-qs[ii]
		curr<-NNIs[[ii]]
		ct<-ct+1
	}
	if(!quiet) 
		message(paste("Solution found after",ct,
			"set of nearest neighbor interchanges."))
	curr
}

rss<-function(D,tree){
	Dp<-cophenetic(tree)[rownames(D),colnames(D)]
	sum((D-Dp)^2)/2
}

minTreeDist<-function(tree,trees,method="quadratic.path.difference",...){
	if(is.null(tree$edge.length)) 
		tree$edge.length<-runif(n=nrow(tree$edge))
	e<-tree$edge.length
	fn<-function(e,tree,trees,m){
		tree$edge.length<-e
		if(method=="quadratic.path.difference")
			d<-qpd(tree,trees)
		else
			d<-sapply(trees,function(x,y) treedist(x,y)[m],
				y=tree)
		sum(d^2)
	}
	fit<-optim(e,fn,tree=tree,trees=trees,m=method,lower=0,
		method="L-BFGS-B",...)
	tree$edge.length<-fit$par
	attr(tree,"SQD")<-fit$value
	tree
}

averageTree<-function(trees,start=NULL,method="quadratic.path.difference",
	tol=1e-12,quiet=FALSE,...){
	if(!quiet)
		message(paste(
			"\n  Function is attempting to find the phylogeny with ",
			"\n  minimum distance to all trees in set under the ",
			"\n  \"",
			paste(strsplit(method,"[.]")[[1]],collapse=" "),
			"\" criterion...\n",sep=""))
	if(is.null(start)){
		if(method%in%c("symmetric.difference","path.difference"))
			start<-multi2di(consensus(trees,p=0.5))
		else start<-ls.consensus(trees,quiet=TRUE)
	}
	if(method%in%c("branch.score.difference","quadratic.path.difference")){
		D<-Reduce("+",lapply(trees,function(x,t) cophenetic(x)[t,t], 
			t=trees[[1]]$tip.label))/length(trees)
		rt<-all(sapply(trees,is.ultrametric))
	} else if(method%in%c("symmetric.difference","path.difference")){
		rt<-all(sapply(trees,is.rooted))
		if(!rt) start<-unroot(start)
	}
	curr<-start
	SS<-Inf
	SSp<-sum(sapply(trees,function(x,y,m) treedist(x,y)[m],
		y=curr,m=method)^2)
	ct<-0
	while((SS-SSp)>tol){
		SS<-SSp
		NNIs<-.uncompressTipLabel(nni(curr))
		curr<-list(curr)
		class(curr)<-"multiPhylo"
		NNIs<-c(NNIs,curr)
		if(method=="symmetric.difference")
			SSp<-colSums(sapply(NNIs,RF.dist,tree2=trees)^2)
		else if(method=="path.difference")
			SSp<-colSums(sapply(NNIs,function(x,y,m) sapply(y,
				function(y,x,m) setNames(treedist(y,x)[m],NULL),
				y=x,m=m),y=trees,m=method)^2)
		else {
			if(!rt) NNIs<-lapply(NNIs,unroot)
			NNIs<-lapply(NNIs,nnls.tree,dm=D,rooted=rt,trace=0)
			NNIs<-lapply(NNIs,minTreeDist,trees=trees,method=method,
				...)
			SSp<-sapply(NNIs,function(x) attr(x,"SQD"))
		}
		ii<-which(SSp==min(SSp))[1]
		if(!quiet) message(paste("  Best SS so far =",SSp[ii]))
		SSp<-SSp[ii]
		curr<-NNIs[[ii]]
		ct<-ct+1
	}
	if(!quiet) message(paste("\n  Solution found after",ct,
			"set of nearest neighbor interchanges.\n"))
	curr
}

qpd<-function(t1,t2){
	if(inherits(t1,"phylo")) t1<-list(t1)
	if(inherits(t2,"phylo")) t2<-list(t2)
	tips<-t1[[1]]$tip.label
	P1<-lapply(t1,function(x) cophenetic(x)[tips,tips])
	P2<-lapply(t2,function(x) cophenetic(x)[tips,tips])
	foo<-function(x,y) sqrt(sum((x-y)^2/2))
	sapply(P2,function(x,y) sapply(y,foo,y=x),y=P1)
}
