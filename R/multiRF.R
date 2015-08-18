## functions computes Robinson-Foulds distance between all trees in a list of class "multiPhylo"
## works only for unrooted trees (if trees are rooted, them will be unrooted)
## written by Liam J. Revell 2013, 2015

multiRF<-function(trees){
	if(!inherits(trees,"multiPhylo")) stop("trees should be an object of class \"multiPhylo\".")
	N<-length(trees)
	RF<-matrix(0,N,N)
	if(any(sapply(unclass(trees),is.rooted))){
		cat("Some trees are rooted. Unrooting all trees.\n")
		trees<-lapply(unclass(trees),unroot)
	}
	if(any(sapply(unclass(trees),function(x) !is.binary.tree(x)))){
		stop("Some trees are not binary. This implementation only works for binary trees.")
	}
	foo<-function(pp) lapply(pp,function(x,pp) sort(attr(pp,"labels")[x]),pp=pp)
	xx<-lapply(unclass(trees),function(x) foo(prop.part(x))[2:x$Nnode])
	for(i in 1:(N-1)) for(j in (i+1):N) RF[i,j]<-RF[j,i]<-2*sum(!xx[[i]]%in%xx[[j]])
	RF
}
