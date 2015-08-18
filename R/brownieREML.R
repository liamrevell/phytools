# This function is a simplified REML version of brownie.lite()
# written by Liam J. Revell 2011, 2013

brownieREML<-function(tree,x,maxit=2000){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	# bookkeeping
	if(!is.binary.tree(tree)) tree<-multi2di(tree)
	x<-x[tree$tip.label] # order in tip.label order
	n<-length(x) # number of species
	p<-ncol(tree$mapped.edge) # number of states
	# fit the single rate model
	lik1<-function(sig1,tree,x){
		tt<-scaleByMap(tree,setNames(rep(sig1,p),colnames(tree$mapped.edge)))
		picX<-pic(x,tt,scaled=FALSE,var.contrasts=TRUE)
		logL<-sum(dnorm(picX[,1],sd=sqrt(picX[,2]),log=TRUE))
		return(-logL)
	}
	sig1<-mean(pic(x,tree)^2)
	logL1<--lik1(sig1,tree,x)
	# fit the multiple rate model
	lik2<-function(sig2,tree,x){
		tt<-scaleByMap(tree,sig2)
		picX<-pic(x,tt,scaled=F,var.contrasts=T)
		logL<-sum(dnorm(picX[,1],sd=sqrt(picX[,2]),log=TRUE))
		return(-logL)
	}
	YY<-optim(setNames(rep(1,p)*runif(n=p),colnames(tree$mapped.edge)),lik2,tree=tree,x=x,method="L-BFGS-B",lower=rep(1e-8,p))
	sig2<-YY$par
	logL2<--YY$value
	convergence=(YY$convergence==0)
	return(list(sig2.single=sig1,logL1=logL1,sig2.multiple=sig2,logL2=logL2,convergence=convergence))
}

# This function scales a mapped tree by sig2
# written by Liam J. Revell 2011
scaleByMap<-function(mtree,sig2){
	edge.length<-mtree$mapped.edge[,names(sig2)]%*%sig2
	tree<-list(Nnode=mtree$Nnode,edge=mtree$edge,tip.label=mtree$tip.label,edge.length=edge.length[,1])
	class(tree)<-"phylo"
	return(tree)
}
