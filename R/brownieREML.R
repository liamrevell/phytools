# This function is a simplified REML version of brownie.lite()
# written by Liam J. Revell 2011, 2013

brownieREML<-function(tree,x,maxit=2000,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	## optional arguments
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-8
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
	H<-optimHess(sig1,lik1,tree=tree,x=x)
	v1<-1/H
	# fit the multiple rate model
	lik2<-function(sig2,tree,x){
		tt<-scaleByMap(tree,sig2)
		picX<-pic(x,tt,scaled=F,var.contrasts=T)
		logL<-sum(dnorm(picX[,1],sd=sqrt(picX[,2]),log=TRUE))
		return(-logL)
	}
	YY<-optim(setNames(rep(1,p)*runif(n=p),colnames(tree$mapped.edge)),lik2,tree=tree,x=x,method="L-BFGS-B",lower=rep(1e-8,p))
	sig2<-YY$par
	obj<-optimHess(sig2,lik2,tree=tree,x=x)
	if(any(diag(obj)<tol)){ 
		ii<-which(diag(obj)>0)
		H<-obj[ii,ii]
		v2<-matrix(Inf,nrow(obj),ncol(obj))
		v2[ii,ii]<-if(length(H)>1) solve(H) else 1/H
	} else v2<-if(length(sig2)>1) solve(obj) else 1/obj
	logL2<--YY$value
	if(YY$convergence==0) converged<-"Optimization has converged."
	else converged<-"Optimization may not have converged.  Consider increasing maxit."
	convergence=(YY$convergence==0)
	obj<-list(sig2.single=sig1,var.single=v1,logL1=logL1,k1=1,sig2.multiple=sig2,vcv.multiple=v2,logL.multiple=logL2,
		k2=length(sig2),convergence=converged)
	class(obj)<-"brownieREML"
	obj
}

# This function scales a mapped tree by sig2
# written by Liam J. Revell 2011
scaleByMap<-function(mtree,sig2){
	edge.length<-if(length(sig2)>1) mtree$mapped.edge[,names(sig2)]%*%sig2 else mtree$mapped.edge*sig2
	tree<-list(Nnode=mtree$Nnode,edge=mtree$edge,tip.label=mtree$tip.label,edge.length=edge.length[,1])
	class(tree)<-"phylo"
	return(tree)
}

## S3 print method for "brownieREML"
## S3 print method for object of class "brownie.lite"
## written by Liam J. Revell 2013
print.brownieREML<-function(x, ...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	x<-lapply(x,function(a,b) if(is.numeric(a)) round(a,b) else a,b=digits)
	cat("REML single-rate model:\n")
	cat("\ts^2\tse\tk\tlogL\n")
	cat(paste("value",x$sig2.single,round(sqrt(x$var.single),digits),x$k1,x$logL1,"\n\n",sep="\t"))
	cat("REML multi-rate model:\n")
	cat(paste(c("",paste("s^2(",names(x$sig2.multiple),")","\tse(",names(x$sig2.multiple),")",sep=""),
		"k","logL","\n"),collapse="\t"))
	cat(paste(paste(c("value",paste(x$sig2.multiple,round(sqrt(diag(x$vcv.multiple)),digits),sep="\t"),x$k2,
		x$logL.multiple),collapse="\t"),"\n\n",sep=""))
	if(x$convergence[1]=="Optimization has converged.") cat("R thinks it has found the ML solution.\n\n")
	else cat("Optimization may not have converged.\n\n")
}
