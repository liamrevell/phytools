## method to compare the rate of evolution for a character between trees
## closely related to 'censored' approach of O'Meara et al. (2006; Evolution)
## written by Liam J. Revell 2017

ratebytree<-function(trees,x,...){
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-8
	## check trees & x
	if(!inherits(trees,"multiPhylo")) 
		stop("trees should be object of class \"multiPhylo\".")
	if(!is.list(x)) stop("x should be a list of vectors.")
	N<-length(trees)
	## reorder the trait vectors in x
	x<-mapply(function(x,t) x<-x[t$tip.label],x=x,t=trees,SIMPLIFY=FALSE)
	## first, fit multi-rate model
	lik.multi<-function(theta,trees,x){
		n<-sapply(trees,Ntip)
		N<-length(trees)
		sig<-theta[1:N]
		a<-theta[(N+1):(2*N)]
		C<-lapply(trees,vcv)
		V<-mapply("*",C,sig,SIMPLIFY=FALSE)
		logL<-0
		for(i in 1:N) 
			logL<-logL-t(x[[i]]-a[i])%*%solve(V[[i]])%*%(x[[i]]-
				a[i])/2-n[i]*log(2*pi)/2-determinant(V[[i]])$modulus[1]/2
		-logL
	}
	foo<-function(tree,x){ 
		pvcv<-phyl.vcv(as.matrix(x),vcv(tree),1)
		c(pvcv$R[1,1],pvcv$a[1,1])
	}
	PP<-mapply(foo,trees,x)
	p<-as.vector(t(PP))
	fit.multi<-optim(p,lik.multi,trees=trees,x=x,method="L-BFGS-B",
		lower=c(rep(tol,N),rep(-Inf,N)),upper=c(rep(Inf,N),rep(Inf,N)))
	## now fit single-rate model
	lik.onerate<-function(theta,trees,x){
		n<-sapply(trees,Ntip)
		N<-length(trees)
		sig<-theta[1]
		a<-theta[1:N+1]
		C<-lapply(trees,vcv)
		V<-lapply(C,"*",sig)
		logL<-0
		for(i in 1:N) 
			logL<-logL-t(x[[i]]-a[i])%*%solve(V[[i]])%*%(x[[i]]-
				a[i])/2-n[i]*log(2*pi)/2-determinant(V[[i]])$modulus[1]/2
		-logL
	}
	p<-c(mean(fit.multi$par[1:N]),fit.multi$par[(N+1):(2*N)])
	fit.onerate<-optim(p,lik.onerate,trees=trees,x=x,method="L-BFGS-B",
		lower=c(tol,rep(-Inf,N)),upper=c(Inf,rep(Inf,N)))
	## compare models:
	LR<-2*(-fit.multi$value+fit.onerate$value)
	P.chisq<-pchisq(LR,df=N-1,lower.tail=FALSE)
	obj<-list(
		multi.rate.model=list(sig2=fit.multi$par[1:N],
			a=fit.multi$par[(N+1):(2*N)],
			k=2*N,
			logL=-fit.multi$value,
			counts=fit.multi$counts,convergence=fit.multi$convergence,
			message=fit.multi$message),
		common.rate.model=list(sig2=fit.onerate$par[1],
			a=fit.onerate$par[1:N+1],
			k=N+1,
			logL=-fit.onerate$value,
			counts=fit.onerate$counts,convergence=fit.onerate$convergence,
			message=fit.onerate$message),
		N=N,likelihood.ratio=LR,P.chisq=P.chisq)
	class(obj)<-"ratebytree"
	obj
}

## S3 print method for ratebytree
print.ratebytree<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	N<-x$N
	cat("ML common-rate model:\n")
	cat(paste("\ts^2\t",paste(paste("a[",1:N,"]",sep=""),collapse="\t")),
		"\tk\tlogL\n")
   	cat(paste("value",round(x$common.rate.model$sig2,digits),
		paste(round(x$common.rate.model$a,digits),collapse="\t"),
		x$common.rate.model$k,round(x$common.rate.model$logL,digits),
		"\n\n",sep="\t"))
	cat("ML multi-rate model:\n")
	cat(paste("\t",paste(paste("s^2[",1:N,"]",sep=""),collapse="\t"),"\t",
		paste(paste("a[",1:N,"]",sep=""),collapse="\t")),
		"\tk\tlogL\n")
   	cat(paste("value",paste(round(x$multi.rate.model$sig2,digits),collapse="\t"),
		paste(round(x$multi.rate.model$a,digits),collapse="\t"),
		x$multi.rate.model$k,round(x$multi.rate.model$logL,digits),
		"\n\n",sep="\t"))
	cat(paste("Likelihood ratio:",round(x$likelihood.ratio,digits),"\n"))
	cat(paste("P-value (based on X^2):",round(x$P.chisq,digits),"\n\n"))
	if(x$multi.rate.model$convergence==0&&x$common.rate.model$convergence==0) 
        cat("R thinks it has found the ML solution.\n\n")
	else cat("One or the other optimization may not have converged.\n\n")
}
