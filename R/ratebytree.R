## method to compare the rate of evolution for a character between trees
## closely related to 'censored' approach of O'Meara et al. (2006; Evolution)
## written by Liam J. Revell 2017

ratebytree<-function(trees,x,...){
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-8
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-FALSE
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-6
	if(hasArg(test)) test<-list(...)$test
	else test<-"chisq"
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	## check trees & x
	if(!inherits(trees,"multiPhylo")) 
		stop("trees should be object of class \"multiPhylo\".")
	if(!is.list(x)) stop("x should be a list of vectors.")
	if(hasArg(se)) se<-list(...)$se
	else {
		se<-x
		for(i in 1:length(x)) se[[i]][1:length(se[[i]])]<-0
	}
	N<-length(trees)
	## reorder the trait vectors in x & SEs in se
	x<-mapply(function(x,t) x<-x[t$tip.label],x=x,t=trees,SIMPLIFY=FALSE)
	se<-mapply(function(x,t) x<-x[t$tip.label],x=se,t=trees,SIMPLIFY=FALSE)
	## first, fit multi-rate model
	lik.multi<-function(theta,trees,x,se,trace=FALSE){
		n<-sapply(trees,Ntip)
		N<-length(trees)
		sig<-theta[1:N]
		a<-theta[(N+1):(2*N)]
		C<-lapply(trees,vcv)
		E<-lapply(se,diag)
		V<-mapply("+",mapply("*",C,sig,SIMPLIFY=FALSE),E,SIMPLIFY=FALSE)
		logL<-0
		for(i in 1:N) 
			logL<-logL-t(x[[i]]-a[i])%*%solve(V[[i]])%*%(x[[i]]-
				a[i])/2-n[i]*log(2*pi)/2-determinant(V[[i]])$modulus[1]/2
		if(trace){
			cat(paste(paste(round(sig,digits),collapse="\t"),round(logL,digits),"\n",sep="\t"))
			flush.console()
		}
		-logL
	}
	foo<-function(tree,x){ 
		pvcv<-phyl.vcv(as.matrix(x),vcv(tree),1)
		c(pvcv$R[1,1],pvcv$a[1,1])
	}
	PP<-mapply(foo,trees,x)
	if(hasArg(init)){ 
		init<-list(...)$init
	 	if(!is.null(init$sigm)) PP[1,]<-init$sigm
		if(!is.null(init$am)) PP[2,]<-init$am
	}
	p<-as.vector(t(PP))
	if(trace){
		cat("\nOptimizing multi-rate model....\n")
		cat(paste(paste("sig[",1:N,"]   ",sep="",collapse="\t"),"logL\n",sep="\t"))
	}
	fit.multi<-optim(p,lik.multi,trees=trees,x=x,se=se,trace=trace,method="L-BFGS-B",
		lower=c(rep(tol,N),rep(-Inf,N)),upper=c(rep(Inf,N),rep(Inf,N)))
	## now fit single-rate model
	lik.onerate<-function(theta,trees,x,se,trace=FALSE){
		n<-sapply(trees,Ntip)
		N<-length(trees)
		sig<-theta[1]
		a<-theta[1:N+1]
		C<-lapply(trees,vcv)
		E<-lapply(se,diag)
		V<-mapply("+",lapply(C,"*",sig),E,SIMPLIFY=FALSE)
		logL<-0
		for(i in 1:N) 
			logL<-logL-t(x[[i]]-a[i])%*%solve(V[[i]])%*%(x[[i]]-
				a[i])/2-n[i]*log(2*pi)/2-determinant(V[[i]])$modulus[1]/2
		if(trace){ 
			cat(paste(round(sig,digits),round(logL,digits),"\n",sep="\t"))
			flush.console()
		}
		-logL
	}
	p<-c(mean(fit.multi$par[1:N]),fit.multi$par[(N+1):(2*N)])
	if(hasArg(init)){
		if(!is.null(init$sigc)) p[1]<-init$sigc
		if(!is.null(init$ac)) p[1:N+1]<-init$ac
	}
	if(trace){
		cat("\nOptimizing common-rate model....\n")
		cat(paste("sig      ","logL\n",sep="\t"))
	}
	fit.onerate<-optim(p,lik.onerate,trees=trees,x=x,se=se,trace=trace,method="L-BFGS-B",
		lower=c(tol,rep(-Inf,N)),upper=c(Inf,rep(Inf,N)))
	## compare models:
	LR<-2*(-fit.multi$value+fit.onerate$value)
	if(test=="chisq") P.chisq<-pchisq(LR,df=N-1,lower.tail=FALSE)
	else if(test=="simulation"){
		if(!quiet) cat("Generating null distribution via simulation -> |")
		flush.console()
		if(hasArg(nsim)) nsim<-list(...)$nsim
		else nsim<-100
		X<-mapply(fastBM,tree=trees,a=as.list(fit.onerate$par[1:N+1]),
			MoreArgs=list(sig2=fit.onerate$par[1],nsim=nsim),
			SIMPLIFY=FALSE)
		P.sim<-1/(nsim+1)
		pct<-0.1
		for(i in 1:nsim){
			x.sim<-lapply(X,function(x,ind) x[,ind],ind=i)
			foo<-function(x,se) sampleFrom(xbar=x,xvar=se^2,n=rep(1,length(x)))
			x.sim<-mapply(foo,x=x.sim,se=se,SIMPLIFY=FALSE)
			fit.sim<-ratebytree(trees,x.sim,se=se)
			P.sim<-P.sim+(fit.sim$likelihood.ratio>=LR)/(nsim+1)
			if(i/nsim>=pct){
				if(!quiet) cat(".")
				flush.console()
				pct<-pct+0.1
			}
		}
		if(!quiet) cat(".|\nDone!\n")
		flush.console()
	}
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
		N=N,likelihood.ratio=LR,
		P.chisq=if(test=="chisq") P.chisq else NULL,
		P.sim=if(test=="simulation") P.sim else NULL)
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
	if(!is.null(x$P.chisq))
		cat(paste("P-value (based on X^2):",round(x$P.chisq,digits),"\n\n"))
	else if(!is.null(x$P.sim))
		cat(paste("P-value (based on simulation):",round(x$P.sim,digits),"\n\n"))
	if(x$multi.rate.model$convergence==0&&x$common.rate.model$convergence==0) 
        cat("R thinks it has found the ML solution.\n\n")
	else cat("One or the other optimization may not have converged.\n\n")
}
