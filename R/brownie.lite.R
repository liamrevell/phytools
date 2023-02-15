## this function fits two or more evolutionary rates for a continuous trait on the tree
## based on O'Meara et al. (2006)
## written by Liam J. Revell 2011/2012, 2019, 2022

brownie.lite<-function(tree,x,maxit=2000,test="chisq",nsim=100,se=NULL,...){
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	# some minor error checking
	if(!inherits(tree,"phylo")) 
		stop("tree should be an object of class \"phylo\".")
	if(!inherits(tree,"simmap")) 
		tree<-paintSubTree(tree,Ntip(tree)+1,"1")
	x<-matchDatatoTree(tree,x,"x")
	x<-x[tree$tip.label]
	if(!is.null(se)){ 
		se<-matchDatatoTree(tree,se,"se")
		se<-se[tree$tip.label]
	} else {
		se<-rep(0,length(x))
		names(se)<-names(x)
	}
	n<-length(x) # number of species
	p<-ncol(tree$mapped.edge) # number of states
	C1<-vcv.phylo(tree)
	a<-as.numeric(colSums(solve(C1))%*%x/sum(solve(C1)))
	sig<-as.numeric(t(x-a)%*%solve(C1)%*%(x-a)/n)
	# single rate model
	model1<-optim(c(sig,a),fn=lik.single,y=x,C=C1,se=se,
		control=list(maxit=maxit),hessian=TRUE,
		method="L-BFGS-B",
		lower=c(0,-Inf))
	logL1<--model1$value
	sig1<-model1$par[1]
	a1<-model1$par[2]
	vcv1<-solve(model1$hessian)
	rownames(vcv1)<-c("sig","a")
	colnames(vcv1)<-rownames(vcv1)
	# multiple rate model
	C2<-multiC(tree)
	s<-c(rep(sig1,p),a1)
	l<-c(rep(0.0001*sig1,p),-Inf)
	model2<-optim(s,fn=lik.multiple,y=x,C=C2,se=se,
		control=list(maxit=maxit),hessian=TRUE,
		method="L-BFGS-B",lower=l)	
	logL2<--model2$value
	while(logL2<logL1){
		if(!quiet){ 
			message("False convergence on first try; trying again with new starting values.")
		} 
		model2<-optim(s*2*runif(n=length(s)),fn=lik.multiple,
			y=x,C=C2,se=se,control=list(maxit=maxit),
			hessian=TRUE,method="L-BFGS-B",lower=l)
		logL2<--model2$value
	}
	sig.i<-model2$par[1:p]
	names(sig.i)<-colnames(tree$mapped.edge)
	a2<-model2$par[p+1]
	vcv2<-solve(model2$hessian)
	rownames(vcv2)<-c(colnames(tree$mapped.edge),"a")
	colnames(vcv2)<-rownames(vcv2)
	if(model2$convergence==0) converged<-"Optimization has converged."
	else converged<-"Optimization may not have converged.  Consider increasing maxit."
	if(test=="chisq")
		xx<-list(sig2.single=sig1,a.single=a1,var.single=vcv1[1,1],
			logL1=logL1,k1=2,sig2.multiple=sig.i,a.multiple=a2,
			vcv.multiple=vcv2[1:p,1:p],logL.multiple=logL2,k2=p+1,
			P.chisq=pchisq(2*(logL2-as.numeric(logL1)),p-1,
			lower.tail=FALSE),convergence=converged)
	else if(test=="simulation"){
		LR<-2*(logL2-logL1)
		X<-fastBM(tree,a=a1,sig2=sig1,nsim=(nsim-1))
		Xe<-matrix(rnorm(n=length(X),mean=X,sd=rep(se,nsim-1)),
			nrow(X),ncol(X),dimnames=dimnames(X))
		# now compute the P-value based on simulation
		Psim<-1/nsim
		for(i in 1:(nsim-1)){
			sim<-brownie.lite(tree,Xe[,i],se=se)
			while(sim$convergence!="Optimization has converged."||
				(sim$logL.multiple-sim$logL1)<0){
				Xe[,i]<-rnorm(n=length(x),fastBM(tree,a=a1,sig2=sig1),sd=se)
				sim<-brownie.lite(tree,Xe[,i])
			}
			Psim<-Psim+(LR<=2*(sim$logL.multiple-sim$logL1))/nsim
		}
		xx<-list(sig2.single=sig1,a.single=a1,var.single=vcv1[1,1],logL1=logL1,
			k1=2,sig2.multiple=sig.i,a.multiple=a2,vcv.multiple=vcv2[1:p,1:p],
			logL.multiple=logL2,k2=p+1,P.sim=Psim,convergence=converged)
	}
	class(xx)<-"brownie.lite"
	return(xx)
}

## S3 print method for object of class "brownie.lite"
## written by Liam J. Revell 2013, 2020

print.brownie.lite<-function(x, ...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-getOption("digits")
	x<-lapply(x,function(a,b) if(is.numeric(a)) round(a,b) else a,b=digits)
	cat("ML single-rate model:\n")
	obj<-matrix(c(x$sig2.single,sqrt(x$var.single),
		x$a.single,x$k1,x$logL1),1,5,
		dimnames=list("value",c("s^2","se","a","k","logL")))
	print(obj,digits=digits)
	cat("\nML multi-rate model:\n")
	nn<-c(unlist(strsplit(paste("s^2(",names(x$sig2.multiple),")__",
		"se(",names(x$sig2.multiple),")",sep=""),"__")),"a","k",
		"logL")
	obj<-matrix(c(as.vector(rbind(x$sig2.multiple,sqrt(diag(x$vcv.multiple)))),
		x$a.multiple,x$k2,x$logL.multiple),1,2*length(x$sig2.multiple)+3,
		dimnames=list("value",nn))
	print(obj,digits=digits)
	if(!is.null(x$P.chisq)) cat(paste("\nP-value (based on X^2):",x$P.chisq,"\n\n"))
	else if(!is.null(x$P.sim)) cat(paste("\nP-value (based on simulation):",
		x$P.sim,"\n\n"))
	if(x$convergence[1]=="Optimization has converged.") 
		cat("R thinks it has found the ML solution.\n\n")
	else cat("Optimization may not have converged.\n\n")
}

# function computes the likelihood for a single rate with sampling error
# written by Liam J. Revell 2012

lik.single<-function(theta,y,C,se){
	n<-length(y)
	sig<-theta[1]
	a<-theta[2]
	E<-diag(se^2)
	logL<-as.numeric(-t(y-a)%*%solve(sig*C+E)%*%(y-a)/2-
		n*log(2*pi)/2-determinant(sig*C+E)$modulus[1]/2)
	return(-logL)
}

# function computes the likelihood for multiple rates
# written by Liam J. Revell 2012

lik.multiple<-function(theta,y,C,se=NULL){
	n<-length(y); p<-length(C)
	sig<-theta[1:p]
	a<-theta[p+1]
	V<-matrix(0,length(y),length(y))
	for(i in 1:p) V<-V+sig[i]*C[[i]]
	E<-diag(se^2)
	logL<--t(y-a)%*%solve(V+E)%*%(y-a)/2-n*log(2*pi)/2-
		determinant(V+E)$modulus[1]/2
	return(-logL)
}

## S3 logLik method for object class

logLik.brownie.lite<-function(object,...){
	lik<-setNames(
		c(object$logL1,object$logL.multiple),
		c("single-rate","multi-rate"))
	attr(lik,"df")<-c(object$k1,object$k2)
	lik
}
