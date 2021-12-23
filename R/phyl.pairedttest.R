## function for phylogenetic paired t-test (Lindenfors et al. 2010)
## written by Liam Revell 2011, 2013, 2015, 2021

phyl.pairedttest<-function(tree,x1,x2=NULL,se1=NULL,se2=NULL,lambda=1.0,
	h0=0.0,fixed=FALSE,...){
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-12
	## check tree
	if(!inherits(tree,"phylo")) 
		stop("tree should be an object of class \"phylo\".")
	## convert x1 to a matrix if necessary
	if(is.data.frame(x1)) x1<-as.matrix(x1)
	## if x2 is NULL
	if(is.null(x2)){
		## check to see that x1 has two columns
		if(dim(x1)[2]!=2) 
			stop("user must provide two data vectors or matrix with two variables")
		else {
			x2<-x1[,2]
			x1<-x1[,1]
		}
	}
	if(is.data.frame(x2)) x2<-as.matrix(x2) # to be safe
	if(is.matrix(x1)) x1<-x1[,1]
	if(is.matrix(x2)) x2<-x2[,1]
	if(is.null(se1)){ 
		v1<-rep(0,length(tree$tip))
		names(v1)<-tree$tip.label
	} else v1<-se1^2
	if(is.null(se2)){
		v2<-rep(0,length(tree$tip))
		names(v2)<-tree$tip.label
	} else v2<-se2^2
	## compute C and sort vectors
	C<-vcv.phylo(tree)
	x1<-x1[rownames(C)]
	x2<-x2[rownames(C)]
	v1<-v1[rownames(C)]
	v2<-v2[rownames(C)]
	V.diff<-diag(v1+v2)
	dimnames(V.diff)<-list(names(v1),names(v1))
	## compute difference
	d<-x1-x2
	## lambda transformation
	lambda.transform<-function(C,lambda) lambda*(C-diag(diag(C)))+diag(diag(C))
	## likelihood function
	likelihood<-if(fixed) function(theta,d,C,V.diff,fixed.lambda){
		theta[1]->sig2
		theta[2]->dbar
		V<-sig2*lambda.transform(C,fixed.lambda)+V.diff
		logL<--t(d-dbar)%*%solve(V,d-dbar)/2-determinant(V)$modulus[1]/2-
			length(d)*log(2*pi)/2
		-logL[1,1]
	} else function(theta,d,C,V.diff){
		theta[1]->sig2
		theta[2]->lambda
		theta[3]->dbar
		V<-sig2*lambda.transform(C,lambda)+V.diff
		logL<--t(d-dbar)%*%solve(V,d-dbar)/2-determinant(V)$modulus[1]/2-
			length(d)*log(2*pi)/2
		-logL[1,1]
	}
	## rescale for optimization:
	dscale<-1/sqrt(mean(pic(d,multi2di(tree,random=FALSE))^2))
	d<-d*dscale
	V.diff<-V.diff*dscale^2
	## maximize the likelihood
	if(!fixed) fit<-optim(c(mean(pic(d,multi2di(tree,random=FALSE))^2),lambda,mean(d)),likelihood,
		d=d,C=C,V.diff=V.diff,method="L-BFGS-B",lower=c(tol,0,-Inf),upper=c(Inf,1,Inf),
		hessian=TRUE)
	else fit<-optim(c(mean(pic(d,multi2di(tree,random=FALSE))^2),mean(d)),likelihood,
		d=d,C=C,V.diff=V.diff,fixed.lambda=lambda,method="L-BFGS-B",lower=c(tol,-Inf),
		upper=c(Inf,Inf),hessian=TRUE)
	## run t-test
	se.dbar<-if(fixed) sqrt(1/fit$hessian[2,2]) else sqrt(1/fit$hessian[3,3])
	t<-if(fixed) (fit$par[2]-h0)/se.dbar else (fit$par[3]-h0)/se.dbar
	P<-2*pt(abs(t),df=Ntip(tree)-if(fixed) 2 else 3,lower.tail=FALSE)
	obj<-list(dbar=if(fixed) fit$par[2]/dscale else fit$par[3]/dscale,
		se=se.dbar/dscale,sig2=fit$par[1]/(dscale^2),
		lambda=if(fixed) lambda else fit$par[2],
		logL=if(fixed) -likelihood(c(fit$par[1]/(dscale^2),fit$par[2]/dscale),
		d/dscale,C,V.diff/(dscale^2),lambda) else -likelihood(c(fit$par[1]/(dscale^2),
		fit$par[2],fit$par[3]/dscale),d/dscale,C,V.diff/(dscale^2)),t.dbar=t,
		P.dbar=P,df=Ntip(tree)-if(fixed) 2 else 3,h0=h0)
	class(obj)<-"phyl.pairedttest"
	obj
}

print.phyl.pairedttest<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-6
	rnd<-function(x,digits){ 
		x<-if(abs(x)>1e5) format(x,scientific=TRUE,digits=digits) else
			if(abs(x)<=1e5&&abs(x)>1e-5) format(x,digits=digits) else
			if(x<=1e-5&&x!=0) format(x,scientific=TRUE,digits=digits) else
			format(x,digits=digits)
			x
	}
	cat("\nPhylogenetic paired t-test:\n\n")
	cat(paste("   t = ",rnd(x$t.dbar,digits),", df = ",x$df,
		", p-value = ",rnd(x$P.dbar,digits),sep=""))
	cat(paste(
		"\n\nalternative hypothesis:\n   true difference in means is not equal to",
		rnd(x$h0,digits)))
	cat("\n\n95 percent confidence interval on the phylogenetic\ndifference in mean:\n")
	cat(paste("   [",paste(sapply(x$dbar+c(-1.96,1.96)*x$se,rnd,digits=digits),
		collapse=", "),"]\n",sep=""))
	cat("\nestimates:\n")
	cat("   phylogenetic mean difference =",rnd(x$dbar,digits),"\n")
	cat("   sig^2 of BM model =",rnd(x$sig2,digits),"\n")
	cat("   lambda (fixed or estimated) =",rnd(x$lambda,digits),"\n\n")
	cat("log-likelihood:\n")
	cat(paste("   ",rnd(x$logL,digits),"\n\n"))
}
