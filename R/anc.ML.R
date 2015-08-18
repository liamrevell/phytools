## lightweight version of ace(...,method="ML") for continuous traits
## also allows missing data in x, in which case missing data are also estimated
## written by Liam J. Revell 2011, 2013, 2014, 2015

anc.ML<-function(tree,x,maxit=2000,model=c("BM","OU"),...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(model[1]=="BM") X<-anc.BM(tree,x,maxit,...)
	else if(model[1]=="OU") X<-anc.OU(tree,x,maxit,...)
	else stop(paste("Do not recognize method",model))
	X
}

## internal to estimate ancestral states under a BM model
## written by Liam J. Revell 2011, 2013, 2014
anc.BM<-function(tree,x,maxit,...){
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-FALSE
	## check to see if any tips are missing data
	xx<-setdiff(tree$tip.label,names(x))
	## function returns the log-likelihood
	likelihood<-function(par,C,invC,detC,x,msp,trace){
		sig2<-par[1]
		a<-par[2]
		y<-par[1:(tree$Nnode-1)+2]
		x<-c(x,setNames(par[1:length(msp)+tree$Nnode+1],msp))
		x<-x[rownames(C)[1:length(tree$tip.label)]]
		z<-c(x,y)-a
		logLik<-(-z%*%invC%*%z/(2*sig2)-nrow(C)*log(2*pi)/2-nrow(C)*log(sig2)/2-
			detC/2)[1,1]
		if(trace) print(c(sig2,logLik))
		-logLik
	}
	## compute C
	C<-vcvPhylo(tree)
	invC<-solve(C)
	detC<-determinant(C,logarithm=TRUE)$modulus[1]
	## assign starting values
	zz<-fastAnc(tree,c(x,setNames(rep(mean(x),length(xx)),xx)))
	y<-zz[2:length(zz)]
	a<-zz[1]
	bb<-c(c(x,setNames(rep(mean(x),length(xx)),xx))[tree$tip.label],y)
	sig2<-((bb-a)%*%invC%*%(bb-a)/nrow(C))[1,1]
	fit<-optim(c(sig2,a,y,rep(mean(x),length(xx))),fn=likelihood,C=C,invC=invC,
		detC=detC,x=x,msp=xx,trace=trace,method="L-BFGS-B",
		lower=c(10*.Machine$double.eps,rep(-Inf,tree$Nnode+length(xx))),
		control=list(maxit=maxit))
	states<-fit$par[1:tree$Nnode+1]
	names(states)<-c(length(tree$tip)+1,rownames(C)[(length(tree$tip)+1):nrow(C)])
	X<-list(sig2=fit$par[1],ace=states,logLik=-fit$value,counts=fit$counts,
		convergence=fit$convergence,message=fit$message,model="BM")
	if(length(xx)>0) X$missing.x<-setNames(fit$par[1:length(xx)+tree$Nnode+1],xx)
	class(X)<-"anc.ML"
	X
}

## internal to estimate ancestral states under an OU model (this may not work)
## written by Liam J. Revell 2014
anc.OU<-function(tree,x,maxit=2000,...){
	## check to see if any tips are missing data
	xx<-setdiff(tree$tip.label,names(x))
	if(length(xx)>0) stop("Some tips of the tree do not have data. Try model=\"BM\".")
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-8
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-FALSE
	likOU<-function(par,tree,x,trace){
		sig2<-par[1]
		alpha<-par[2]
		a0<-par[3]
		a<-par[1:(tree$Nnode-1)+3]
		logLik<-sum(dmnorm(c(x,a),mean=rep(a0,length(tree$tip.label)+tree$Nnode-1),
			varcov=sig2*vcvPhylo(tree,model="OU",alpha=alpha),log=TRUE))
		if(trace) print(c(sig2,alpha,logLik))
		-logLik
	}
	x<-x[tree$tip.label]
	pp<-rep(NA,tree$Nnode+2)
	pp[1:tree$Nnode+2]<-fastAnc(tree,x)
	pp[1]<-phyl.vcv(as.matrix(c(x,pp[2:tree$Nnode+2])),vcvPhylo(tree),lambda=1)$R[1,1]
	pp[2]<-tol ## arbitrarily
	fit<-optim(pp,likOU,tree=tree,x=x,trace=trace,method="L-BFGS-B",
		lower=c(tol,tol,rep(-Inf,tree$Nnode)),upper=rep(Inf,length(pp)),
		control=list(maxit=maxit))
	X<-list(sig2=fit$par[1],alpha=fit$par[2],
		ace=setNames(fit$par[1:tree$Nnode+2],1:tree$Nnode+length(tree$tip.label)),
		logLik=-fit$value,counts=fit$counts,convergence=fit$convergence,
		message=fit$message,model="OU")
	class(X)<-"anc.ML"
	X
}
