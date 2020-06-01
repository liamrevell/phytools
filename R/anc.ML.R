## lightweight version of ace(...,method="ML") for continuous traits
## also allows missing data in x, in which case missing data are also estimated
## written by Liam J. Revell 2011, 2013, 2014, 2015

anc.ML<-function(tree,x,maxit=2000,model=c("BM","OU","EB"),...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(model[1]=="BM") obj<-anc.BM(tree,x,maxit,...)
	else if(model[1]=="OU") obj<-anc.OU(tree,x,maxit,...)
	else if(model[1]=="EB") obj<-anc.EB(tree,x,maxit,...)
	else stop(paste("Do not recognize method",model))
	obj
}

## internal to estimate ancestral states under a BM model
## written by Liam J. Revell 2011, 2013, 2014, 2016, 2018
anc.BM<-function(tree,x,maxit,...){
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-FALSE
	if(hasArg(vars)) vars<-list(...)$vars
	else vars<-FALSE
	if(hasArg(CI)) CI<-list(...)$CI
	else CI<-FALSE
	if(hasArg(se)) se<-list(...)$se
	else se<-setNames(rep(0,length(x)),names(x))
	SE<-setNames(rep(0,Ntip(tree)),tree$tip.label)
	SE[names(se)]<-se
	E<-diag(SE)
	colnames(E)<-rownames(E)<-names(SE)
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-10*.Machine$double.eps
	## check to see if any tips are missing data
	xx<-setdiff(tree$tip.label,names(x))
	## function returns the log-likelihood
	likelihood<-function(par,C,invC,detC,xvals,msp,trace,E=0){
		sig2<-par[1]
		a<-par[2]
		y<-par[1:(tree$Nnode-1)+2]
		xvals<-c(xvals,setNames(par[1:length(msp)+tree$Nnode+1],msp))
		xvals<-xvals[rownames(C)[1:length(tree$tip.label)]]
		z<-c(xvals,y)-a
		if(trace) cat(paste(round(sig2,6)," --- ",sep=""))
		if(sum(E)>0){
			C<-sig2*C
			C[rownames(E),colnames(E)]<-C[rownames(E),colnames(E)]+E
			invC<-solve(C)
			detC<-determinant(C,logarithm=TRUE)$modulus[1]
		}
		logLik<-(-z%*%invC%*%z/(2*sig2)-nrow(C)*log(2*pi)/2-nrow(C)*log(sig2)/2-
			detC/2)[1,1]
		if(trace) cat(paste(round(logLik,6),"\n"))
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
		detC=detC,xvals=x,msp=xx,trace=trace,E=E,method="L-BFGS-B",
		lower=c(tol,rep(-Inf,tree$Nnode+length(xx))),
		control=list(maxit=maxit))
	if(vars||CI){ 
		H<-hessian(likelihood,fit$par,C=C,invC=invC,detC=detC,
			xvals=x,msp=xx,trace=trace,E=E)
		vcv<-solve(H)
	}
	states<-fit$par[1:tree$Nnode+1]
	names(states)<-c(length(tree$tip)+1,rownames(C)[(length(tree$tip)+1):nrow(C)])
	obj<-list(sig2=fit$par[1],ace=states,logLik=-fit$value,counts=fit$counts,
		convergence=fit$convergence,message=fit$message,model="BM")
	if(vars) obj$var<-setNames(diag(vcv)[1:tree$Nnode+1],
		c(length(tree$tip)+1,rownames(C)[(length(tree$tip)+1):nrow(C)]))
	if(CI){
		obj$CI95<-cbind(obj$ace-1.96*sqrt(diag(vcv)[1:tree$Nnode+1]),
			obj$ace+1.96*sqrt(diag(vcv)[1:tree$Nnode+1]))
		rownames(obj$CI95)<-c(length(tree$tip)+1,
			rownames(C)[(length(tree$tip)+1):nrow(C)])
	}
	if(length(xx)>0){
		obj$missing.x<-setNames(fit$par[1:length(xx)+tree$Nnode+1],xx)
		if(vars) obj$missing.var<-setNames(diag(vcv)[1:length(xx)+
			tree$Nnode+1],xx)
		if(CI){ 
			obj$missing.CI95<-cbind(obj$missing.x-
				1.96*sqrt(diag(vcv)[1:length(xx)+tree$Nnode+1]),
				obj$missing.x+1.96*sqrt(diag(vcv)[1:length(xx)+tree$Nnode+1]))
			rownames(obj$missing.CI95)<-xx
		}
	}
	class(obj)<-"anc.ML"
	obj
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
	if(hasArg(a.init)) a.init<-list(...)$a.init
	else a.init<-2*log(2)/max(nodeHeights(tree))
	likOU<-function(par,tree,x,trace){
		sig2<-par[1]
		alpha<-par[2]
		a0<-par[3]
		a<-par[1:(tree$Nnode-1)+3]
		logLik<-logMNORM(c(x,a),rep(a0,Ntip(tree)+tree$Nnode-1),sig2*vcvPhylo(tree,model="OU",alpha=alpha))
		if(trace) print(c(sig2,alpha,logLik))
		-logLik
	}
	x<-x[tree$tip.label]
	pp<-rep(NA,tree$Nnode+2)
	pp[1:tree$Nnode+2]<-fastAnc(tree,x)
	pp[1]<-phyl.vcv(as.matrix(c(x,pp[2:tree$Nnode+2])),vcvPhylo(tree),lambda=1)$R[1,1]
	pp[2]<-a.init ## arbitrarily
	fit<-optim(pp,likOU,tree=tree,x=x,trace=trace,method="L-BFGS-B",
		lower=c(tol,tol,rep(-Inf,tree$Nnode)),upper=rep(Inf,length(pp)),
		control=list(maxit=maxit))
	obj<-list(sig2=fit$par[1],alpha=fit$par[2],
		ace=setNames(fit$par[1:tree$Nnode+2],1:tree$Nnode+length(tree$tip.label)),
		logLik=-fit$value,counts=fit$counts,convergence=fit$convergence,
		message=fit$message,model="OU")
	class(obj)<-"anc.ML"
	obj
}

## EB is the Early-burst model (Harmon et al. 2010) and also called the ACDC model 
## (accelerating-decelerating; Blomberg et al. 2003). Set by the a rate parameter, EB fits a model where 
## the rate of evolution increases or decreases exponentially through time, under the model 
## r[t] = r[0] * exp(a * t), where r[0] is the initial rate, a is the rate change parameter, and t is 
## time. The maximum bound is set to -0.000001, representing a decelerating rate of evolution. The minimum 
## bound is set to log(10^-5)/depth of the tree.

## internal to estimate ancestral states under an EB model
## written by Liam J. Revell 2017
anc.EB<-function(tree,x,maxit=2000,...){
	## check to see if any tips are missing data
	xx<-setdiff(tree$tip.label,names(x))
	if(length(xx)>0) stop("Some tips of the tree do not have data. Try model=\"BM\".")
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-8
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-FALSE
	if(hasArg(vars)) vars<-list(...)$vars
	else vars<-FALSE
	if(hasArg(CI)) CI<-list(...)$CI
	else CI<-FALSE
	if(hasArg(r.init)){ 
		r.init<-list(...)$r.init
		obj<-phyl.vcv(as.matrix(x[tree$tip.label]),
			vcv(ebTree(tree,r.init)),1)
		s2.init<-obj$R[1,1]
		a0.init<-obj$alpha[1,1]
	} else { 
		## optimize r.init
		lik<-function(p,tree,x)
			logLik<--logMNORM(x,rep(p[3],Ntip(tree)),
				p[1]*vcvPhylo(tree,model="EB",r=p[2],anc.nodes=F))
		obj<-phyl.vcv(as.matrix(x[tree$tip.label]),vcv(tree),1)
		fit.init<-optim(c(obj$R[1,1],0,obj$alpha[1,1]),
			lik,tree=tree,x=x,method="L-BFGS-B",lower=c(tol,-Inf,-Inf),
			upper=rep(Inf,3))
		r.init<-fit.init$par[2]
		s2.init<-fit.init$par[1]
		a0.init<-fit.init$par[3]
	}
	likEB<-function(par,tree,x,trace){
		sig2<-par[1]
		r<-par[2]
		obj<-fastAnc(ebTree(tree,r),x)
		a0<-obj[1]
		a<-obj[2:length(obj)]
		logLik<-logMNORM(c(x,a),rep(a0,Ntip(tree)+tree$Nnode-1),
			sig2*vcvPhylo(tree,model="EB",r=r))
		if(trace) print(c(sig2,r,logLik))
		-logLik
	}
	x<-x[tree$tip.label]
	pp<-rep(NA,2)
	pp[1]<-s2.init
	pp[2]<-r.init
	fit<-optim(pp,likEB,tree=tree,x=x,trace=trace,method="L-BFGS-B",
		lower=c(tol,-Inf),upper=rep(Inf,2),control=list(maxit=maxit))
	obj<-list(sig2=fit$par[1],r=fit$par[2],
		ace=unclass(fastAnc(ebTree(tree,fit$par[2]),x)),
		logLik=-fit$value,counts=fit$counts,convergence=fit$convergence,
		message=fit$message,model="EB")
	if(vars||CI){
		likEB.hessian<-function(par,tree,y){
			sig2<-par[1]
			r<-par[2]
			a<-par[3:length(par)]
			logLik<-logMNORM(c(y,a[2:length(a)]),rep(a[1],Ntip(tree)+tree$Nnode-1),
				sig2*vcvPhylo(tree,model="EB",r=r))
			-logLik
		}
		H<-hessian(likEB.hessian,c(fit$par,obj$ace),tree=tree,y=x)
		vcv<-solve(H)
		if(vars) obj$var<-setNames(diag(vcv)[1:tree$Nnode+1],1:tree$Nnode+Ntip(tree))
		if(CI){
			obj$CI95<-cbind(obj$ace-1.96*sqrt(diag(vcv)[1:tree$Nnode+2]),
				obj$ace+1.96*sqrt(diag(vcv)[1:tree$Nnode+2]))
			rownames(obj$CI95)<-1:tree$Nnode+Ntip(tree)
		}
	}
	class(obj)<-"anc.ML"
	obj
}

logMNORM<-function(x,x0,vcv)
	-t(x-x0)%*%solve(vcv)%*%(x-x0)/2-length(x)*log(2*pi)/2-determinant(vcv,logarithm=TRUE)$modulus[1]/2

## print method for "anc.ML"
## written by Liam J. Revell 2015, 2016
print.anc.ML<-function(x,digits=6,printlen=NULL,...){
	cat(paste("Ancestral character estimates using anc.ML under a(n)",
		x$model,"model:\n"))
	Nnode<-length(x$ace)
	if(is.null(printlen)||printlen>=Nnode) print(round(x$ace,digits))
	else printDotDot(x$ace,digits,printlen)
	if(!is.null(x$var)){
		cat("\nVariances on ancestral states:\n")
		if(is.null(printlen)||printlen>=Nnode) print(round(x$var,digits))
		else printDotDot(x$var,digits,printlen)
	}
	if(!is.null(x$CI95)){
		cat("\nLower & upper 95% CIs:\n")
		colnames(x$CI95)<-c("lower","upper")
		if(is.null(printlen)||printlen>=Nnode) print(round(x$CI95,digits))
		else printDotDot(x$CI95,digits,printlen)
	}
	cat("\nFitted model parameters & likelihood:\n")
	if(x$model=="BM"){
		obj<-data.frame(round(x$sig2,digits),round(x$logLik,digits))
		colnames(obj)<-c("sig2","log-likelihood")
		rownames(obj)<-""
		print(obj)
	} else if(x$model=="OU"){
		obj<-data.frame(round(x$sig2,digits),round(x$alpha,digits),
			round(x$logLik,digits))
		colnames(obj)<-c("sigma^2","alpha","logLik")
		rownames(obj)<-""
		print(obj)
	} else if(x$model=="EB"){
		obj<-data.frame(round(x$sig2,digits),round(x$r,digits),
			round(x$logLik,digits))
		colnames(obj)<-c("sigma^2","r","logLik")
		rownames(obj)<-""
		print(obj)
	}
	if(x$convergence==0) cat("\nR thinks it has found the ML solution.\n\n") 
	else cat("\nOptimization may not have converged.\n\n")	
}

## S3 logLik method for "anc.ML" object class
logLik.anc.ML<-function(object,...){
	lik<-object$logLik
	if(object$model=="BM") attr(lik,"df")<-length(object$ace)+1
	else if(object$model=="EB") attr(lik,"df")<-length(object$ace)+2
	else if(object$model=="OU") attr(lik,"df")<-length(object$ace)+2
	lik
}
