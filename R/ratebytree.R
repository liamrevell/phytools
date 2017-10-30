## method to compare the rate of evolution for a character between trees
## closely related to 'censored' approach of O'Meara et al. (2006; Evolution)
## written by Liam J. Revell 2017

ratebytree<-function(trees,x,...){
	if(hasArg(type)) type<-list(...)$type
	else if(!missing(x)&&!is.null(x)) {
		if(is.factor(unlist(x))||is.character(unlist(x))) 
			type<-"discrete"
		else type<-"continuous"
	} else type<-"diversification"
	if(type=="continuous") obj<-rbt.cont(trees,x,...)
	else if(type=="discrete") obj<-rbt.disc(trees,x,...)
	else if(type=="diversification") obj<-rbt.div(trees,...)
	else {
		cat(paste("type =",type,"not recognized.\n"))
		obj<-NULL
	}
	obj
}

rbt.div<-function(trees,...){
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-FALSE
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	if(hasArg(test)) test<-list(...)$test
	else test<-"chisq"
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(hasArg(model)) model<-list(...)$model
	else model<-"birth-death"
	if(hasArg(rho)) rho<-list(...)$rho
	else rho<-rep(1,length(trees))
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-12
	if(!inherits(trees,"multiPhylo")) 
		stop("trees should be object of class \"multiPhylo\".")
	if(any(!sapply(trees,is.ultrametric))){
		cat("One or more trees fails check is.ultrametric.\n")
		cat("If you believe your tree to be ultrametric ")
		cat("use force.ultrametric.\n")
		stop()
	}
	t<-lapply(trees,function(phy) sort(branching.times(phy),
		decreasing=TRUE))
	if(model=="birth-death"){
		fit.multi<-mapply(fit.bd,tree=trees,rho=rho,SIMPLIFY=FALSE)
		logL.multi<-sum(sapply(fit.multi,logLik))
	} else if(model=="equal-extinction"){
		lik.eqmu<-function(theta,t,rho,trace=FALSE){
			lam<-theta[1:length(t)]
			mu<-theta[length(t)+1]
			logL<-0
			for(i in 1:length(t)) logL<-logL-lik.bd(c(lam[i],mu),t[[i]],rho[i])
			if(trace) cat(paste(paste(c(lam,mu,logL),collapse="\t"),"\n",sep=""))
			-logL
		}
		init.b<-sapply(trees,qb)
		obj<-nlminb(c(init.b,0),lik.eqmu,t=t,rho=rho,trace=trace,
			lower=rep(0,length(trees)+1),upper=rep(Inf,length(trees)+1))
		count<-0
		while(!is.finite(obj$objective)&&count<10){
			obj<-nlminb(runif(n=2,0,2)*c(init.b,0),lik.eqmu,t=t,rho=rho,
				trace=trace,lower=rep(0,length(trees)+1),
				upper=rep(Inf,length(trees)+1))
			count<-count+1
		}
		fit.multi<-vector(mode="list",length=length(t))
		for(i in 1:length(fit.multi)){ 
			fit.multi[[i]]$b<-obj$par[i]
			fit.multi[[i]]$d<-obj$par[length(obj$par)]
		}
		logL.multi<--obj$objective
	} else if(model=="equal-speciation"){
		lik.eqlam<-function(theta,t,rho,trace=FALSE){
			lam<-theta[1]
			mu<-theta[1:length(t)+1]
			logL<-0
			for(i in 1:length(t)) logL<-logL-lik.bd(c(lam,mu[i]),t[[i]],rho[i])
			if(trace) cat(paste(paste(c(lam,mu,logL),collapse="\t"),"\n",sep=""))
			-logL
		}
		init.b<-mean(sapply(trees,qb))
		obj<-nlminb(c(init.b,rep(0,length(trees))),lik.eqlam,t=t,rho=rho,
			trace=trace,lower=rep(0,length(trees)+1),
			upper=rep(Inf,length(trees)+1))
		count<-0
		while(!is.finite(obj$objective)&&count<10){
			obj<-nlminb(runif(n=2,0,2)*c(init.b,rep(0,length(trees))),lik.eqlam,
				t=t,rho=rho,trace=trace,lower=rep(0,length(trees)+1),
				upper=rep(Inf,length(trees)+1))
			count<-count+1
		}
		fit.multi<-vector(mode="list",length=length(t))
		for(i in 1:length(fit.multi)){ 
			fit.multi[[i]]$b<-obj$par[1]
			fit.multi[[i]]$d<-obj$par[i+1]
		}
		logL.multi<--obj$objective
	}
	lik.onerate<-function(theta,t,rho,model,trace=FALSE){
		logL<-sum(-mapply(lik.bd,t=t,rho=rho,
			MoreArgs=list(theta=theta)))
		if(trace) cat(paste(theta[1],"\t",theta[2],"\t",logL,"\n"))
		-logL
	}
	init.b<-mean(sapply(fit.multi,function(x) x$b))
	init.d<-mean(sapply(fit.multi,function(x) x$d))
	fit.onerate<-nlminb(c(init.b,init.d),lik.onerate,t=t,rho=rho,
		model=model,trace=trace,lower=c(0,0),upper=rep(Inf,2))
	count<-0
	while(!is.finite(fit.onerate$objective)&&count<10){
		fit.onerate<-nlminb(runif(n=2,0,2)*c(init.b,init.d),lik.onerate,
			t=t,rho=rho,model=model,trace=trace,lower=c(0,0),
			upper=rep(Inf,2))
		count<-count+1
	}
	rates.multi<-cbind(sapply(fit.multi,function(x) x$b),
			sapply(fit.multi,function(x) x$d))
	if(!is.null(names(trees))) rownames(trees)<-names(trees)
	else rownames(rates.multi)<-paste("tree",1:length(trees),sep="")
	colnames(rates.multi)<-c("b","d")
	LR<-2*(logL.multi+fit.onerate$objective)
	km<-if(model=="birth-death") 2*length(trees) else 
		if(model=="equal-extinction") length(trees)+1 else
		if(model=="equal-speciation") length(trees)+1 else
		if(model=="Yule") length(trees)
	k1<-if(model=="Yule") 1 else 2
	P.chisq<-pchisq(LR,df=km-k1,lower.tail=FALSE)
	obj<-list(
		multi.rate.model=list(
			logL=logL.multi,
			rates=rates.multi,
			k=km,
			method="nlminb"),
		common.rate.model=list(
			logL=-fit.onerate$objective,
			rates=setNames(fit.onerate$par,c("b","d")),
			k=k1,
			method="nlminb"),
		model=model,N=length(trees),
		n=sapply(trees,Ntip),
		likelihood.ratio=LR,P.chisq=P.chisq,
		type="diversification")
	class(obj)<-"ratebytree"
	obj
}

## used (for now) to get a starting value for optimization in type="diversification"
qb<-function(tree) (log(Ntip(tree))-log(2))/max(nodeHeights(tree))

## discrete character ratebytree
rbt.disc<-function(trees,x,...){
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-FALSE
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	if(hasArg(test)) test<-list(...)$test
	else test<-"chisq"
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(hasArg(model)) model<-list(...)$model
	else model<-"ER"
	if(!inherits(trees,"multiPhylo")) 
		stop("trees should be object of class \"multiPhylo\".")
	if(!is.list(x)) stop("x should be a list of vectors.")
	N<-length(trees)
	fit.multi<-mapply(fitMk,tree=trees,x=x,MoreArgs=list(model=model),
		SIMPLIFY=FALSE)
	logL.multi<-sum(sapply(fit.multi,logLik))
	lik.onerate<-function(theta,trees,x,model,trace=FALSE){
		ss<-sort(unique(unlist(x)))
		m<-length(ss)
		if(is.character(model)){
			rate<-matrix(NA,m,m)
			if(model=="ER"){
				k<-rate[]<-1
				diag(rate)<-NA
			}
			else if(model=="ARD") {
				k<-m*(m-1)
				rate[col(rate)!=row(rate)]<-1:k
			}
			else if(model=="SYM") {
				k<-m*(m-1)/2
				ii<-col(rate)<row(rate)
				rate[ii]<-1:k
				rate<-t(rate)
				rate[ii]<-1:k
			}
		} else {
			rate<-model
			k<-max(rate)
		}
		Q <- matrix(0, m, m)
		index.matrix <- rate
		tmp <- cbind(1:m, 1:m)
		rate[tmp] <- 0
		rate[rate == 0] <- k + 1
		Q[]<-c(theta,0)[rate]
		diag(Q)<--rowSums(Q)
		colnames(Q)<-rownames(Q)<-ss
		fit.onerate<-mapply(fitMk,tree=trees,x=x,MoreArgs=list(fixedQ=Q),
			SIMPLIFY=FALSE)
		-sum(sapply(fit.onerate,logLik))
	}
	pp<-sapply(fit.multi,function(x) x$rates)
	pp<-if(is.matrix(pp)) rowMeans(pp) else mean(pp)
	if(length(pp)>1){
		fit.onerate<-optim(pp,lik.onerate,trees=trees,x=x,model=model)
	} else {
		fit.onerate<-optimize(lik.onerate,c(0,1000*pp),trees=trees,x=x,
			model=model)
		names(fit.onerate)<-c("par","value")
	}
	rates.multi<-t(sapply(fit.multi,function(x) x$rates))
	if(is.character(model)){
		m<-length(fit.multi[[1]]$states)
		if(model=="ER"){
			rates.multi<-t(rates.multi)
			colnames(rates.multi)<-"q"
		} else if(model=="SYM"){
			if(length(fit.multi[[1]]$states)==2) rates.multi<-t(rates.multi)
			colnames(rates.multi)<-
				sapply(fit.multi[[1]]$states,function(x,y) 
					sapply(y,function(y,x) paste(x,"<->",y,sep=""),x=x),
					y=fit.multi[[1]]$states)[lower.tri(matrix(0,m,m))]
		} else if(model=="ARD"){
			ii<-(upper.tri(matrix(0,m,m))+lower.tri(matrix(0,m,m))==TRUE)
			colnames(rates.multi)<-
				sapply(fit.multi[[1]]$states,function(x,y) 
					sapply(y,function(y,x) paste(y,"->",x,sep=""),x=x),
					y=fit.multi[[1]]$states)[ii]
		}			
	} else {
		k<-max(fit.multi[[1]]$index.matrix,na.rm=TRUE)
		if(k==1) rates.multi<-t(rates.multi)
		foo<-function(i,index.matrix,states){
			rc<-which(index.matrix==i,arr.ind=TRUE)
			lab<-apply(rc,1,function(ind,ss) paste(ss[ind],collapse="->"),
				ss=states)
			if(length(lab)>1) paste(lab,collapse=",") else lab
		}
		colnames(rates.multi)<-sapply(1:k,
			foo,fit.multi[[1]]$index.matrix,fit.multi[[1]]$states)
	}
	if(!is.null(names(trees))) rownames(rates.multi)<-names(trees)
	else rownames(rates.multi)<-paste("tree",1:length(trees),sep="")
	LR<-2*(logL.multi+fit.onerate$value)
	km<-sum(sapply(fit.multi,function(x) max(x$index.matrix,na.rm=TRUE)))
	k1<-length(pp)
	P.chisq<-pchisq(LR,df=km-k1,lower.tail=FALSE)
	obj<-list(
		multi.rate.model=list(
			logL=logL.multi,
			rates=rates.multi,
			method=fit.multi[[1]]$method),
		common.rate.model=list(
			logL=-fit.onerate$value,
			rates=setNames(fit.onerate$par,colnames(rates.multi)),
			method=if(length(pp)>1) "optim" else "optimize"),
		index.matrix=fit.multi[[1]]$index.matrix,
		states=fit.multi[[1]]$states,
		pi=fit.multi[[1]]$pi,
		model=model,N=N,n=sapply(trees,Ntip),
		likelihood.ratio=LR,P.chisq=P.chisq,
		type="discrete")
	class(obj)<-"ratebytree"
	obj
}

## continuous character ratebytree
rbt.cont<-function(trees,x,...){
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-8
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-FALSE
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	if(hasArg(test)) test<-list(...)$test
	else test<-"chisq"
	if(hasArg(quiet)) quiet<-list(...)$quiet	
	else quiet<-FALSE
	if(hasArg(maxit)) maxit<-list(...)$maxit
	else maxit<-500
	if(hasArg(regimes)){ 
		regimes<-list(...)$regimes
		if(!is.factor(regimes)) regimes<-as.factor(regimes)
	} else regimes<-as.factor(1:length(trees))
	if(hasArg(model)) model<-list(...)$model
	else model<-"BM"
	if(!(model%in%c("BM","OU","EB"))){
		cat(paste("model =",model,"not recognized. using model = \"BM\"\n"))
		model<-"BM"
	}
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
	m<-length(levels(regimes))
	## reorder the trait vectors in x & SEs in se
	x<-mapply(function(x,t) x<-x[t$tip.label],x=x,t=trees,SIMPLIFY=FALSE)
	se<-mapply(function(x,t) x<-x[t$tip.label],x=se,t=trees,SIMPLIFY=FALSE)
	## likelihood functions
	lik.multi<-function(theta,trees,y,se,model,regimes,trace=FALSE){
		m<-length(unique(regimes))
		ind<-setNames(lapply(levels(regimes),function(x,y) which(y==x),y=regimes),
			levels(regimes))
		THETA<-list()
		for(i in 1:m){
			THETA[[i]]<-c(theta[i],theta[ind[[i]]+m])
			if(model%in%c("OU","EB")) THETA[[i]]<-c(THETA[[i]],theta[i+m+N])
		}
		logL<-0
		for(i in 1:m)
			logL<-logL-lik.onerate(THETA[[i]],trees[ind[[i]]],y[ind[[i]]],se[ind[[i]]],
				model=model,trace=FALSE)
		if(trace){
			cat(paste(paste(round(theta[1:m],digits),collapse="\t"),
				if(model=="OU") paste(round(theta[1:m+m+N],digits),collapse="\t"),
				if(model=="EB") paste(round(theta[1:m+m+N],digits),collapse="\t"),
				round(logL,digits),"\n",sep="\t"))
			flush.console()
		}
		-logL
	}
	lik.onerate<-function(theta,trees,y,se,model,trace=FALSE){
		n<-sapply(trees,Ntip)
		N<-length(trees)
		sig<-theta[1]
		a<-theta[1:N+1]
		if(model=="OU") alpha<-theta[N+2]
		if(model=="EB") r<-theta[N+2]
		if(model=="BM") C<-lapply(trees,vcv)
		else if(model=="OU") C<-lapply(trees,vcvPhylo,model="OU",alpha=alpha,
			anc.nodes=FALSE)
		else if(model=="EB") C<-lapply(trees,vcvPhylo,model="EB",r=r,
			anc.nodes=FALSE)
		E<-lapply(se,diag)
		V<-mapply("+",lapply(C,"*",sig),E,SIMPLIFY=FALSE)
		logL<-0
		for(i in 1:N) 
			logL<-logL-t(y[[i]]-a[i])%*%solve(V[[i]])%*%(y[[i]]-
				a[i])/2-n[i]*log(2*pi)/2-determinant(V[[i]])$modulus[1]/2
		if(trace){ 
			cat(paste(round(sig,digits),
				if(model=="OU") round(alpha,digits),
				if(model=="EB") round(r,digits),
				round(logL,digits),"\n",sep="\t"))
			flush.console()
		}
		-logL
	}
	## first, fit multi-rate model
	f1<-function(tree,x){ 
		pvcv<-phyl.vcv(as.matrix(x),vcv(tree),1)
		c(pvcv$R[1,1],pvcv$a[1,1])
	}
	P<-mapply(f1,trees,x,SIMPLIFY=FALSE)
	PP<-list()
	PP[[1]]<-vector()
	for(i in 1:m)
		PP[[1]][i]<-mean(sapply(P,function(x) x[1])[which(regimes==levels(regimes)[i])])
	PP[[2]]<-sapply(P,function(x) x[2])
	if(model%in%c("OU","EB")) PP[[3]]<-rep(0,m)
	if(hasArg(init)){ 
		init<-list(...)$init
	 	if(!is.null(init$sigm)) PP[[1]]<-init$sigm
		if(!is.null(init$am)) PP[[2]]<-init$am
		if(model=="OU") if(!is.null(init$alpham)) PP[[3]]<-init$alpham
		if(model=="EB") if(!is.null(init$rm)) PP[[3]]<-init$rm
	}
	p<-unlist(PP)
	if(trace){
		if(model=="BM"){
			cat("\nOptimizing multi-rate model....\n")
			cat(paste(paste("sig[",1:m,"]",sep="",collapse="\t"),"logL\n",sep="\t"))
		} else if(model=="OU"){
			cat("\nOptimizing multi-regime model....\n")
			cat(paste(paste("sig[",1:m,"]",sep="",collapse="\t"),
				paste("alpha[",1:m,"]",sep="",collapse="\t"),"logL\n",sep="\t"))
		} else if(model=="EB"){
			cat("\nOptimizing multi-regime model....\n")
			cat(paste(paste("sig[",1:m,"]",sep="",collapse="\t"),
				paste("r[",1:m,"]",sep="",collapse="\t"),"logL\n",sep="\t"))
		}
	}
	fit.multi<-optim(p,lik.multi,trees=trees,y=x,se=se,model=model,
		regimes=regimes,trace=trace,
		method="L-BFGS-B",lower=c(rep(tol,m),rep(-Inf,N),
			if(model%in%c("OU","EB")) rep(-Inf,m)),upper=c(rep(Inf,m),rep(Inf,N),
			if(model%in%c("OU","EB")) rep(Inf,m)),control=list(maxit=maxit))
	## compute covariance matrix
	H.multi<-hessian(lik.multi,fit.multi$par,trees=trees,y=x,
		se=se,model=model,regimes=regimes,trace=FALSE)
	Cov.multi<-if(qr(H.multi)$rank!=ncol(H.multi)) ginv(H.multi) else solve(H.multi)
	## now fit single-rate model
	p<-c(mean(fit.multi$par[1:m]),fit.multi$par[1:N+m])
	if(model%in%c("OU","EB")) p<-c(p,mean(fit.multi$par[1:m+m+N]))
	if(hasArg(init)){
		if(!is.null(init$sigc)) p[1]<-init$sigc
		if(!is.null(init$ac)) p[1:N+1]<-init$ac
		if(model=="OU") if(!is.null(init$alphac)) p[N+2]<-init$alphac
		if(model=="EB") if(!is.null(init$rc)) p[N+2]<-init$rc
	}
	if(trace){
		if(model=="BM"){
			cat("\nOptimizing common-rate model....\n")
			cat(paste("sig  ","logL\n",sep="\t"))
		} else if(model=="OU"){
			cat("\nOptimizing common-regime model....\n")
			cat(paste("sig  ","alpha ","logL\n",sep="\t"))
		} else if(model=="EB"){
			cat("\nOptimizing common-regime mode.....\n")
			cat(paste("sig  ","r    ","logL\n",sep="\t"))
		}
	}
	fit.onerate<-optim(p,lik.onerate,trees=trees,y=x,se=se,model=model,
		trace=trace,method="L-BFGS-B",
		lower=c(tol,rep(-Inf,N),if(model=="OU") tol else if(model=="EB") -Inf),
		upper=c(Inf,rep(Inf,N),if(model%in%c("OU","EB")) Inf),
		control=list(maxit=maxit))
	## compute covariance matrix
	H.onerate<-hessian(lik.onerate,fit.onerate$par,trees=trees,y=x,
		se=se,model=model,trace=FALSE)
	Cov.onerate<-if(qr(H.onerate)$rank!=ncol(H.onerate)) ginv(H.onerate) else solve(H.onerate)
	## compare models:
	LR<-2*(-fit.multi$value+fit.onerate$value)
	km<-N+m+if(model=="BM") 0 else if(model%in%c("OU","EB")) m
	k1<-N+if(model=="BM") 1 else if(model%in%c("OU","EB")) 2
	if(test=="simulation"&&model%in%c("OU","EB")){
		cat("Simulation test not yet available for chosen model. Using chi-square test.\n")
		test<-"chisq"
	}
	if(test=="chisq") P.chisq<-pchisq(LR,df=km-k1,lower.tail=FALSE)
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
			f2<-function(x,se) sampleFrom(xbar=x,xvar=se^2,n=rep(1,length(x)))
			x.sim<-mapply(f2,x=x.sim,se=se,SIMPLIFY=FALSE)
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
		multi.rate.model=list(sig2=fit.multi$par[1:m],
			SE.sig2=sqrt(diag(Cov.multi)[1:m]),
			a=fit.multi$par[1:N+m],
			SE.a=sqrt(diag(Cov.multi)[1:N+m]),
			alpha=if(model=="OU") fit.multi$par[1:m+m+N] else NULL,
			SE.alpha=if(model=="OU") sqrt(diag(Cov.multi)[1:m+m+N]) else NULL,
			r=if(model=="EB") fit.multi$par[1:N+m+N] else NULL,
			SE.r=if(model=="EB") sqrt(diag(Cov.multi)[1:N+m+N]) else NULL,
			k=km,
			logL=-fit.multi$value,
			counts=fit.multi$counts,convergence=fit.multi$convergence,
			message=fit.multi$message),
		common.rate.model=list(sig2=fit.onerate$par[1],
			SE.sig2=sqrt(diag(Cov.onerate)[1]),
			a=fit.onerate$par[1:N+1],
			SE.a=sqrt(diag(Cov.onerate)[1:N+1]),
			alpha=if(model=="OU") fit.onerate$par[N+2] else NULL,
			SE.alpha=if(model=="OU") sqrt(diag(Cov.onerate)[N+2]) else NULL,
			r=if(model=="EB") fit.onerate$par[N+2] else NULL,
			SE.r=if(model=="EB") sqrt(diag(Cov.onerate)[N+2]) else NULL,
			k=k1,
			logL=-fit.onerate$value,
			counts=fit.onerate$counts,convergence=fit.onerate$convergence,
			message=fit.onerate$message),
		model=model,
		N=N,n=sapply(trees,Ntip),likelihood.ratio=LR,
		regimes=regimes,m=m,
		P.chisq=if(test=="chisq") P.chisq else NULL,
		P.sim=if(test=="simulation") P.sim else NULL,
		type="continuous")
	class(obj)<-"ratebytree"
	obj
}

## S3 print method for ratebytree
print.ratebytree<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	if(x$type=="continuous") prbt.cont(x,digits=digits)
	else if(x$type=="discrete") prbt.disc(x,digits=digits)
	else if(x$type=="diversification") prbt.div(x,digits=digits)
	else print(x)
}

prbt.div<-function(x,digits=digits){
	cat("ML common diversification-rate model:")
	cat("\n\tb\td\tk\tlog(L)")
	cat(paste("\nvalue",round(x$common.rate.model$rates[1],digits),
		round(x$common.rate.model$rates[2],digits),
		x$common.rate.model$k,round(x$common.rate.model$logL,digits),
		sep="\t"))
	cat("\n\nML multi diversification-rate model:")
	cat(paste("\n",paste(paste("b[",1:x$N,"]",sep=""),collapse="\t"),
		paste(paste("d[",1:x$N,"]",sep=""),collapse="\t"),"k\tlog(L)",
		sep="\t"))
	cat(paste("\nvalue",paste(round(x$multi.rate.model$rates[,"b"],digits),
		collapse="\t"),paste(round(x$multi.rate.model$rates[,"d"],digits),
		collapse="\t"),x$multi.rate.model$k,round(x$multi.rate.model$logL,
		digits),sep="\t"))
	cat(paste("\n\nDiversification model was \"",x$model,"\".\n",sep=""))
	cat(paste("Model fitting method was \"",x$multi.rate.model$method,
		"\".\n",sep=""))
	cat(paste("\nLikelihood ratio:",round(x$likelihood.ratio,digits),"\n"))
	cat(paste("P-value (based on X^2):",round(x$P.chisq,digits),"\n\n"))
}

prbt.cont<-function(x,digits=digits){
	N<-x$N
	m<-x$m
	if(x$model=="BM"){
		cat("ML common-rate model:\n")
		cat(paste("\ts^2\t",paste(paste("a[",1:N,"]",sep=""),collapse="\t")),
			"\tk\tlogL\n")
		cat(paste("value",round(x$common.rate.model$sig2,digits),
			paste(round(x$common.rate.model$a,digits),collapse="\t"),
			x$common.rate.model$k,round(x$common.rate.model$logL,digits),
			"\n",sep="\t"))
		cat(paste("SE   ",round(x$common.rate.model$SE.sig2,digits),
			paste(round(x$common.rate.model$SE.a,digits),collapse="\t"),
			"\n\n",sep="\t"))
		cat("ML multi-rate model:\n")
		cat(paste("\t",paste(paste("s^2[",levels(x$regimes),"]",sep=""),collapse="\t"),"\t",
			paste(paste("a[",1:N,"]",sep=""),collapse="\t")),
			"\tk\tlogL\n")
		cat(paste("value",paste(round(x$multi.rate.model$sig2,digits),collapse="\t"),
			paste(round(x$multi.rate.model$a,digits),collapse="\t"),
			x$multi.rate.model$k,round(x$multi.rate.model$logL,digits),
			"\n",sep="\t"))
		cat(paste("SE   ",paste(round(x$multi.rate.model$SE.sig2,digits),collapse="\t"),
			paste(round(x$multi.rate.model$SE.a,digits),collapse="\t"),
			"\n\n",sep="\t"))
	} else if(x$model=="OU"){
		cat("ML common-regime OU model:\n")
		cat(paste("\ts^2\t",paste(paste("a[",1:N,"]",sep=""),collapse="\t")),
			"\talpha\tk\tlogL\n")
		cat(paste("value",round(x$common.rate.model$sig2,digits),
			paste(round(x$common.rate.model$a,digits),collapse="\t"),
			round(x$common.rate.model$alpha,digits),
			x$common.rate.model$k,round(x$common.rate.model$logL,digits),
			"\n",sep="\t"))
		cat(paste("SE   ",round(x$common.rate.model$SE.sig2,digits),
			paste(round(x$common.rate.model$SE.a,digits),collapse="\t"),
			round(x$common.rate.model$SE.alpha,digits),
			"\n\n",sep="\t"))
		cat("ML multi-regime OU model:\n")
		cat(paste("\t",paste(paste("s^2[",levels(x$regimes),"]",sep=""),collapse="\t"),"\t",
			paste(paste("a[",1:N,"]",sep=""),collapse="\t"),"\t",
			paste(paste("alp[",levels(x$regimes),"]",sep=""),collapse="\t")),
			"\tk\tlogL\n")
		cat(paste("value",paste(round(x$multi.rate.model$sig2,digits),collapse="\t"),
			paste(round(x$multi.rate.model$a,digits),collapse="\t"),
			paste(round(x$multi.rate.model$alpha,digits),collapse="\t"),
			x$multi.rate.model$k,round(x$multi.rate.model$logL,digits),
			"\n",sep="\t"))
		cat(paste("SE   ",paste(round(x$multi.rate.model$SE.sig2,digits),collapse="\t"),
			paste(round(x$multi.rate.model$SE.a,digits),collapse="\t"),
			paste(round(x$multi.rate.model$SE.alpha,digits),collapse="\t"),
			"\n\n",sep="\t"))
	} else if(x$model=="EB"){
		cat("ML common-regime EB model:\n")
		cat(paste("\ts^2\t",paste(paste("a[",1:N,"]",sep=""),collapse="\t")),
			"\tr\tk\tlogL\n")
		cat(paste("value",round(x$common.rate.model$sig2,digits),
			paste(round(x$common.rate.model$a,digits),collapse="\t"),
			round(x$common.rate.model$r,digits),
			x$common.rate.model$k,round(x$common.rate.model$logL,digits),
			"\n",sep="\t"))
		cat(paste("SE   ",round(x$common.rate.model$SE.sig2,digits),
			paste(round(x$common.rate.model$SE.a,digits),collapse="\t"),
			round(x$common.rate.model$SE.r,digits),
			"\n\n",sep="\t"))
		cat("ML multi-regime EB model:\n")
		cat(paste("\t",paste(paste("s^2[",levels(x$regimes),"]",sep=""),collapse="\t"),"\t",
			paste(paste("a[",1:N,"]",sep=""),collapse="\t"),"\t",
			paste(paste("r[",levels(x$regimes),"]",sep=""),collapse="\t")),
			"\tk\tlogL\n")
		cat(paste("value",paste(round(x$multi.rate.model$sig2,digits),collapse="\t"),
			paste(round(x$multi.rate.model$a,digits),collapse="\t"),
			paste(round(x$multi.rate.model$r,digits),collapse="\t"),
			x$multi.rate.model$k,round(x$multi.rate.model$logL,digits),
			"\n",sep="\t"))
		cat(paste("SE   ",paste(round(x$multi.rate.model$SE.sig2,digits),collapse="\t"),
			paste(round(x$multi.rate.model$SE.a,digits),collapse="\t"),
			paste(round(x$multi.rate.model$SE.r,digits),collapse="\t"),
			"\n\n",sep="\t"))
	}
	cat(paste("Likelihood ratio:",round(x$likelihood.ratio,digits),"\n"))
	if(!is.null(x$P.chisq))
		cat(paste("P-value (based on X^2):",round(x$P.chisq,digits),"\n\n"))
	else if(!is.null(x$P.sim))
		cat(paste("P-value (based on simulation):",round(x$P.sim,digits),"\n\n"))
	if(x$multi.rate.model$convergence==0&&x$common.rate.model$convergence==0) 
        cat("R thinks it has found the ML solution.\n\n")
	else cat("One or the other optimization may not have converged.\n\n")
}

prbt.disc<-function(x,digits=digits){
	cat("ML common-rate model:\n")
	
	cat(paste("\t",paste(names(x$common.rate.model$rates),collapse="\t"),
		"\tk\tlogL\n",sep=""))
	cat(paste("value\t",paste(round(x$common.rate.model$rates,digits),
		collapse="\t"),"\t",max(x$index.matrix,na.rm=T),"\t",
		round(x$common.rate.model$logL,digits),"\n",sep=""))
	cat(paste("\nModel fitting method was \"",x$common.rate.model$method,
		"\".\n",sep=""))
	cat("\nML multi-rate model:\n")
	cat(paste("\t",paste(colnames(x$multi.rate.model$rates),collapse="\t"),
		"\tk\tlogL\n",sep=""))
	for(i in 1:nrow(x$multi.rate.model$rates)){
		if(i>1) cat(paste(rownames(x$multi.rate.model$rates)[i],"\t",
			paste(round(x$multi.rate.model$rates[i,],
			digits),collapse="\t"),"\n",sep=""))
		else if(i==1) cat(paste(rownames(x$multi.rate.model$rates)[i],"\t",
			paste(round(x$multi.rate.model$rates[i,],
			digits),collapse="\t"),"\t",x$N*max(x$index.matrix,na.rm=T),"\t",
			round(x$multi.rate.model$logL,digits),"\n",sep=""))
	}
	cat(paste("\nModel fitting method was \"",x$multi.rate.model$method,
		"\".\n",sep=""))
	cat(paste("\nLikelihood ratio:",round(x$likelihood.ratio,digits),"\n"))
	cat(paste("P-value (based on X^2):",round(x$P.chisq,digits),"\n\n"))
}

## posthoc comparison S3 method

posthoc<-function(x, ...) UseMethod("posthoc")

posthoc.ratebytree<-function(x,...){
	if(hasArg(p.adjust.method)) p.adjust.method<-list(...)$p.adjust.method
	else p.adjust.method<-"none"
	if(x$type!="continuous"){
		cat("Sorry. No posthoc method yet implemented for this data type.\n\n")
	} else {
		if(x$model=="BM") k<-2
		else if(x$model%in%c("EB","OU")) k<-3
		t<-df<-P<-matrix(0,x$m,x$m)
		for(i in 1:x$m){
			for(j in 1:x$m){
				x1<-if(x$model=="BM") x$multi.rate.model$sig2[i]
					else if(x$model=="OU") x$multi.rate.model$alpha[i]
					else if(x$model=="EB") x$multi.rate.model$r[i]
				x2<-if(x$model=="BM") x$multi.rate.model$sig2[j]
					else if(x$model=="OU") x$multi.rate.model$alpha[j]
					else if(x$model=="EB") x$multi.rate.model$r[j]
				s1<-x$multi.rate.model$SE.sig2[i]^2
				s2<-x$multi.rate.model$SE.sig2[j]^2
				n1<-x$n[i]
				n2<-x$n[j]
				se<-sqrt(s1+s2)
				df[i,j]<-(s1/n1+s2/n2)^2/
					((s1/n1)^2/(n1-k)+(s2/n2)^2/(n2-k))
				t[i,j]<-(x1-x2)/se
				P[i,j]<-2*pt(abs(t[i,j]),df[i,j],lower.tail=FALSE)		
			}
		}
		p<-P[upper.tri(P)]
		p<-p.adjust(p,method=p.adjust.method)
		P[upper.tri(P)]<-p
		P[lower.tri(P)]<-p		
	}
	obj<-list(t=t,df=df,P=P,model=x$model,type=x$type,
		p.adjust.method=p.adjust.method)
	class(obj)<-"posthoc.ratebytree"
	obj
}

print.posthoc.ratebytree<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	t<-x$t[upper.tri(x$t)]
	df<-x$df[upper.tri(x$df)]
	P<-x$P[upper.tri(x$P)]
	N<-nrow(x$P)
	paste("reg.",1:(N-1),"vs.",2:N)
	nn<-vector(mode="character",length=N*(N-1)/2)
	k<-1
	for(i in 1:N) for(j in i:N){
		if(i!=j){
			nn[k]<-paste("reg.",i," vs. ",j,sep="")
			k<-k+1
		}
	}
	X<-data.frame(t=t,df=df,P=P)
	rownames(X)<-nn
	cat(paste("\nPost-hoc test for \"",x$model,"\" model.\n",sep=""))
	cat(paste("(Comparison is of estimated values of",
		if(x$model=="BM") "sigma^2.)\n\n"
		else if(x$model=="OU") "alpha.)\n\n"
		else if(x$model=="EB") "r.)\n\n"))
	print(round(X,digits))
	cat(paste("\nP-values adjusted using method=\"",x$p.adjust.method,
		"\".\n\n",sep=""))
}

AIC.ratebytree<-function(object,...,k=2){
	aic<-data.frame(AIC=c(k*object$common.rate.model$k-2*object$common.rate.model$logL,
		k*object$multi.rate.model$k-2*object$multi.rate.model$logL),
		df=c(object$common.rate.model$k,object$multi.rate.model$k))
	model.names<-if(is.null(object$model)) 1 else object$model
	addtl.obj<-list(...)
	if(length(addtl.obj)>0){
		for(i in 1:length(addtl.obj)){ 
			aic<-rbind(aic,c(k*addtl.obj[[i]]$multi.rate.model$k-
				2*addtl.obj[[i]]$multi.rate.model$logL,
				addtl.obj[[i]]$multi.rate.model$k))
			model.names<-c(model.names,
				if(is.null(addtl.obj[[i]]$model)) 1+i else addtl.obj[[i]]$model)
		}
		rownames(aic)<-c("common-rate",paste("multi-rate:",model.names,sep=""))
	} else rownames(aic)<-c("common-rate","multi-rate")
	aic
}
