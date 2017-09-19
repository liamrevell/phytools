## method to compare the rate of evolution for a character between trees
## closely related to 'censored' approach of O'Meara et al. (2006; Evolution)
## written by Liam J. Revell 2017

ratebytree<-function(trees,x,...){
	if(hasArg(type)) type<-list(...)$type
	else {
		if(is.factor(unlist(x))||is.character(unlist(x))) 
			type<-"discrete"
		else type<-"continuous"
	}
	if(type=="continuous") obj<-rbt.cont(trees,x,...)
	else if(type=="discrete") obj<-rbt.disc(trees,x,...)
	else {
		cat(paste("type =",type,"not recognized.\n"))
		obj<-NULL
	}
	obj
}

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
	if(!is.null(names(trees))) rownames(trees)<-names(trees)
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
	## reorder the trait vectors in x & SEs in se
	x<-mapply(function(x,t) x<-x[t$tip.label],x=x,t=trees,SIMPLIFY=FALSE)
	se<-mapply(function(x,t) x<-x[t$tip.label],x=se,t=trees,SIMPLIFY=FALSE)
	## first, fit multi-rate model
	lik.multi<-function(theta,trees,y,se,model,trace=FALSE){
		n<-sapply(trees,Ntip)
		N<-length(trees)
		sig<-theta[1:N]
		a<-theta[1:N+N]
		if(model=="OU") alpha<-theta[1:N+2*N]
		if(model=="EB") r<-theta[1:N+2*N]
		if(model=="BM") C<-lapply(trees,vcv)
		else if(model=="OU") C<-mapply(vcvPhylo,tree=trees,alpha=alpha,MoreArgs=list(model="OU",
			anc.nodes=FALSE),SIMPLIFY=FALSE)
		else if(model=="EB") C<-mapply(vcvPhylo,tree=trees,r=r,MoreArgs=list(model="EB",
			anc.nodes=FALSE),SIMPLIFY=FALSE)
		E<-lapply(se,diag)
		V<-mapply("+",mapply("*",C,sig,SIMPLIFY=FALSE),E,SIMPLIFY=FALSE)
		logL<-0
		for(i in 1:N) 
			logL<-logL-t(y[[i]]-a[i])%*%solve(V[[i]])%*%(y[[i]]-
				a[i])/2-n[i]*log(2*pi)/2-determinant(V[[i]])$modulus[1]/2
		if(trace){
			cat(paste(paste(round(sig,digits),collapse="\t"),
				if(model=="OU") paste(round(alpha,digits),collapse="\t"),
				if(model=="EB") paste(round(r,digits),collapse="\t"),
				round(logL,digits),"\n",sep="\t"))
			flush.console()
		}
		-logL
	}	
	f1<-function(tree,x){ 
		pvcv<-phyl.vcv(as.matrix(x),vcv(tree),1)
		c(pvcv$R[1,1],pvcv$a[1,1])
	}
	PP<-mapply(f1,trees,x)
	if(hasArg(init)){ 
		init<-list(...)$init
	 	if(!is.null(init$sigm)) PP[1,]<-init$sigm
		if(!is.null(init$am)) PP[2,]<-init$am
		if(model=="OU") if(!is.null(init$alpham)) PP<-rbind(PP,init$alpham)
		if(model=="EB") if(!is.null(init$rm)) PP<-rbind(PP,init$rm)
	}
	if((model%in%c("OU","EB"))&&nrow(PP)==2) PP<-rbind(PP,rep(0,ncol(PP)))
	p<-as.vector(t(PP))
	if(trace){
		if(model=="BM"){
			cat("\nOptimizing multi-rate model....\n")
			cat(paste(paste("sig[",1:N,"]",sep="",collapse="\t"),"logL\n",sep="\t"))
		} else if(model=="OU"){
			cat("\nOptimizing multi-regime model....\n")
			cat(paste(paste("sig[",1:N,"]",sep="",collapse="\t"),
				paste("alpha[",1:N,"]",sep="",collapse="\t"),"logL\n",sep="\t"))
		} else if(model=="EB"){
			cat("\nOptimizing multi-regime model....\n")
			cat(paste(paste("sig[",1:N,"]",sep="",collapse="\t"),
				paste("r[",1:N,"]",sep="",collapse="\t"),"logL\n",sep="\t"))
		}
	}
	fit.multi<-optim(p,lik.multi,trees=trees,y=x,se=se,model=model,trace=trace,
		method="L-BFGS-B",lower=c(rep(tol,N),rep(-Inf,N),
			if(model%in%c("OU","EB")) rep(-Inf,N)),upper=c(rep(Inf,N),rep(Inf,N),
			if(model%in%c("OU","EB")) rep(Inf,N)))
	## compute covariance matrix
	H.multi<-hessian(lik.multi,fit.multi$par,trees=trees,y=x,
		se=se,model=model,trace=FALSE)
	Cov.multi<-if(qr(H.multi)$rank!=ncol(H.multi)) ginv(H.multi) else solve(H.multi)
	## now fit single-rate model
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
	p<-c(mean(fit.multi$par[1:N]),fit.multi$par[1:N+N])
	if(model%in%c("OU","EB")) p<-c(p,mean(fit.multi$par[1:N+2*N]))
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
		upper=c(Inf,rep(Inf,N),if(model%in%c("OU","EB")) Inf))
		## compute covariance matrix
	H.onerate<-hessian(lik.onerate,fit.onerate$par,trees=trees,y=x,
		se=se,model=model,trace=FALSE)
	Cov.onerate<-if(qr(H.onerate)$rank!=ncol(H.onerate)) ginv(H.onerate) else solve(H.onerate)
	## compare models:
	LR<-2*(-fit.multi$value+fit.onerate$value)
	km<-2*N+if(model=="BM") 0 else if(model%in%c("OU","EB")) N
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
		multi.rate.model=list(sig2=fit.multi$par[1:N],
			SE.sig2=sqrt(diag(Cov.multi)[1:N]),
			a=fit.multi$par[1:N+N],
			SE.a=sqrt(diag(Cov.multi)[1:N+N]),
			alpha=if(model=="OU") fit.multi$par[1:N+2*N] else NULL,
			SE.alpha=if(model=="OU") sqrt(diag(Cov.multi)[1:N+2*N]) else NULL,
			r=if(model=="EB") fit.multi$par[1:N+2*N] else NULL,
			SE.r=if(model=="EB") sqrt(diag(Cov.multi)[1:N+2*N]) else NULL,
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
}

prbt.cont<-function(x,digits=digits){
	N<-x$N
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
		cat(paste("\t",paste(paste("s^2[",1:N,"]",sep=""),collapse="\t"),"\t",
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
		cat(paste("\t",paste(paste("s^2[",1:N,"]",sep=""),collapse="\t"),"\t",
			paste(paste("a[",1:N,"]",sep=""),collapse="\t"),"\t",
			paste(paste("alp[",1:N,"]",sep=""),collapse="\t")),
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
		cat(paste("\t",paste(paste("s^2[",1:N,"]",sep=""),collapse="\t"),"\t",
			paste(paste("a[",1:N,"]",sep=""),collapse="\t"),"\t",
			paste(paste("r[",1:N,"]",sep=""),collapse="\t")),
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
		t<-df<-P<-matrix(0,x$N,x$N)
		for(i in 1:x$N){
			for(j in 1:x$N){
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
	paste("tree",1:(N-1),"vs.",2:N)
	nn<-vector(mode="character",length=N*(N-1)/2)
	k<-1
	for(i in 1:N) for(j in i:N){
		if(i!=j){
			nn[k]<-paste("tree",i," vs. ",j,sep="")
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
