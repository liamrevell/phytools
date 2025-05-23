## bounded Brownian model based on Boucher & Demery 2016 and
## wrapped (circular) Brownian model based on Juhn et al. (in review)

to_binned<-function(x,bins){
	xx<-setNames(
		sapply(x,function(y,bins) which((y>=bins[,1])+(y<bins[,2])==2),
		bins=bins),names(x))
	to.matrix(xx,1:nrow(bins))
}

bounded_bm<-function(tree,x,lims=NULL,...){
	if(hasArg(wrapped)) wrapped<-list(...)$wrapped
	else wrapped<-FALSE
	if(hasArg(absorbing)) absorbing<-list(...)$absorbing
	else absorbing<-FALSE
	if(wrapped&&absorbing){
		cat("wrapped = TRUE and absorbing = TRUE not permitted.\n")
		cat("Setting absorbing = FALSE.\n\n")
		absorbing<-FALSE
	}
	if(hasArg(bins)) bins<-list(...)$bins
	else bins<-NULL
	if(is.null(lims)){ 
		lims<-expand.range(x,p=2)
		df<-2
	} else if(wrapped){ 
		df<-2
	} else df<-4
	if(lims[1]==-Inf){ 
		lims[1]<-expand.range(x,p=2)[1]
		df<-df-1
	}
	if(lims[2]==Inf){
		lims[2]<-expand.range(x,p=2)[2]
		df<-df-1
	}
	## overrides above if user-supplied
	if(hasArg(df)) df<-list(...)$df
	if(hasArg(levs)) levs<-list(...)$levs
	else levs<-200
	if(hasArg(expm.method)) expm.method<-list(...)$expm.method
	else expm.method<-"R_Eigen"
	if(hasArg(lik.func)) lik.func<-list(...)$lik.func
	else lik.func<-if(absorbing) "pruning" else "eigen"
	if(absorbing&&lik.func=="eigen") lik.func<-"pruning"
	if(hasArg(parallel)) parallel<-list(...)$parallel
	else parallel<-FALSE
	if(hasArg(root)) root=list(...)$root
	else root<-"mle"
	if(root=="nuisance") pi<-"fitzjohn"
	else if(root=="mle") pi<-"mle"
	else if(root=="flat") pi<-rep(1/levs,levs)
	dd<-diff(lims)
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-8*dd/levs
	if(hasArg(CI)) CI<-list(...)$CI
	else CI<-FALSE
	if(CI==TRUE) CI<-0.95
	pic_x<-pic(x,multi2di(tree))
	nn<-max(c(floor(0.2*multi2di(tree)$Nnode),1))
	if(hasArg(max.sigsq)) max.sigsq<-list(...)$max.sigsq
	else max.sigsq<-2*mean(sort(pic_x^2,decreasing=TRUE)[1:nn])
	bins<-if(is.null(bins)) cbind(seq(from=lims[1]-tol,
		by=(dd+2*tol)/levs,length.out=levs),seq(to=lims[2]+tol,
		by=(dd+2*tol)/levs,length.out=levs)) else bins
	X<-to_binned(x,bins)
	MODEL<-matrix(0,levs,levs,dimnames=list(1:levs,1:levs))
	for(i in 1:(levs-1)) MODEL[i,i+1]<-MODEL[i+1,i]<-1
	if(wrapped) MODEL[1,levs]<-MODEL[levs,1]<-1
	else if(absorbing) MODEL[1,2]<-MODEL[levs,levs-1]<-0
	max.q=(1/2)*max.sigsq*(levs/dd)^2
	if(lik.func%in%c("pruning","parallel")){
		q.init<-runif(n=1,0,2)*(1/2)*mean(sort(pic_x^2,decreasing=TRUE)[1:nn])*
			(levs/dd)^2
		fit<-fitMk(tree,X,model=MODEL,lik.func=lik.func,pi=pi,
			expm.method=expm.method,logscale=TRUE,max.q=max.q,
			q.init=q.init)
	} else if(lik.func=="eigen"){
		QQ<-MODEL
		diag(QQ)<--rowSums(MODEL)
		eQQ<-eigen(QQ)
		pw<-reorder(tree,"postorder")
		if(parallel){
			if(hasArg(ncores)) ncores<-list(...)$ncores
			else ncores<-min(nrow(tree$edge),detectCores()-1)
			mc<-makeCluster(ncores,type="PSOCK")
			registerDoParallel(cl=mc)
		}
		fit<-optimize(eigen_pruning,c(tol,max.q),tree=pw,
			x=X,eigenQ=eQQ,parallel=parallel,pi=pi,maximum=TRUE)
		fit<-list(
			logLik=fit$objective,
			rates=fit$maximum,
			index.matrix=MODEL,
			states=colnames(X),
			pi=eigen_pruning(fit$maximum,pw,X,eQQ,pi=pi,return="pi",
				parallel=parallel),
			method="optimize",
			root.prior=if(pi[1]=="fitzjohn") "nuisance" else pi,
			opt_results=list(convergence=0),
			data=X,
			tree=pw)
		class(fit)<-"fitMk"
	}
	lik<-logLik(fit)-Ntip(tree)*log(dd/levs)
	attr(lik,"df")<-df
	class(lik)<-"logLik"
	sigsq<-2*fit$rates*(dd/levs)^2
	ff<-if(lik.func%in%c("pruning","parallel","eigen")) lik.func else "pruning"
	x0<-sum(ancr(fit,lik.func=ff,expm.method=expm.method,parallel=parallel)$ace[1,]*
		rowMeans(bins))
	if(parallel) stopCluster(cl=mc)
	lfunc<-function(sig2,x0="nuisance",...){
		if(hasArg(lik.func)) lik.func<-list(...)$lik.func
		else lik.func<-"pruning"
		if(hasArg(parallel)) parallel<-list(...)$parallel
		else parallel<-FALSE
		q<-(sig2/2)*(levs/dd)^2
		if(x0=="nuisance") pi<-"fitzjohn"
		else if(is.numeric(x0)) pi<-to_binned(x0,bins)[1,]
		if(lik.func=="pruning"){
			lnL<-pruning(q,tree,X,MODEL,pi=pi,expm.method=expm.method)-
				Ntip(tree)*log(dd/levs)
		} else if(lik.func=="parallel"){
			if(!exists("ncores")) ncores<-min(nrow(tree$edge),detectCores()-1)
			mc<-makeCluster(ncores,type="PSOCK")
			registerDoParallel(cl=mc)
			lnL<-parallel_pruning(q,tree,X,MODEL,pi=pi,
				expm.method=expm.method)-Ntip(tree)*log(dd/levs)
			stopCluster(cl=mc)
		} else if(lik.func=="eigen"){
			if(parallel){
				mc<-makeCluster(ncores,type="PSOCK")
				registerDoParallel(cl=mc)
			}
			Q<-MODEL
			diag(Q)<--rowSums(MODEL)
			eQ<-eigen(Q)
			lnL<-eigen_pruning(q,tree,X,eQ,parallel=parallel,pi=pi)-
				Ntip(tree)*log(dd/levs)
			if(parallel) stopCluster(cl=mc)
		}
		lnL		
	}
	if(CI!=FALSE){
		foo<-function(ln_sigsq) lfunc(exp(ln_sigsq))
		hh<-hessian(foo,log(sigsq))
		vv<--1/hh
		if(vv<=0) ci<-c(0,Inf)
		else ci_sigsq<-exp(qnorm(c((1-CI)/2,1-(1-CI)/2),
			mean=log(sigsq),sd=sqrt(vv)))
		attr(ci_sigsq,"prob")<-CI
	} else ci_sigsq<-NULL
	object<-list(
		wrapped=wrapped,
		absorbing=absorbing,
		sigsq=sigsq,
		CI=ci_sigsq,
		x0=x0,
		bounds=lims,
		ncat=levs,
		logLik=lik,
		opt_results=fit$opt_results,
		at_bounds=(sigsq>=(max.sigsq-tol)),
		mk_fit=fit,
		lik=lfunc)
	class(object)<-"bounded_bm"
	object
}

eigen_pruning<-function(q,tree,x,eigenQ,...){
	if(hasArg(return)) return<-list(...)$return
	else return<-"likelihood"
	if(hasArg(parallel)) parallel<-list(...)$parallel
	else parallel<-FALSE
	pw<-if(!is.null(attr(tree,"order"))&&
		attr(tree,"order")=="postorder") tree else 
		reorder(tree,"postorder")
	k<-ncol(x)
	if(hasArg(pi)) pi<-list(...)$pi
	else pi<-rep(1/k,k)
	L<-rbind(x[pw$tip.label,],
		matrix(0,pw$Nnode,k,
		dimnames=list(1:pw$Nnode+Ntip(pw))))
	nn<-unique(pw$edge[,1])
	pp<-vector(mode="numeric",length=length(nn))
	root<-min(nn)
	V<-eigenQ$vectors
	Vi<-t(V)
	Vals<-eigenQ$values
	Expm<-function(t,q) Re(V%*%(exp(q*t*Vals)*Vi))
	if(parallel){
		P.all<-foreach(i=1:nrow(pw$edge))%dopar%{ 
			Expm(pw$edge.length[i],q)
		}
	}
	for(i in 1:length(nn)){
		ee<-which(pw$edge[,1]==nn[i])
		PP<-matrix(NA,length(ee),k)
		for(j in 1:length(ee)){
			if(parallel) P<-P.all[[ee[j]]]
			else P<-Expm(q=q,t=pw$edge.length[ee[j]])
			PP[j,]<-P%*%L[pw$edge[ee[j],2],]
		}
		L[nn[i],]<-apply(PP,2,prod)
		if(nn[i]==root){
			if(pi[1]=="fitzjohn") pi<-L[nn[i],]/sum(L[nn[i],])
			else if(pi[1]=="mle") pi<-as.numeric(L[nn[i],]==max(L[nn[i],]))
			L[nn[i],]<-pi*L[nn[i],]
		}
		pp[i]<-sum(L[nn[i],])
		L[nn[i],]<-L[nn[i],]/pp[i]
	}
	prob<-sum(log(pp))
	if(return=="likelihood") 
		if(is.na(prob)||is.nan(prob)) 
			return(-Inf) else return(prob)
	else if(return=="conditional") L
	else if(return=="pi") pi
}

print.bounded_bm<-function(x,digits=6,...){
	cat(paste("Object of class \"bounded_bm\" based on\n",
		"  a discretization with k =",
		x$ncat,"levels.\n"))
	if(x$wrapped) cat("\nWrapped (i.e., circular) model.\n\n")
	if(!x$wrapped){ 
		if(!x$absorbing) cat("\nUnwrapped (i.e., bounded) model\n\n")
		if(x$absorbing) cat("\nAbsorbing bounded model\n\n")
	}
	cat(paste("Set or estimated bounds: [",round(x$bounds[1],digits),
		",",round(x$bounds[2],digits),"]\n\n"))
	cat("Fitted model parameters:\n")
	if(is.null(x$CI)) cat(paste("  sigsq:",round(x$sigsq,6),"\n"))
	else cat(paste("  sigsq:",round(x$sigsq,6),"  [",round(x$CI[1],4),
		",",round(x$CI[2],4),"]\n"))
	cat(paste("     x0:",round(x$x0,6),"\n\n"))
	cat(paste("Log-likelihood:",round(x$logLik,digits),"\n\n"))
	if(x$opt_results$convergence == 0 && !x$at_bounds) 
		cat("R thinks it has found the ML solution.\n\n")
	else if(x$at_bounds) cat("Optimization may be at bounds.\n\n")
	else cat("R thinks optimization may not have converged.\n\n")
}

logLik.bounded_bm<-function(object,...) object$logLik

ancr.bounded_bm<-function(object,...){
	if(hasArg(lik.func)) lik.func<-list(...)$lik.func
	else lik.func<-"pruning"
	if(hasArg(expm.method)) expm.method<-list(...)$expm.method
	else expm.method<-"R_Eigen"
	if(hasArg(parallel)) parallel<-list(...)$parallel
	else parallel<-FALSE
	dd<-diff(object$bounds)
	tol<-1e-8*dd/object$ncat
	bins<-cbind(seq(from=object$bounds[1]-tol,
		by=(dd+2*tol)/object$ncat,length.out=object$ncat),
		seq(to=object$bounds[2]+tol,by=(dd+2*tol)/object$ncat,
			length.out=object$ncat))
	mids<-rowMeans(bins)
	Anc<-ancr(object$mk_fit,lik.func=lik.func,parallel=parallel,
		expm.method=expm.method)
	ace<-colSums(apply(Anc$ace,1,function(x,y) x*y,
		y=mids))
	ci<-matrix(NA,length(ace),2,dimnames=list(names(ace),
		c("lower","upper")))
	for(i in 1:nrow(Anc$ace)){
		cumprob<-cumsum(Anc$ace[i,])
		ci[i,1]<-mids[which(cumprob>0.025)][1]
		ci[i,2]<-mids[which(cumprob>0.975)][1]
	}
	result<-list(
		ace=ace,
		CI95=ci)
	class(result)<-"ancr.bounded_bm"
	result
}

print.ancr.bounded_bm<-function(x,digits=6,printlen=6,...){
	cat("Ancestral character estimates from \"bounded_bm\" object:\n")
	Nnode<-length(x$ace)
	if(is.null(printlen)||printlen>=Nnode) print(round(x$ace,digits))
	else printDotDot(x$ace,digits,printlen)
	cat("\nLower & upper 95% CIs:\n")
	colnames(x$CI95)<-c("lower","upper")
	if(is.null(printlen)||printlen>=Nnode) print(round(x$CI95,digits))
	else printDotDot(x$CI95,digits,printlen)
	cat("\n")
}

