## bounded Brownian model based on Boucher & Demery 2016 and
## wrapped (circular) Brownian model based on Kuhn et al

bounded_bm<-function(tree,x,lims=NULL,...){
	if(is.null(lims)) lims<-range(x)
	if(hasArg(levs)) levs<-list(...)$levs
	else levs<-100
	if(hasArg(wrapped)) wrapped<-list(...)$wrapped
	else wrapped<-FALSE
	if(hasArg(expm.method)) expm.method<-list(...)$expm.method
	else expm.method<-"R_Eigen"
	if(hasArg(lik.func)) lik.func<-list(...)$lik.func
	else lik.func<-"pruning"
	if(hasArg(root)) root=list(...)$root
	else root<-"nuisance"
	if(root=="nuisance") pi<-"fitzjohn"
	else if(root=="mle") pi<-"mle"
	dd<-diff(lims)
	tol<-1e-8*dd/levs
	bins<-cbind(seq(from=lims[1]-tol,by=(dd+2*tol)/levs,length.out=levs),
		seq(to=lims[2]+tol,by=(dd+2*tol)/levs,length.out=levs))
	xx<-setNames(
		sapply(x,function(y,bins) which((y>=bins[,1])+(y<bins[,2])==2),
		bins=bins),names(x))
	X<-to.matrix(xx,1:levs)
	MODEL<-matrix(0,levs,levs,dimnames=list(1:levs,1:levs))
	for(i in 1:(levs-1)) MODEL[i,i+1]<-MODEL[i+1,i]<-1
	if(wrapped) MODEL[1,levs]<-MODEL[levs,1]<-1
	pic_x<-pic(x,multi2di(tree))
	nn<-max(c(floor(0.2*multi2di(tree)$Nnode),1))
	max.q=20*(1/2)*mean(sort(pic_x^2,decreasing=TRUE)[1:nn])*
		(levs/dd)^2
	q.init<-runif(n=1,0,2)*(1/2)*mean(sort(pic_x^2,decreasing=TRUE)[1:nn])*
		(levs/dd)^2
	fit<-fitMk(tree,X,model=MODEL,lik.func=lik.func,pi=pi,
		expm.method=expm.method,logscale=TRUE,max.q=max.q,q.init=q.init)
	lik<-logLik(fit)-Ntip(tree)*log(dd/levs)
	attr(lik,"df")<-2
	class(lik)<-"logLik"
	sigsq<-2*fit$rates*(dd/levs)^2
	ff<-if(lik.func%in%c("pruning","parallel")) lik.func else "pruning"
	x0<-sum(ancr(fit,lik.func=ff,expm.method=expm.method)$ace[1,]*
			rowMeans(bins))
	object<-list(
		wrapped=wrapped,
		sigsq=sigsq,
		x0=x0,
		bounds=lims,
		ncat=levs,
		logLik=lik,
		opt_results=fit$opt_results,
		mk_fit=fit)
	class(object)<-"bounded_bm"
	object
}

print.bounded_bm<-function(x,digits=6,...){
	cat(paste("Object of class \"bounded_bm\" based on\n",
		"	a discretization with k =",
		x$ncat,"levels.\n"))
	if(x$wrapped) cat("\nWrapped (i.e., circular) model.\n\n")
	if(!x$wrapped) cat("\nUnwrapped (i.e., bounded) model\n\n")
	cat(paste("Set or estimated bounds: [",round(x$bounds[1],digits),
		",",round(x$bounds[2],digits),"]\n\n"))
	cat("Fitted model parameters:\n")
	cat(paste("	 sigsq:",round(x$sigsq,6),"\n"))
	cat(paste("     x0:",round(x$x0,6),"\n\n"))
	cat(paste("Log-likelihood:",round(x$logLik,digits),"\n\n"))
	if(x$opt_results$convergence == 0) 
		cat("R thinks it has found the ML solution.\n\n")
	else cat("R thinks optimization may not have converged.\n\n")
}

logLik.bounded_bm<-function(x) x$logLik

ancr.bounded_bm<-function(object,...){
	if(hasArg(lik.func)) lik.func<-list(...)$lik.func
	else lik.func<-"pruning"
	if(hasArg(expm.method)) expm.method<-list(...)$expm.method
	else expm.method<-"R_Eigen"
	dd<-diff(object$bounds)
	tol<-1e-8*dd/object$ncat
	bins<-cbind(seq(from=object$bounds[1]-tol,
		by=(dd+2*tol)/object$ncat,length.out=object$ncat),
		seq(to=object$bounds[2]+tol,by=(dd+2*tol)/object$ncat,
			length.out=object$ncat))
	mids<-rowMeans(bins)
	Anc<-ancr(object$mk_fit,lik.func=lik.func,
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
