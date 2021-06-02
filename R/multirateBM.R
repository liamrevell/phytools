## function
## by Liam J. Revell 2021

VCV<-function(tree){
	H<-nodeHeights(tree)
	h<-c(H[1,1],H[,2])[order(c(tree$edge[1,1],tree$edge[,2]))]
	M<-mrca(tree,full=TRUE)
	M[lower.tri(M)]<-t(M)[lower.tri(M)]
	ROOT<-Ntip(tree)+1
	M<-M[-ROOT,-ROOT]
	matrix(h[M],nrow(M),ncol(M),
		dimnames=list(c(tree$tip.label,2:tree$Nnode+Ntip(tree)),
		c(tree$tip.label,2:tree$Nnode+Ntip(tree))))
}

p.func<-function(x,a,C) dmnorm(x,rep(a,nrow(C)),C,log=TRUE)

	
ln.mean<-function(x){
	if(x[1]==x[2]) return(x[1])
	else {
		a<-x[2]
		b<-log(x[1])-log(x[2])
		return(a/b*exp(b)-a/b)
	}
}

log_lik<-function(lnsig2,tree,x,lambda=1,trace=0){
	sig2<-exp(lnsig2)
	tt<-tree
	tt$edge.length<-tree$edge.length*apply(tree$edge,
		1,function(e,x) ln.mean(x[e]),x=sig2)
	Tips<-phyl.vcv(as.matrix(x),vcv(tt),1)
	root<-Ntip(tree)+1
	n<-nrow(Tips$C)
	m<-tree$Nnode
	ln.p1<-p.func(x,Tips$alpha[1,1],Tips$C)
	ln.p2<--p.func(log(sig2)[-root],log(sig2)[root],
		VCV(tree))
	logL<-ln.p1-lambda*ln.p2
	if(trace>0){
		cat(paste("log(L) =",round(logL,4),"\n"))
		flush.console()
	}
	-logL
}

log_relik<-function(lnsig2,tree,x,trace=0){
	sig2<-exp(lnsig2[1:(length(lnsig2)-1)])
	SCALE<-exp(lnsig2[length(lnsig2)])
	tt<-tree
	tt$edge.length<-SCALE*tree$edge.length*apply(tree$edge,
		1,function(e,x) mean(x[e]),x=sig2)
	px<-pic(x,tt)
	dLNSIG<-apply(tree$edge,1,function(edge,lnsig2) diff(lnsig2[edge]),
		lnsig2=lnsig2)/sqrt(tree$edge.length)
	logL<-sum(dnorm(px,log=TRUE))+sum(dnorm(dLNSIG,log=TRUE))
	if(trace>0){
		cat(paste("log(reL) =",round(logL,4),"\n"))
		flush.console()
	}
	-logL
}

multirateBM<-function(tree,x,method=c("ML","REML"),
	optim=c("L-BFGS-B","Nelder-Mead","BFGS","CG"),
	maxit=NULL,n.iter=1,lambda=1,...){
	method<-method[1]
	optim.method<-optim
	if(!is.null(maxit)) control<-list(maxit=maxit)
	else control<-list()
	if(method=="REML"){
		cat("Sorry! method=\"REML\" doesn't work. Switching to \"ML\".\n")
		method<-"ML"
	}
	if("L-BFGS-B"%in%optim.method){
		vv<-log(var(x)/max(nodeHeights(tree)))
		if(hasArg(lower)){ 
			lower<-list(...)$lower
			lower<-log(lower)
		} else lower<-rep(-10+vv,Ntip(tree)+tree$Nnode)
		if(length(lower)!=(Ntip(tree)+tree$Nnode)) 
			lower<-rep(lower[1],Ntip(tree)+tree$Nnode)
		if(hasArg(upper)){ 
			upper<-list(...)$upper
			lower<-log(lower)
		} else upper<-rep(10+vv,Ntip(tree)+tree$Nnode)
		if(length(upper)!=(Ntip(tree)+tree$Nnode))
			upper<-rep(upper[1],Ntip(tree)+tree$Nnode)
	}	
	lik<-if(method=="ML") log_lik else log_relik
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-0
	x<-x[tree$tip.label]
	fit<-list()
	fit$convergence<-99
	init<-log(mean(pic(x,multi2di(collapse.singles(tree)))^2)*(Ntip(tree)-1)/Ntip(tree))
	fit$par<-rep(init,Ntip(tree)+tree$Nnode)
	class(fit)<-"try-error"
	ii<-1
	cat("Beginning optimization....\n")
	while(inherits(fit,"try-error")||fit$convergence!=0||ii<=n.iter){
		if(length(optim.method)==1) OPTIM<-optim.method[1]
		else OPTIM<-optim.method[if(ii!=length(optim.method)) ii%%length(optim.method) else 
			length(optim.method)]
		cat(paste("Optimization iteration ",ii,". Using \"",
			OPTIM,"\" optimization method.\n",sep=""))
		if(OPTIM=="L-BFGS-B"){
			fit$par[which(fit$par<lower)]<-lower
			fit$par[which(fit$par>upper)]<-upper
		}
		flush.console()
		cur.vals<-fit$par
		fit<-try(optim(fit$par,
			lik,tree=tree,x=x,lambda=lambda,trace=trace,
			control=control,
			method=OPTIM,
			lower=if(OPTIM=="L-BFGS-B") lower else -Inf,
			upper=if(OPTIM=="L-BFGS-B") upper else Inf))
		if(inherits(fit,"try-error")){
			fit<-list()
			fit$convergence<-99
			fit$par<-cur.vals
			class(fit)<-"try-error"
			cat("Caught error without failing. Trying again....\n")
		} else {
			cat(paste("Best (penalized) log-likelihood so far:",signif(-fit$value,6),"\n"))
		}
		ii<-ii+1
	}
	cat("Done optimization.\n")
	LIK<-function(sig2) -lik(log(sig2),tree=tree,x=x,
		lambda=0)
	object<-list(
		sig2=setNames(exp(fit$par),
		c(tree$tip.label,1:tree$Nnode+Ntip(tree))),
		lambda=lambda,
		logLik=LIK(exp(fit$par)),
		k=length(fit$par)+1,
		tree=tree,
		x=x,
		convergence=fit$convergence,
		method=method,
		lik=LIK)
	class(object)<-"multirateBM"
	object
}

print.multirateBM<-function(x,digits=6,printlen=NULL,...){
	if(is.null(printlen)) printlen<-8
	cat("Multi-rate Brownian model using multirateBM.\n\n")
	cat("Fitted rates:\n")
	if(printlen>=length(x$sig2)) print(round(x$sig2,digits))
	else printDotDot(x$sig2,digits,printlen)
	cat(paste("\nlambda penalty term:",
		round(x$lambda,digits)))
	cat(paste("\nlog-likelihood: ",round(x$logLik,digits),"\n"))
	cat(paste("AIC: ",round(AIC(x),digits),"\n\n"))
	if(x$convergence==0) cat("R thinks it has found a solution.\n\n") 
	else cat("R may not have found a solution.\n\n")
}

logLik.multirateBM<-function(object,...){
	lik<-object$logLik
	attr(lik,"df")<-object$k
	lik
}

plot.multirateBM<-function(x,digits=2,...){
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	cols<-setNames(rainbow(1000,start=0.7,end=0),
		1:1000)
	est.sig2<-apply(x$tree$edge,1,function(e,x) 
		ln.mean(x[e]),x=x$sig2)
	ln.sig2<-log(est.sig2)
	min.sig2<-min(ln.sig2)
	max.sig2<-max(ln.sig2)
	edge.states<-vector()
	for(i in 1:length(est.sig2)){
		edge.states[i]<-round((ln.sig2[i]-min.sig2)/
			(max.sig2-min.sig2)*999)+1
	}
	tree<-paintBranches(x$tree,edge=x$tree$edge[1,2],
		state=edge.states[1])
	for(i in 2:length(edge.states))
		tree<-paintBranches(tree,edge=tree$edge[i,2],
			state=edge.states[i])
	nticks<-10
	ticks<-exp(seq(min(log(x$sig2)),max(log(x$sig2)),
		length.out=nticks))
	object<-list(tree=tree,cols=cols,ticks=ticks)
	class(object)<-"multirateBM_plot"
	if(plot) plot(object,digits=digits,...)
	invisible(object)
}

plot.multirateBM_plot<-function(x,digits=2,...){
	args<-list(...)
	if(!is.null(args$type)) if(args$type=="fan") args$type<-"phylogram"
	if(is.null(args$lwd)) args$lwd<-3
	args$split.vertical<-TRUE
	args$tree<-x$tree
	args$plot<-FALSE
	do.call(plotTree,args)
	lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	args$plot<-TRUE
	args$colors<-x$cols
	args$xlim<-c(-0.3,1.05)*diff(lastPP$x.lim)
	args$add<-TRUE
	do.call(plotSimmap,args)
	lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	h<-max(nodeHeights(x$tree))
	LWD<-diff(par()$usr[1:2])/dev.size("px")[1]
	lines(x=rep(-0.25*h+LWD*15/2,2),y=c(2,Ntip(x$tree)-1))
	nticks<-length(x$ticks)
	Y<-cbind(seq(2,Ntip(x$tree)-1,length.out=nticks),
    		seq(2,Ntip(x$tree)-1,length.out=nticks))
	X<-cbind(rep(-0.25*h+LWD*15/2,nticks),
    		rep(-0.25*h+LWD*15/2+0.02*h,nticks))
	for(i in 1:nrow(Y)) lines(X[i,],Y[i,])
	add.color.bar(Ntip(x$tree)-3,
		x$cols,
		title=expression(paste("evolutionary rate ( ",sigma^2,")")),
		lims=NULL,digits=3,
		direction="upwards",
    	subtitle="",lwd=15,
		x=-0.25*h,
		y=2,prompt=FALSE)
	text(x=X[,2],y=Y[,2],signif(x$ticks,digits),pos=4,cex=0.7)
}

print.multirateBM_plot<-function(x,...){
	cat("Object of class \"multirateBM_plot\" containing:\n\n")
	cat(paste("(1) A phylogenetic tree with ",length(x$tree$tip.label)," tips and ",x$tree$Nnode," internal nodes.\n\n",sep=""))
	cat("(2) A mapped set of Brownian evolution rates.\n\n") 
}

setMap.multirateBM_plot<-function(x,...){
	if(hasArg(invert)) invert<-list(...)$invert
	else invert<-FALSE
	n<-length(x$cols)
	if(invert) x$cols<-setNames(rev(x$cols),names(x$cols))
	else x$cols[1:n]<-colorRampPalette(...)(n)
	x
}
