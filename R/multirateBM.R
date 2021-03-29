## function
## by Liam J. Revell 2021

p.func<-function(x,a,C) dmnorm(x,rep(a,nrow(C)),C,log=TRUE)

log_lik<-function(lnsig2,tree,x,trace=0){
	sig2<-exp(lnsig2)
	tt<-tree
	tt$edge.length<-tree$edge.length*apply(tree$edge,1,
		function(e,x) mean(x[e]),x=sig2)
	Tips<-phyl.vcv(as.matrix(x),vcv(tt),1)
	root<-Ntip(tree)+1
	n<-nrow(Tips$C)
	m<-tree$Nnode
	logL<-p.func(x,Tips$alpha[1,1],Tips$C)+
		p.func(lnsig2[-root],lnsig2[root],vcvPhylo(tree))
	if(trace>0){
		cat(paste("log(L) =",round(logL,4),"\n"))
		flush.console()
	}
	-logL
}

log_relik<-function(lnsig2,tree,x,trace=0){
	sig2<-exp(lnsig2)
	tt<-tree
	tt$edge.length<-tree$edge.length*apply(tree$edge,1,
		function(e,x) mean(x[e]),x=sig2)
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
	optim=c("Nelder-Mead","BFGS","CG"),
	maxit=2000,n.iter=1,...){
	method<-method[1]
	if(method=="REML"){
		cat("Sorry! method=\"REML\" doesn't work. Switching to \"ML\".\n")
		method<-"ML"
	}
	lik<-if(method=="ML") log_lik else log_relik
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-0
	x<-x[tree$tip.label]
	fit<-list()
	fit$convergence<-99
	fit$par<-rep(log(mean(pic(x,multi2di(tree))^2)),Ntip(tree)+tree$Nnode)
	ii<-0
	cat("Beginning optimization....\n")
	while(fit$convergence!=0||ii<n.iter){
		cat(paste("Optimization iteration ",ii+1,".\n",sep=""))
		flush.console()
		OPTIM<-optim[if(ii%%length(optim)) ii%%length(optim) else 
			length(optim)]
		fit<-optim(fit$par,
			lik,tree=tree,x=x,trace=trace,
			control=list(maxit=maxit),
			method=OPTIM)
		ii<-ii+1
	}
	cat("Done optimization.\n")
	LIK<-function(sig2) -lik(log(sig2),tree=tree,x=x)
	object<-list(
		sig2=setNames(exp(fit$par),c(tree$tip.label,
		1:tree$Nnode+Ntip(tree))),
		logLik=-fit$value,
		k=length(fit$par)+1,
		tree=tree,
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

plot.multirateBM<-function(x,digits=1,...){
	cols<-setNames(rainbow(1000,start=0.7,end=0),
		1:1000)
	est.sig2<-apply(x$tree$edge,1,function(e,x) 
		mean(x[e]),x=x$sig2)
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
	plotTree(x$tree,plot=FALSE,...)
	lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	obj<-plot(tree,colors=cols,lwd=3,split.vertical=TRUE,
		xlim=c(-0.3,1.05)*diff(lastPP$x.lim),
		add=TRUE,...)
	lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	h<-max(nodeHeights(x$tree))
	LWD<-diff(par()$usr[1:2])/dev.size("px")[1]
	lines(x=rep(-0.25*h+LWD*15/2,2),y=c(2,Ntip(x$tree)-1))
	nticks<-10
	Y<-cbind(seq(2,Ntip(x$tree)-1,length.out=nticks),
    		seq(2,Ntip(x$tree)-1,length.out=nticks))
	X<-cbind(rep(-0.25*h+LWD*15/2,nticks),
    		rep(-0.25*h+LWD*15/2+0.02*h,nticks))
	for(i in 1:nrow(Y)) lines(X[i,],Y[i,])
	ticks<-exp(seq(min(log(x$sig2)),max(log(x$sig2)),
		length.out=nticks))
	add.color.bar(Ntip(x$tree)-3,
		rainbow(1000,start=0.7,end=0),
		title=expression(paste("evolutionary rate ( ",sigma^2,")")),
		lims=NULL,digits=3,
		direction="upwards",
    		subtitle="",lwd=15,
		x=-0.25*h,
		y=2,prompt=FALSE)
	text(x=X[,2],y=Y[,2],signif(ticks,digits),pos=4,cex=0.7)
}

		
