fitcontMk<-function(tree,x,y,model="sigmoidal",...){

	## reorder tree post-order
	pw<-reorder(tree,"postorder")
	
	## discretize y
	if(hasArg(levs)) levs<-list(...)$levs
	else levs<-100
	y<-y[pw$tip.label]
	lims<-expand.range(y,p=1.5)
	dd<-diff(lims)
	delta<-dd/levs
	tol<-1e-8*delta
	bins<-cbind(seq(from=lims[1]-tol,by=(dd+2*tol)/levs,
		length.out=levs),seq(to=lims[2]+tol,by=(dd+2*tol)/levs,
			length.out=levs))
	Y<-to_binned(y,bins)
	
	## convert x to matrix
	if(is.matrix(x)){
		X<-x[pw$tip.label,]
		m<-ncol(X)
		states<-colnames(X)
	} else {
		X<-to.matrix(x,sort(unique(x)))
		X<-X[pw$tip.label,]
		m<-ncol(X)
		states<-colnames(X)
	}
	
	## combined discretized y & x
	nn<-x_by_y(colnames(X),colnames(Y))
	XY<-matrix(0,nrow=nrow(X),ncol=ncol(X)*ncol(Y),
		dimnames=list(rownames(X),nn))
	for(i in 1:ncol(X)){
		for(j in 1:nrow(Y)){
			XY[j,1:levs+(i-1)*levs]<-X[j,i]*Y[j,]
		}
	}
	
	## likelihood function
	lfunc<-function(theta,plot_model=FALSE,trace=0,parallel=TRUE){
		sig2<-theta[1]
		q0.f<-theta[2]
		q1.f<-theta[3]
		B.f<-theta[4]
		M.f<-theta[5]
		q0.b<-theta[6]
		q1.b<-theta[7]
		B.b<-theta[8]
		M.b<-theta[9]
		INDEX<-matrix(0,nrow=ncol(XY),ncol=ncol(XY),
			dimnames=list(nn,nn))
		for(i in 1:(levs-1))
			for(j in 0:(ncol(X)-1))
				INDEX[i+j*levs,i+1+j*levs]<-
					INDEX[i+1+j*levs,i+j*levs]<-1
		xx<-rowMeans(bins)
		q.f<-q0.f+(q1.f-q0.f)/(1+exp(-B.f*(xx-M.f)))
		q.b<-q0.b+(q1.b-q0.b)/(1+exp(-B.b*(xx-M.b)))
		for(i in 1:levs){
			INDEX[i,i+levs]<-i+1
			INDEX[i+levs,i]<-i+levs+1
		}
		if(parallel)
			lik<-parallel_pruning(
				q=c(sig2/(2*delta^2),q.f,q.b),
				tree=pw,x=XY,model=INDEX,
				pi="mle")-Ntip(pw)*log(delta)
		else 
			lik<-pruning(
				q=c(sig2/(2*delta^2),q.f,q.b),
				tree=pw,x=XY,model=INDEX,
				pi="mle")-Ntip(pw)*log(delta)
		if(plot_model){
			par(mfrow=c(2,1),mar=c(5.1,4.1,1.1,1.1))
			forward<-q.f
			plot(rowMeans(bins),forward,type="l",
				xlab="continuous trait",
				ylab="transition rate (0 to 1)",
				las=1,bty="n",cex.axis=0.6)
			box(col="grey")
			grid()
			backward<-q.b
			plot(rowMeans(bins),backward,type="l",
				xlab="continuous trait",
				ylab="transition rate (1 to 0)",
				las=1,bty="n",cex.axis=0.6)
			box(col="grey")
			grid()
		}
		if(trace==1) cat(sample(".",".\n",prob=c(0.9,0.1)))
		else if(trace==2){
			cat(paste("sig2 ","q0.f ","q1.f ","B.f  ","M.f  ",
				"q0.b ","q1.b ","B.b  ","M.b  ","\n",sep="\t"))
			cat(paste(round(theta,3),collapse="\t"))
			cat(paste("\nlog(L) = ",round(lik,4),"\n\n"))
		}
		return(-lik)
	}
	
	## get reasonable starting values for optimization
	
	## starting value of sig2
	sig2<-sum(pic(y,multi2di(tree))^2)/Ntip(tree)
	
	## subsample x for only the lowest values of y
	nm<-names(which(y<=mean(y)))
	low<-fitMk(keep.tip(tree,nm),x[nm],model="ARD",
		pi="mle")
	## repeat for highest values of y
	nm<-names(which(y>=mean(y)))
	high<-fitMk(keep.tip(tree,nm),x[nm],model="ARD",
		pi="mle")
	
	## set q0 and q1
	q0.f<-low$rates[2]
	q1.f<-high$rates[2]
	q0.b<-low$rates[1]
	q1.b<-high$rates[1]
	
	## set B and M
	B.f<-5.89/diff(range(y))
	M.f<-mean(y)
	B.b<-5.89/diff(range(y))
	M.b<-mean(y)
	
	## starting parameter values for optimization
	theta<-c(sig2,
		q0.f,q1.f,B.f,M.f,
		q0.b,q1.b,B.b,M.b)
	
	## rand_start
	if(hasArg(rand_start)) rand_start<-list(...)$rand_start
	else rand_start<-FALSE
	if(rand_start)
		theta<-theta*runif(n=length(theta),min=0,max=2)
	
	## optimize model
	if(hasArg(parallel)) parallel<-list(...)$parallel
	else parallel<-TRUE
	if(hasArg(plot_model)) plot_model<-list(...)$plot_model
	else plot_model<-FALSE
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-0
	if(hasArg(maxit)) maxit<-list(...)$maxit
	else maxit<-2000
	if(hasArg(opt.method)) opt.method<-list(...)$opt.method
	else opt.method<-"nlminb"

	## set up to parallelize pruning
	if(hasArg(ncores)) ncores<-list(...)$ncores
	else ncores<-max(detectCores()-2,1)
	
	if(ncores==1) parallel<-FALSE
	
	if(parallel){
		mc<-makeCluster(ncores,type="PSOCK")
		registerDoParallel(cl=mc)
	}
	
	## tolerance for starting theta
	tol<-1e-8
	
	## set limits
	lower<-c(tol,
		tol,tol,0,min(y)-dd/2,
		tol,tol,0,min(y)-dd/2)
	upper<-c(Inf,
		Inf,Inf,200/dd,max(y)+dd/2,
		Inf,Inf,200/dd,max(y)+dd/2)
	
	## run optimizaton
	if(opt.method=="optim"){
		fit_optim<-optim(theta,lfunc,
			plot_model=plot_model,
			trace=trace,parallel=parallel,
			method="L-BFGS-B",
			lower=lower,upper=upper,
			control=list(maxit=maxit))
		lnL<--fit_optim$value
		attr(lnL,"df")<-length(theta)
		class(lnL)<-"logLik"
		obj<-list(
			states=colnames(X),
			sigsq=fit_optim$par[1],
			q01=fit_optim$par[c(2,3)],
			q10=fit_optim$par[c(6,7)],
			B=fit_optim$par[c(4,8)],
			M=fit_optim$par[c(5,9)],
			bounds=lims,
			ncat=levs,
			opt_results=list(
				counts=fit_optim$counts,
				convergence=fit_optim$convergence,
				message=fit_optim$message),
			opt.method=opt.method,
			logLik=lnL,
			lik=lfunc)
		class(obj)<-"fitcontMk"
	} else if(opt.method=="nlminb"){
		fit_nlminb<-nlminb(theta,lfunc,
			plot_model=plot_model,
			trace=trace,parallel=parallel,
			control=list(iter.max=maxit),
			lower=lower,upper=upper)
		lnL<--fit_nlminb$objective
		attr(lnL,"df")<-length(theta)
		class(lnL)<-"logLik"
		obj<-list(
			states=colnames(X),
			sigsq=fit_nlminb$par[1],
			q01=fit_nlminb$par[c(2,3)],
			q10=fit_nlminb$par[c(6,7)],
			B=fit_nlminb$par[c(4,8)],
			M=fit_nlminb$par[c(5,9)],
			bounds=lims,
			ncat=levs,
			opt_results=list(
				counts=fit_nlminb$iterations,
				convergence=fit_nlminb$convergence,
				evaluations=fit_nlminb$evaluations,
				message=fit_nlminb$message),
			opt.method=opt.method,
			logLik=lnL,
			lik=lfunc)
		class(obj)<-"fitcontMk"
	}
	if(parallel) stopCluster(cl=mc)
	return(obj)
}

plot.fitcontMk<-function(x,...){
	xx<-seq(x$bounds[1],x$bounds[2],length.out=100)
	q.f<-x$q01[1]+(x$q01[2]-x$q01[1])/(1+exp(-x$B[1]*(xx-x$M[1])))
	q.b<-x$q10[1]+(x$q10[2]-x$q10[1])/(1+exp(-x$B[2]*(xx-x$M[2])))
	par(mfrow=c(2,1),mar=c(5.1,4.1,1.1,1.1))
	ylab<-bquote(q(.(x$states[1])*"\u2192"*.(x$states[2])))
	usr<-list()
	mfg<-list()
	plot(xx,q.f,type="l",
		xlab="continuous trait",
		ylab=ylab,
		las=1,bty="n",cex.axis=0.6,
		ylim=c(0,max(c(x$q01,x$q10))))
	usr[[1]]<-par()$usr
	mfg[[1]]<-par()$mfg
	box(col="grey")
	grid()
	ylab<-bquote(q(.(x$states[2])*"\u2192"*.(x$states[1])))
	plot(xx,q.b,type="l",
		xlab="continuous trait",
		ylab=ylab,
		las=1,bty="n",cex.axis=0.6,
		ylim=c(0,max(c(x$q01,x$q10))))
	usr[[2]]<-par()$usr
	mfg[[2]]<-par()$mfg
	box(col="grey")
	grid()
	invisible(list(usr=usr,mfg=mfg))
}

print.fitcontMk<-function(x,digits=4,...){
	cat(paste(
		"Object of class \"fitcontMk\" based on\n",
		"   a discretization with k =",
		x$ncat,"levels.\n\n"))
	cat(paste("Fitted sigmoidal continuous-trait dependent\n",
		"   M2 model parameters:\n\n"))
	cat(paste(" states: [",
		paste(x$states,collapse=", "),
		"]\n"))
	cat(paste("  sigsq: ",
		round(x$sigsq,digits),"\n"))
	cat(paste("    q01: [",
		paste(round(x$q01,digits),collapse=", "),
		"]\n"))
	cat(paste("    q10: [",
		paste(round(x$q10,digits),collapse=", "),
		"]\n"))
	cat(paste("      B: [",
		paste(round(x$B,digits),collapse=", "),
		"]\n"))
	cat(paste("      M: [",
		paste(round(x$M,digits),collapse=", "),
		"]\n"))
	cat(paste("\nlog(L): ",round(x$logLik,digits), 
		" (df = ",attr(x$logLik,"df"),")\n\n",sep=""))
	if(x$opt_results$convergence == 0) 
		cat("R thinks it has found the ML solution.\n\n")
	else cat(
		"R thinks optimization may not have converged.\n\n")
}
