fitfnMk<-function(tree,x,model="polynomial",degree=2,...){
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(model!="polynomial"){
		stop("Sorry. Only available model (so far) is \"polynomial\". Stopping.\n")
	}
	if(hasArg(niter)) niter<-list(...)$niter
	else niter<-10
	if(niter>1){
		if(hasArg(parallel)) parallel<-list(...)$parallel
		else parallel<-FALSE
		if(hasArg(opt.method)) opt.method<-list(...)$opt.method
		else opt.method<-"nlminb"
	} else parallel<-FALSE
	if(niter>1){
		args<-list(...)
		args$niter<-1
		args$tree<-tree
		args$x<-x
		args$degree<-degree
		if(parallel&&(opt.method!="DEoptim")){
			if(hasArg(ncores)) ncores<-list(...)$ncores
			else ncores<-min(c(detectCores()-1,niter))
			mc<-makeCluster(ncores,type="PSOCK")
			registerDoParallel(cl=mc)
			if(!quiet){
				cat(paste("Opened cluster with",ncores,"cores.\n"))
				cat("Running optimization iterations in parallel.\n")
				cat("Please wait....\n")
				flush.console()
			}
			fits<-foreach(i = 1:niter)%dopar%{
				do.call(fitfnMk,args)
			}
			stopCluster(cl=mc)
			lnL<-sapply(fits,logLik)
		} else {
			fits<-list()
			for(i in 1:niter){ 
				fits[[i]]<-do.call(fitfnMk,args)
				lnL<-sapply(fits,logLik)
				if(!quiet){
					cat(paste("log-likelihood from current iteration:", 
						round(logLik(fits[[i]]),4),"\n"))
					cat(paste(" --- Best log-likelihood so far:", 
						round(max(lnL), 4),"---\n"))
						flush.console()
				}
			}
		}
		object<-fits[[which(lnL==max(lnL))[1]]]
		object$all.fits<-fits
	} else {
		if(hasArg(trace)) trace<-list(...)$trace
		else trace<-0
		if(hasArg(start)) start<-list(...)$start
		else start<-NULL
		if(hasArg(opt.method)) opt.method<-list(...)$opt.method
		else opt.method<-"nlminb"
		if(length(degree)==1) degree<-rep(degree,2)
		if(is.matrix(x)){
			levs<-colnames(x)
		} else if(is.numeric(x)){
			levs<-min(x):max(x)
			x<-to.matrix(x,levs)
		} else if(is.factor(x)){
			if(suppressWarnings(all(!is.na(as.numeric(levels(x)))))){
				levs<-min(as.numeric(levels(x))):max(as.numeric(levels(x)))
				x<-to.matrix(x,levs)
			} else {
				levs<-sort(levels(x))
				x<-to.matrix(x,levs)
			}
		} else if(is.character(x)){
			if(suppressWarnings(all(!is.na(as.numeric(x))))){
				levs<-min(as.numeric(x)):max(as.numeric(x))
				x<-to.matrix(x,levs)
			} else {
				levs<-sort(unique(x))
				x<-to.matrix(x,levs)
			}
		}
		x<-x[tree$tip.label,]
		k<-ncol(x)	
		if(hasArg(pi)) pi<-list(...)$pi
		else pi<-"equal"
		if(is.numeric(pi)) root.prior<-"given"
		if(pi[1]=="equal"){ 
			pi<-setNames(rep(1/k,k),levs)
			root.prior<-"flat"
		} else if(pi[1]=="fitzjohn"){ 
			root.prior<-"nuisance"
		} else if(pi[1]=="mle") root.prior<-"it's MLE"	
		lik<-function(par,pw,Y,pi=pi,degree=degree){
			k<-ncol(Y)
			m<-1:(k-1)-0.5
			q1<-rep(0,length(m))
			for(i in 0:degree[1]) q1<-q1+par[i+1]*m^(degree[1]-i)
			q2<-rep(0,length(m))
			for(i in 0:degree[2]) q2<-q2+par[degree[1]+i+2]*m^(degree[2]-i)
			q1[q1<0]<-0
			q2[q2<0]<-0
			if(any(is.nan(c(q1,q2)))){
				return(Inf)
			} else if(all(q1<0)||all(q2<0)){ 
				return(Inf)
			} else {
				MODEL<-matrix(0,k,k,dimnames=list(colnames(Y),colnames(Y)))
				MODEL[cbind(1:(k-1),2:k)]<-1:(k-1)
				MODEL[cbind(2:k,1:(k-1))]<-k:(2*k-2)
				return(-pruning(c(q1,q2),pw,Y,model=MODEL,pi=pi))
			}
		}
		pw<-reorder(tree,"postorder")
		xx<-0:(k-2)+0.5
		if(is.null(start)) q1_start<-q2_start<--1
		else if(start[1]=="smart"){
			cat("\nNumerically optimizing simple equal-rates ordered model\n")
			cat("  to get better random starting values....\n")
			MODEL<-matrix(0,k,k,dimnames=list(colnames(x),colnames(x)))
			MODEL[cbind(1:(k-1),2:k)]<-1
			MODEL[cbind(2:k,1:(k-1))]<-2
			RATES<-fitMk(pw,x,model=MODEL,pi=pi,lik.func="pruning")$rates
			start<-rep(0,sum(degree)+2)
			start[c(degree[1]+1,sum(degree)+2)]<-RATES
			start<-start+runif(n=sum(degree)+2,min=-0.0001*mean(RATES),
				max=0.0001*mean(RATES))
			q1_start<-q2_start<-rep(0,length(xx))
			for(i in 0:degree[1]) q1_start<-q1_start+start[i+1]*xx^(degree[1]-i)
			for(i in 0:degree[2]) q2_start<-q2_start+start[degree[1]+i+2]*xx^(degree[2]-i)
			q1_start[q1_start<0]<-0
			q2_start[q2_start<0]<-0
			if(all(q1_start==0)&&all(q2_start==0)) q1_start<-q2_start<--1	
		} else if(start[1]=="really smart"){
			cat("\nNumerically optimizing ordered all-rates-different model\n")
			cat("  to get better random starting values....\n")
			MODEL<-matrix(0,k,k,dimnames=list(colnames(x),colnames(x)))
			MODEL[cbind(1:(k-1),2:k)]<-1:(k-1)
			MODEL[cbind(2:k,1:(k-1))]<-1:(k-1)+(k-1)
			RATES<-fitMk(pw,x,model=MODEL,pi=pi,lik.func="pruning",rand_start=TRUE)$rates
			lm1<-lm(RATES[1:(k-1)]~poly(xx,degree=degree[1],raw=TRUE))
			lm2<-lm(RATES[1:(k-1)+(k-1)]~poly(xx,degree=degree[2],raw=TRUE))
			start<-c(unclass(coef(lm1))[(degree[1]+1):1],
				unclass(coef(lm2))[(degree[2]+1):1])
			names(start)<-NULL
			q1_start<-q2_start<-rep(0,length(xx))
			for(i in 0:degree[1]) q1_start<-q1_start+start[i+1]*xx^(degree[1]-i)
			for(i in 0:degree[2]) q2_start<-q2_start+start[degree[1]+i+2]*xx^(degree[2]-i)
			q1_start[q1_start<0]<-0
			q2_start[q2_start<0]<-0
		} else if(is.numeric(start)){
			q1_start<-q2_start<-rep(0,length(xx))
			for(i in 0:degree[1]) q1_start<-q1_start+start[i+1]*xx^(degree[1]-i)
			for(i in 0:degree[2]) q2_start<-q2_start+start[degree[1]+i+2]*xx^(degree[2]-i)
			q1_start[q1_start<0]<-0
			q2_start[q2_start<0]<-0
			if(all(q1_start==0)&&all(q2_start==0)) q1_start<-q2_start<--1
		} else q1_start<-q2_start<--1
		while(any(q1_start<0)||any(q2_start<0)){
			start<-runif(n=sum(degree)+2)
			q1_start<-q2_start<-rep(0,length(xx))
			for(i in 0:degree[1]) q1_start<-q1_start+start[i+1]*xx^(degree[1]-i)
			for(i in 0:degree[2]) q2_start<-q2_start+start[degree[1]+i+2]*xx^(degree[2]-i)
		}
		if(opt.method=="nlminb"){
			fit<-nlminb(start,lik,pw=pw,Y=x,pi=pi,degree=degree,
				control=list(trace=trace))
		} else if(opt.method=="DEoptim"){
			if(hasArg(maxit)) maxit<-list(...)$maxit
			else maxit<-1000
			if(hasArg(prompt)) prompt<-list(...)$prompt
			else prompt<-FALSE
			if(hasArg(parallel)) parallel<-list(...)$parallel
			else parallel<-FALSE
			if(hasArg(trace)) trace<-list(...)$trace
			else trace<-if(maxit>100) 100 else 1
			parallelType<-if(parallel) 1 else 0
			if(hasArg(NP)) NP<-list(...)$NP
			else NP<-100*length(start)
			if(hasArg(F)) F<-list(...)$F
			else F<-0.8
			if(hasArg(CR)) CR<-list(...)$CR
			else CR<-0.9
			initialpop<-matrix(rep(start,NP)*runif(NP*length(start)),
				NP,length(start),byrow=TRUE)
			if(parallel){
				if(hasArg(ncores)) ncores<-list(...)$ncores
				else ncores<-detectCores()-1
				mc<-makeCluster(ncores,type="PSOCK")
				registerDoParallel(cl=mc)
			}
			if(trace>0) cat("\n")
			fit<-DEoptim(lik,
				lower=rep(-10*max(start),length(start)),
				upper=rep(10*max(start),length(start)),
				pw=pw,Y=x,pi=pi,degree=degree,
				control=list(itermax=maxit,
				initialpop=initialpop,
				cluster=if(parallel) mc else NULL,
				trace=trace,F=F,CR=CR,NP=NP))
			if(prompt){ 
				cat("\nMax iterations (itermax) reached. Continue? (yes/no):")
				continue<-readline()
			} else continue<-"no"
			while(continue=="yes"&&prompt){
				fit<-DEoptim(lik,
					lower=rep(-10*max(start),length(start)),
					upper=rep(10*max(start),length(start)),
					pw=pw,Y=x,pi=pi,degree=degree,
					control=list(itermax=maxit,
					initialpop=fit$member$pop,
					cluster=if(parallel) mc else NULL,
					trace=trace,F=F,CR=CR,NP=NP))
				cat("\nMax iterations (itermax) reached. Continue? (yes/no):")
				continue<-readline()
			}
			if(trace>0) cat("\n")
			stopCluster(cl=mc)
			fit<-list(
				par=fit$optim$bestmem,
				objective=fit$optim$bestval,
				iterations=fit$optim$iter,
				evaluations=fit$optim$nfeval)
		}
		q1_est<-rep(0,length(xx))
		for(i in 0:degree[1]) q1_est<-q1_est+fit$par[i+1]*xx^(degree[1]-i)
		q2_est<-rep(0,length(xx))
		for(i in 0:degree[2]) q2_est<-q2_est+fit$par[degree[1]+i+2]*xx^(degree[2]-i)
		q1_est[q1_est<0]<-0
		q2_est[q2_est<0]<-0
		index.matrix<-matrix(0,k,k,dimnames=list(colnames(x),colnames(x)))
		index.matrix[cbind(1:(k-1),2:k)]<-1:(k-1)
		index.matrix[cbind(2:k,1:(k-1))]<-k:(2*k-2)
		lik.f<-function(par) -lik(par,pw=pw,Y=x,pi=pi,degree=degree)
		object<-list(
			logLik=-fit$objective,
			rates=c(q1_est,q2_est),
			index.matrix=index.matrix,
			states=levs,
			pi=pi,
			method=opt.method,
			root.prior=root.prior,
			opt_results=fit[c("convergence","iterations","evaluations","message")],
			par=fit$par,
			degree=degree,
			data=x,
			tree=tree,
			lik=lik.f)
		class(object)<-c("fitfnMk","fitMk")
	}	
	object
}

plot.fitfnMk<-function(x,...){
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-c(5.1,4.1,4.1,4.1)
	if(hasArg(log)) log<-list(...)$log
	else log=""
	if(hasArg(leg_pos)) leg_pos<-list(...)$leg_pos
	else leg_pos<-"topright"
	k<-length(x$states)
	q1<-x$rates[1:(k-1)]
	q2<-x$rates[k:(2*k-2)]
	xx<-seq(0.5,(k-2)+0.5,length.out=100)
	## compute actual polynomial function
	X<-seq(0.5,(k-0.5),length.out=100)
	qq1<-rep(0,100)
	for(i in 0:x$degree[1]) qq1<-qq1+x$par[i+1]*xx^(x$degree[1]-i)
	qq1[qq1<0]<-0
	qq2<-rep(0,100)
	for(i in 0:x$degree[2]) qq2<-qq2+x$par[x$degree[1]+i+2]*xx^(x$degree[2]-i)
	qq2[qq2<0]<-0
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-range(pretty(c(0,max(x$rates))))
	par(mar=mar)
	plot(qq1,xx,type="l",col="blue",bty="n",las=1,
		axes=FALSE,xlab="transition rate (q)",ylab="",
		xlim=xlim,ylim=c(min(xx),1.05*max(xx)),log=log)
	points(q1,0:(k-2)+0.5,pch=16,
		col=if(par()$bg=="transparent") "white" else par()$bg,cex=2)
	lines(qq2,xx,type="l",col="red")
	points(q2,0:(k-2)+0.5,pch=16,
		col=if(par()$bg=="transparent") "white" else par()$bg,cex=2)
	points(q1,0:(k-2)+0.5,pch=16,col="blue",cex=1)
	points(q2,0:(k-2)+0.5,pch=16,col="red",cex=1)	
	labs<-mapply(function(x,y) bquote(.(x) %->% .(y)),
		x=x$states[1:(k-1)],y=x$states[2:k])
	axis(2,at=seq(0.5,k-1.5,by=1),labels=rep("",k-1))
	nulo<-mapply(mtext,text=labs,at=seq(0.5,k-1.5,by=1),
		MoreArgs=list(side=2,line=1,las=3,cex=0.7))
	axis(4,at=seq(0.5,k-1.5,by=1),labels=rep("",k-1))
	labs<-mapply(function(x,y) bquote(.(x) %<-% .(y)),
		x=x$states[1:(k-1)],y=x$states[2:k])
	nulo<-mapply(mtext,text=labs,at=seq(0.5,k-1.5,by=1),
		MoreArgs=list(side=4,line=1,las=3,cex=0.7))
	axis(1,las=1,cex.axis=0.8,
		at=pretty(c(0,max(x$rates))))
	clip(par()$usr[1],par()$usr[2],par()$usr[3],max(xx))
	grid()
	do.call("clip",as.list(par()$usr))
	legend(leg_pos,c("forward","backward"),
		col=c("blue","red"),
		lty="solid",pch=16,cex=0.8,pt.cex=1,
		bg=make.transparent("white",0.75))
}

logLik.fitfnMk<-function(object,...){
	lik<-object$logLik
	attr(lik,"df")<-length(object$par)
	lik
}

## print method for objects of class "fitMk"
print.fitfnMk<-function(x,digits=6,...){
	cat("Object of class \"fitfnMk\".\n\n")
	cat("Fitted (or set) value of Q:\n")
	Q<-matrix(NA,length(x$states),length(x$states))
	Q[]<-c(0,x$rates)[x$index.matrix+1]
	diag(Q)<-0
	diag(Q)<--rowSums(Q)
	colnames(Q)<-rownames(Q)<-x$states
	print(round(Q,digits))
	cat("\nFitted (or set) value of pi:\n")
	print(round(x$pi,digits))
	cat(paste("due to treating the root prior as (a) ",x$root.prior,".\n",
		sep=""))
	cat(paste("\nLog-likelihood:",round(x$logLik,digits),"\n"))
	cat(paste("\nOptimization method used was \"",x$method,"\"\n\n",
		sep=""))
	if(!is.null(x$opt_results$convergence)){
		if(x$opt_results$convergence==0) 
			cat("R thinks it has found the ML solution.\n\n")
		else cat("R thinks optimization may not have converged.\n\n")
	}
}
