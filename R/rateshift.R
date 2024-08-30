## find the temporal position of a rate shift using ML
## written by Liam J. Revell 2013, 2014, 2015, 2020, 2023, 2024

rateshift<-function(tree,x,nrates=1,niter=10,method="ML",...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-8
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-FALSE
	if(hasArg(print)) print<-list(...)$print
	else print<-FALSE
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(hasArg(minL)) minL<-list(...)$minL
	else minL<--1e12
	if(hasArg(fixed.shift)) fixed.shift<-list(...)$fixed.shift
	else fixed.shift<-FALSE
	if(hasArg(compute.se)) compute.se<-list(...)$compute.se
	else compute.se<-TRUE
	if(hasArg(parallel)) parallel<-list(...)$parallel
	else parallel<-FALSE
	if(niter==1) parallel<-FALSE
	fn<-if(method=="ML") brownie.lite else brownieREML
	if(!fixed.shift[1]){
		if(print){
			if(!parallel) cat("Optimization progress:\n\n")
			if(nrates>1&&!parallel) cat(paste(c("iter",paste("shift",1:(nrates-1),sep=":"),"logL\n"),collapse="\t"))
			else if(!parallel) cat("iter\ts^2(1)\tlogL\n")
		} else if(niter==1) {
			if(!quiet) cat("Optimizing. Please wait.\n\n")
			flush.console()
		} else { 
			if(!quiet&&!parallel) cat("Optimization progress:\n|")
			flush.console()
		}
	} else {
		if(!quiet) cat("Estimating rates conditioned on input shift points...\n\n")
		nrates<-if(fixed.shift[1]!=TRUE) length(fixed.shift)+1 else 1
		if(nrates>2) fixed.shift<-sort(fixed.shift)
		names(fixed.shift)<-NULL
	}
	lik<-function(par,tree,y,nrates,plot,print,iter,Tol,maxh,minL){
		shift<-sort(c(setNames(0,1),setNames(par,2:nrates)))
		if((any(shift[2:length(shift)]<=0)||any(shift>=maxh))) logL<-minL
		else {
			tree<-make.era.map(tree,shift,tol=Tol/10)
			if(plot){ 
				plotSimmap(tree,setNames(rainbow(nrates),1:nrates),lwd=3,ftype="off",mar=c(0.1,0.1,4.1,0.1))
				title(main=paste("Optimizing rate shift(s), round",iter,"....",sep=" "))
				for(i in 2:(length(shift))) lines(rep(shift[i],2),c(0,length(tree$tip.label)+1),lty="dashed")
			}
			logL<-fn(tree,y)$logL.multiple
		}
		if(print){
			if(nrates>1) cat(paste(c(iter,round(shift[2:length(shift)],4),round(logL,4),"\n"),collapse="\t"))
			else cat(paste(c(iter,round(par,4),round(logL,4),"\n"),collapse="\t"))
		}
		-logL
	}
	h<-max(nodeHeights(tree))
	N<-length(tree$tip.label)
	x<-x[tree$tip.label]
	if(!fixed.shift[1]){
		fit<-list()
		if(!parallel){
			for(i in 1:niter){
				if(nrates==1){ 
					fit[[i]]<-fn(make.era.map(tree,setNames(0,1)),x)
					fit[[i]]$convergence<-if(fit[[i]]$convergence==
						"Optimization has converged.") 0 else 1
				} else {
					fit[[i]]<-list()
					class(fit[[i]])<-"try-error"
					while(inherits(fit[[i]],"try-error")){
						par<-sort(runif(n=nrates-1)*h)
						suppressWarnings(fit[[i]]<-try(optim(par,lik,tree=tree,
							y=x,nrates=nrates,print=print,plot=plot,iter=i,
							Tol=tol,maxh=h,minL=minL)))
					}
				}
				if(!print&&niter>1){ 
					if(!quiet) cat(".")
					flush.console()
				}
			} 
		} else {
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
			fit<-foreach(i=1:niter)%dopar%{
				if(nrates>1) par<-sort(runif(n=nrates-1)*h)
				if(nrates==1){ 
					obj<-fn(make.era.map(tree,setNames(0,1)),x)
					obj$convergence<-if(obj$convergence=="Optimization has converged.") 0 else 1
				} else {
					obj<-list()
					class(obj)<-"try-error"
					while(inherits(obj,"try-error")){
						par<-sort(runif(n=nrates-1)*h)
						suppressWarnings(obj<-try(optim(par,lik,tree=tree,y=x,nrates=nrates,
							print=FALSE,plot=FALSE,iter=i,Tol=tol,maxh=h,minL=minL)))
					}
				}
				obj
			}
			stopCluster(cl=mc)
		}
		if(!print&&niter>1) if(!quiet) if(!parallel) cat("|\nDone.\n\n") else cat("Done.\n\n")
		ll<-sapply(fit,function(x) if(nrates>1) x$value else -x$logL1)
		fit<-fit[[which(ll==min(ll))[1]]]
		frequency.best<-mean(ll<=(min(ll)+1e-4))
		likHess<-if(method=="ML") function(par,tree,y,nrates,tol,maxh){
			sig2<-par[1:nrates]
			shift<-if(nrates>1) setNames(c(0,par[1:(nrates-1)+nrates]),1:nrates) else shift<-setNames(0,1)
			tree<-make.era.map(tree,shift,tol=tol/10)
			mC<-multiC(tree)
			mC<-mapply("*",mC,sig2,SIMPLIFY=FALSE)
			V<-Reduce("+",mC)
			invV<-solve(V)
			a<-as.numeric(colSums(invV)%*%x/sum(invV))
			logL<-sum(dmnorm(y,rep(a,length(x)),V,log=TRUE))
			-logL
		} else if(method=="REML") function(par,tree,y,nrates,tol,maxh){
			sig2<-par[1:nrates]
			shift<-if(nrates>1) setNames(c(0,par[1:(nrates-1)+nrates]),1:nrates) else shift<-setNames(0,1)
			tree<-make.era.map(tree,shift,tol=tol/10)
			tree<-scaleByMap(tree,setNames(sig2,1:nrates))
			picX<-pic(y,tree,scaled=FALSE,var.contrasts=TRUE)
			logL<-sum(dnorm(picX[,1],sd=sqrt(picX[,2]),log=TRUE))
			-logL
		}
		mtree<-if(nrates>1) make.era.map(tree,c(0,fit$par)) else make.era.map(tree,0)
		obj<-fn(mtree,x)
		if(compute.se){
			H<-optimHess(c(obj$sig2.multiple,fit$par),likHess,tree=tree,y=x,nrates=nrates,tol=tol,maxh=h)
			vcv<-if(nrates>1) solve(H) else 1/H
		} else vcv<-matrix(-1,2*nrates-1,2*nrates-1)
		if(nrates>1)
			rownames(vcv)<-colnames(vcv)<-c(paste("sig2(",1:nrates,")",sep=""),paste(1:(nrates-1),"<->",2:nrates,sep=""))
		else rownames(vcv)<-colnames(vcv)<-"sig2(1)"
		obj<-list(sig2=setNames(obj$sig2.multiple,1:nrates),
			shift=if(nrates>1) setNames(fit$par,paste(1:(nrates-1),"<->",2:nrates,sep="")) else NULL,
			vcv=vcv,tree=mtree,logL=obj$logL.multiple,convergence=fit$convergence,message=fit$message,
			method=method,frequency.best=frequency.best)
	} else {
		mtree<-if(nrates>1) make.era.map(tree,c(0,fixed.shift)) else make.era.map(tree,0)
		fit<-fn(mtree,x)
		if(fit$convergence=="Optimization has converged.") fit$convergence<-0
		obj<-list(sig2=setNames(fit$sig2.multiple,1:nrates),
			shift=if(nrates>1) setNames(fixed.shift,paste(1:(nrates-1),"<->",2:nrates,sep="")) else NULL,
			vcv=matrix(-1,2*nrates-1,2*nrates-1),tree=mtree,logL=fit$logL.multiple,convergence=fit$convergence,
			method=method,message="Fitted rates from a fixed shifts",frequency.best=NA)
	}
	class(obj)<-"rateshift"
	if(plot) plot(obj,ftype="off")	
	obj
}

## S3 print method for object of class "rateshift"
## written by Liam J. Revell 2013, 2020
print.rateshift<-function(x,...){
	sqroot<-function(x){
		if(length(x)==1) if(x>=0) sqrt(x) else NaN
		else sapply(x,sqroot)
	}
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	x<-lapply(x,function(a,b) if(is.numeric(a)) round(a,b) else a,b=digits)
	cat(paste("ML ",length(x$sig2),"-rate model:\n",sep=""))
	tmp<-list()
	nn<-vector()
	for(i in 1:length(x$sig2)){
		tmp[2*i-1]<-x$sig2[i]
		tmp[2*i]<-sqroot(diag(x$vcv)[i])
		nn[2*i-1]<-paste("s^2(",names(x$sig2)[i],")",sep="")
		nn[2*i]<-paste("se(",names(x$sig2)[i],")",sep="")
	}
	tmp[2*i+1]<-2*length(x$sig2)
	tmp[2*i+2]<-x$logL
	nn[2*i+1]<-"k"
	nn[2*i+2]<-"logL"
	tmp<-as.data.frame(tmp)
	colnames(tmp)<-nn
	rownames(tmp)<-"value"
	print(tmp,digits=digits)
	if(!is.null(x$shift)){
		cat("\nShift point(s) between regimes (height above root):\n")
		tmp<-list()
		nn<-vector()
		for(i in 1:length(x$shift)){
			tmp[2*i-1]<-x$shift[i]
			tmp[2*i]<-sqroot(diag(x$vcv)[i+length(x$sig2)])
			nn[2*i-1]<-paste(strsplit(names(x$shift[i]),"<->")[[1]],collapse="|")
			nn[2*i]<-paste("se(",paste(strsplit(names(x$shift[i]),"<->")[[1]],
				collapse="|"),")",sep="")
		}
		tmp<-as.data.frame(tmp)
		colnames(tmp)<-nn
		rownames(tmp)<-"value"
		print(tmp,digits=digits)
	} else cat("\nThis is a one-rate model.\n")
	if(x$method=="ML") cat("\nModel fit using ML.\n\n") 
	else if(x$method=="REML") cat("\nModel fit using REML.\n\n")
	cat(paste("Frequency of best fit:",x$frequency.best,"\n\n"))
	if (x$convergence==0) cat(paste("R thinks it has found the",x$method,"solution.\n\n"))
    	else cat("Optimization may not have converged.\n\n")
}

## S3 logLik method for object of class "rateshift"
## written by Liam J. Revell 2013, 2020
logLik.rateshift<-function(object,...){
	logLik<-object$logL
	class(logLik)<-"logLik"
	attr(logLik,"df")<-2*length(object$sig2)
	logLik
}

## S3 plot method for object of class "rateshift"
## written by Liam J. Revell 2015, 2020, 2021, 2024
plot.rateshift<-function(x,...){
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-c(-0.1*Ntip(x$tree),Ntip(x$tree))
	if(hasArg(legend)) legend<-list(...)$legend
	else legend<-TRUE
	if(hasArg(lims)) lims<-list(...)$lims
	else lims<-range(x$sig2)
	if(length(x$sig2)>1||(lims[1]!=lims[2])){
		if(hasArg(col)) col<-list(...)$col
		else col<-c("blue","purple","red")
		cols<-colorRampPalette(col)(101)
		names(cols)<-seq(lims[1],lims[2],by=diff(lims)/100)
		ii<-sapply(x$sig2,function(x,y) order(abs(y-x))[1],
			y=as.numeric(names(cols)))
		colors<-setNames(cols[ii],names(ii))
		args<-list(x=x$tree,ylim=ylim,colors=colors,...)
		args$col<-NULL
		args$lims<-NULL
		args$legend<-NULL
		do.call(plot,args)
		nulo<-lapply(x$shift,function(x,y) lines(rep(x,2),c(1,Ntip(y)),
			lty="dotted",col="grey"),y=x$tree)
		if(legend){
			add.color.bar(leg=0.5*max(nodeHeights(x$tree)),cols=cols,
				prompt=FALSE,x=0,y=-0.05*Ntip(x$tree),lims=round(lims,3),
				title=expression(sigma^2))
		}
	} else {
		if(hasArg(col)) col<-list(...)$col
		else col<-"blue"
		colors<-setNames(col[1],1)
		args<-list(x=x$tree,ylim=ylim,colors=colors,...)
		args$col<-NULL
		args$lims<-NULL
		args$legend<-NULL
		do.call(plot,args)
		if(legend){
			legend(x=0,y=0,
				legend=bquote(sigma^2 == .(round(x$sig2,3))),
				pch=15,col=colors,pt.cex=2,bty="n")
		}
	}
}

## function to visualize the likelihood surface for 2 & 3 rate models (1 & 2 rate-shifts)
## written by Liam J. Revell 2016

likSurface.rateshift<-function(tree,x,nrates=2,shift.range=NULL,
	density=20,plot=TRUE,...){
	h<-max(nodeHeights(tree))
	if(is.null(shift.range)) shift.range<-c(0.01*h,0.99*h)
	shift<-seq(shift.range[1],shift.range[2],length.out=density)
	if(nrates==2){
		cat("Computing likelihood surface for 2-rate (1 rate-shift) model....\n")
		flush.console()
		logL1<-sapply(shift,function(s,tree,x) logLik(rateshift(tree,x,fixed.shift=s,
			quiet=TRUE)),tree=tree,x=x)
		if(plot) plot(shift,logL1,type="l",lwd=2,xlab="shift point",ylab="log(L)",...)
		cat("Done.\n")
		obj<-list(shift=shift,logL=logL1)
	} else if(nrates==3){
		cat("Computing likelihood surface for 3-rate (2 rate-shift) model....\n")
		flush.console()
		logL2<-sapply(shift,function(s1,s2,tree,x) 
			sapply(s2,function(s2,s1,tree,x) 
			logLik(rateshift(tree,x,fixed.shift=if(s1!=s2) c(s1,s2) else s1,
			quiet=TRUE)),s1=s1,tree=tree,x=x),s2=shift,tree=tree,x=x)
		if(plot) contour(shift,shift,logL2,nlevels=20,xlab="shift 1 (or 2)",
			ylab="shift 2 (or 1)",...)
		cat("Done.\n")
		obj<-list(shift=shift,logL=logL2)
	} else if((nrates%in%c(2,3))==FALSE){
		cat("Method only available for nrates==2 and nrates==3\n")
		obj<-NULL
	}
	invisible(obj)
}
