## find the temporal position of a rate shift using ML
## written by Liam J. Revell 2013, 2014, 2015

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
	fn<-if(method=="ML") brownie.lite else brownieREML
	if(!fixed.shift[1]){
		if(print){
			cat("Optimization progress:\n\n")
			if(nrates>1) cat(paste(c("iter",paste("shift",1:(nrates-1),sep=":"),"logL\n"),collapse="\t"))
			else cat("iter\ts^2(1)\tlogL\n")
		} else if(niter==1) {
			if(!quiet) cat("Optimizing. Please wait.\n\n")
			flush.console()
		} else { 
			if(!quiet) cat("Optimization progress:\n|")
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
		for(i in 1:niter){
			if(nrates>1) par<-sort(runif(n=nrates-1)*h)
			if(nrates==1){ 
				fit[[i]]<-fn(make.era.map(tree,setNames(0,1)),x)
				fit[[i]]$convergence<-if(fit[[i]]$convergence=="Optimization has converged.") 0 else 1
			} else suppressWarnings(fit[[i]]<-optim(par,lik,tree=tree,y=x,nrates=nrates,print=print,plot=plot,
				iter=i,Tol=tol,maxh=h,minL=minL))
			if(!print&&niter>1){ 
				if(!quiet) cat(".")
				flush.console()
			}
		}
		if(!print&&niter>1) if(!quiet) cat("|\nDone.\n\n")
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
		H<-optimHess(c(obj$sig2.multiple,fit$par),likHess,tree=tree,y=x,nrates=nrates,tol=tol,maxh=h)
		vcv<-if(nrates>1) solve(H) else 1/H
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
## written by Liam J. Revell 2013
print.rateshift<-function(x,...){
	sqroot<-function(x){
		if(length(x)==1) if(x>=0) sqrt(x) else NaN
		else sapply(x,sqroot)
	}
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	x<-lapply(x,function(a,b) if(is.numeric(a)) round(a,b) else a,b=digits)
	cat(paste("ML ",length(x$sig2),"-rate model:\n",sep=""))
	cat(paste(c("",paste("s^2(",names(x$sig2),")","\tse(",names(x$sig2),")",sep=""), 
      	"k","logL","\n"),collapse="\t"))
	cat(paste(paste(c("value",paste(x$sig2,round(sqroot(diag(x$vcv)[1:length(x$sig2)]),digits),
		sep="\t"),2*length(x$sig2),x$logL),collapse="\t"),"\n\n",sep=""))
	if(!is.null(x$shift)){
		cat("Shift point(s) between regimes (height above root):\n")
		nn<-sapply(strsplit(names(x$shift),"<->"),paste,collapse="|")
		cat(paste(c("",paste(nn,paste("se(",nn,")",sep=""),sep="\t"),"\n"),
			collapse="\t"))
		cat(paste(paste(c("value",paste(x$shift,
			round(sqroot(diag(x$vcv)[1:length(x$shift)+length(x$sig2)]),digits),
			sep="\t")),collapse="\t"),"\n\n",sep=""))
	} else cat("This is a one-rate model.\n\n")
	if(x$method=="ML") cat("Model fit using ML.\n\n") 
	else if(x$method=="REML") cat("Model fit using REML.\n\n")
	cat(paste("Frequency of best fit:",x$frequency.best,"\n\n"))
	if (x$convergence==0) cat(paste("R thinks it has found the",x$method,"solution.\n\n"))
    	else cat("Optimization may not have converged.\n\n")
}

## S3 logLik method for object of class "rateshift"
## written by Liam J. Revell 2013
logLik.rateshift<-function(object,...){
	logLik<-object$logL
	class(logLik)<-"logLik"
	attr(logLik,"df")<-2*length(object$sig2)
	logLik
}

## S3 plot method for object of class "rateshift"
## written by Liam J. Revell 2015
plot.rateshift<-function(x,...){
	if(length(x$sig2)>1){
		cols<-colorRampPalette(c("blue","purple","red"))(101)
		rr<-range(x$sig2)
		names(cols)<-seq(rr[1],rr[2],by=diff(rr)/100)
		ii<-sapply(x$sig2,function(x,y) order(abs(y-x))[1],
			y=as.numeric(names(cols)))
		colors<-setNames(cols[ii],names(ii))
		plot(x$tree,ylim=c(-0.1*Ntip(x$tree),Ntip(x$tree)),
			colors=colors,...)
		nulo<-lapply(x$shift,function(x,y) lines(rep(x,2),c(1,Ntip(y)),
			lty="dotted",col="grey"),y=x$tree)
		add.color.bar(leg=0.5*max(nodeHeights(x$tree)),cols=cols,
			prompt=FALSE,x=0,y=-0.05*Ntip(x$tree),lims=round(rr,3),
			title=expression(sigma^2))
	} else {
		colors<-setNames("blue",1)
		plot(x$tree,ylim=c(-0.1*Ntip(x$tree),Ntip(x$tree)),
			colors=colors,...)
		txt<-as.character(round(x$sig2,3))
		add.simmap.legend(leg=expression(paste(sigma^2," = ",sep="")),
			colors="blue",prompt=FALSE,x=0,y=-0.05*Ntip(x$tree))
		text(x=5.5*strwidth("W"),y=-0.05*Ntip(x$tree),round(x$sig2,3))
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
