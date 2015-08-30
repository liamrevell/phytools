## find the temporal position of a rate shift using ML
## written by Liam J. Revell 2013, 2014, 2015

rateshift<-function(tree,x,nrates=1,niter=10,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-8
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-FALSE
	if(hasArg(print)) print<-list(...)$print
	else print<-FALSE
	if(hasArg(minL)) minL<-list(...)$minL
	else minL<--1e12
	if(print){
		cat("Optimization progress:\n\n")
		if(nrates>1)
			cat(paste(c("iter",paste("s^2(",1:nrates,")",sep=""),
				paste("shift",1:(nrates-1),sep=":"),"logL\n"),collapse="\t"))
		else cat("iter\ts^2(1)\tlogL\n")
	} else if(niter==1) {
		cat("Optimizing. Please wait.\n\n")
		flush.console()
	} else { 
		cat("Optimization progress:\n|")
		flush.console()
	}
	lik<-function(par,tree,y,nrates,plot,print,iter,tol,maxh,for.hessian=FALSE,minL){
		sig2<-par[1:nrates]
		shift<-if(nrates>1) setNames(sort(c(0,par[1:(nrates-1)+nrates])),1:nrates) else shift<-setNames(0,1)
		if((nrates>1)&&(any(sig2<=0)||any(shift[2:length(shift)]<=0)||any(shift>=maxh))&&(for.hessian==FALSE)) logL<-minL
		else {
			tree<-make.era.map(tree,shift,tol=tol/10)
			if(plot&&!for.hessian){ 
				plotSimmap(tree,setNames(rainbow(nrates),1:nrates),lwd=3,ftype="off",mar=c(0.1,0.1,4.1,0.1))
				title(main=paste("Optimizing rate shift(s), round",iter,"....",sep=" "))
				for(i in 2:(length(shift))) lines(rep(shift[i],2),c(0,length(tree$tip.label)+1),lty="dashed")
			}
			mC<-multiC(tree)
			mC<-mapply("*",mC,sig2,SIMPLIFY=FALSE)
			V<-Reduce("+",mC)
			invV<-solve(V)
			a<-as.numeric(colSums(invV)%*%x/sum(invV))
			logL<-sum(dmnorm(y,rep(a,length(x)),V,log=TRUE))
		}
		if(print&&!for.hessian){
			if(nrates>1) cat(paste(c(iter,round(sig2,4),round(shift[2:length(shift)],4),round(logL,4),"\n"),collapse="\t"))
			else cat(paste(c(iter,round(sig2,4),round(logL,4),"\n"),collapse="\t"))
		}
		-logL
	}
	h<-max(nodeHeights(tree))
	N<-length(tree$tip.label)
	x<-x[tree$tip.label]
	fit<-list()
	for(i in 1:niter){
		if(nrates>1)
			par<-c(rchisq(n=nrates,df=N)/N*rep(phyl.vcv(as.matrix(x),vcv(tree),lambda=1)$R[1,1],nrates),
				sort(runif(n=nrates-1))*h)
		else par<-phyl.vcv(as.matrix(x),vcv(tree),lambda=1)$R[1,1]
		fit[[i]]<-if(nrates==1) optim(par,lik,tree=tree,y=x,nrates=nrates,print=print,plot=plot,iter=i,tol=tol,maxh=h,minL=minL,
				method="Brent",lower=tol,upper=10*phyl.vcv(as.matrix(x),vcv(tree),lambda=1)$R[1,1])
			else optim(par,lik,tree=tree,y=x,nrates=nrates,print=print,plot=plot,iter=i,tol=tol,maxh=h,minL=minL)
		if(!print&&niter>1){ 
			cat(".")
			flush.console()
		}
	}
	if(!print&&niter>1) cat("|\nDone.\n\n")
	ll<-sapply(fit,function(x) x$value)
	fit<-fit[[which(ll==min(ll))[1]]]
	H<-optimHess(fit$par,lik,tree=tree,y=x,nrates=nrates,print=print,plot=plot,iter=NA,tol=tol,maxh=h,for.hessian=TRUE,minL=minL)
	vcv<-if(nrates>1) solve(H) else 1/H
	if(nrates>1)
		rownames(vcv)<-colnames(vcv)<-c(paste("sig2(",1:nrates,")",sep=""),paste(1:(nrates-1),"<->",2:nrates,sep=""))
	else rownames(vcv)<-colnames(vcv)<-"sig2(1)"
	obj<-list(sig2=setNames(fit$par[1:nrates],1:nrates),
		shift=if(nrates>1)setNames(fit$par[1:(nrates-1)+nrates],paste(1:(nrates-1),"<->",2:nrates,sep="")) else NULL,
		vcv=vcv,
		tree=if(nrates>1) make.era.map(tree,c(0,fit$par[1:(nrates-1)+nrates])) else make.era.map(tree,0),
		logL=-fit$value,convergence=fit$convergence,message=fit$message)
	class(obj)<-"rateshift"
	if(plot){
		plotSimmap(obj$tree,setNames(rainbow(nrates),1:nrates),lwd=3,ftype="off",mar=c(0.1,0.1,4.1,0.1))
		title(main="ML optimized rate shift(s)")
		if(nrates>1)
			for(i in 1:(length(obj$shift))) lines(rep(obj$shift[i],2),c(0,length(obj$tree$tip.label)+1),lty="dashed")
	}	
	return(obj)
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
	if (x$convergence==0) cat("R thinks it has found the ML solution.\n\n")
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
