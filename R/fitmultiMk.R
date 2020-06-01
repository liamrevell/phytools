## new function to fit multi-regime Mk model
## written by Liam J. Revell 2017
## adapted from ape::ace by E. Paradis et al. 2013

fitmultiMk<-function(tree,x,model="ER",...){
	if(!inherits(tree,"simmap")){
		stop("tree should be an object of class \"simmap\". Use fitMk.\n")
	} else {
		regimes<-mapped.states(tree)
		nregimes<-length(regimes)
	}
	if(hasArg(q.init)) q.init<-list(...)$q.init
	else q.init<-length(unique(x))/sum(tree$edge.length)
	if(hasArg(opt.method)) opt.method<-list(...)$opt.method
	else opt.method<-"nlminb"
	if(hasArg(min.q)) min.q<-list(...)$min.q
	else min.q<-1e-12
	if(is.matrix(x)){
		x<-x[tree$tip.label,]
		m<-ncol(x)
		states<-colnames(x)
	} else {
		x<-to.matrix(x,sort(unique(x)))
		x<-x[tree$tip.label,]
		m<-ncol(x)
		states<-colnames(x)
	}
	if(hasArg(pi)) pi<-list(...)$pi
	else pi<-"equal"
	if(pi[1]=="equal") pi<-setNames(rep(1/m,m),states)	
	else pi<-pi/sum(pi)
	if(is.character(model)){
		rate<-matrix(NA,m,m)
		if(model=="ER"){ 
			k<-rate[]<-1
			diag(rate)<-NA
		} else if(model=="ARD"){
			k<-m*(m-1)
			rate[col(rate)!=row(rate)]<-1:k
		} else if(model=="SYM"){
			k<-m*(m-1)/2
			ii<-col(rate)<row(rate)
			rate[ii]<-1:k
			rate<-t(rate)
			rate[ii]<-1:k
		}
	} else {
		if(ncol(model)!=nrow(model)) 
			stop("model is not a square matrix")
		if(ncol(model)!=ncol(x)) 
			stop("model does not have the right number of columns")
		rate<-model
		k<-max(rate)
	}
	Q<-replicate(nregimes,matrix(0,m,m),simplify=FALSE)
	names(Q)<-regimes
	index.matrix<-rate
	tmp<-cbind(1:m,1:m)
	rate[tmp]<-0
	rate[rate==0]<-k+1
	pw<-reorder(map.to.singleton(tree),"pruningwise")
	N<-Ntip(pw)
	M<-pw$Nnode
	liks<-rbind(x,matrix(0,M,m,dimnames=list(1:M+N,states)))
	lik<-function(pp,output.liks=FALSE,pi){
		if(any(is.nan(pp))||any(is.infinite(pp))) return(1e50)
		comp<-vector(length=N+M,mode="numeric")
		for(i in 1:nregimes){ 
			Q[[i]][]<-c(pp[1:k+(i-1)*k],0)[rate]
			diag(Q[[i]])<--rowSums(Q[[i]])
		}
		parents<-unique(pw$edge[,1])
		root<-min(parents)
		for(i in 1:length(parents)){
			anc<-parents[i]
			ii<-which(pw$edge[,1]==parents[i])
			desc<-pw$edge[ii,2]
			el<-pw$edge.length[ii]
			v<-vector(length=length(desc),mode="list")
			reg<-names(pw$edge.length)[ii]
			for(j in 1:length(v))
				v[[j]]<-EXPM(Q[[reg[j]]]*el[j])%*%liks[desc[j],]
			vv<-if(anc==root) Reduce('*',v)[,1]*pi else 
				Reduce('*',v)[,1]
			comp[anc]<-sum(vv)
			liks[anc,]<-vv/comp[anc]
		}
		logL<--sum(log(comp[1:M+N]))
		return(if(is.na(logL)) Inf else logL)
	}
	if(length(q.init)!=(k*nregimes)) q.init<-rep(q.init[1],k*nregimes)
	if(opt.method=="optim")
		fit<-optim(q.init,function(p) lik(p,pi=pi),method="L-BFGS-B",
			lower=rep(min.q,k))
	else	
		fit<-nlminb(q.init,function(p) lik(p,pi=pi),
			lower=rep(0,k*nregimes),upper=rep(1e50,k*nregimes))
	obj<-list(logLik=
		if(opt.method=="optim") -fit$value else -fit$objective,
		rates=fit$par,
		index.matrix=index.matrix,
		states=states,
		regimes=regimes,
		pi=pi,
		method=opt.method)
	class(obj)<-"fitmultiMk"
	return(obj)
}

## print method for objects of class "fitMk"
print.fitmultiMk<-function(x,digits=6,...){
	k<-max(x$index.matrix,na.rm=TRUE)
	cat("Object of class \"fitmultiMk\".\n\n")
	for(i in 1:length(x$regimes)){
		cat(paste("Fitted value of Q[",x$regimes[i],"]:\n",sep=""))
		Q<-matrix(NA,length(x$states),length(x$states))
		Q[]<-c(0,x$rates[1:k+(i-1)*k])[x$index.matrix+1]
		diag(Q)<-0
		diag(Q)<--rowSums(Q)
		colnames(Q)<-rownames(Q)<-x$states
		print(round(Q,digits))
		cat("\n")
	}
	cat("Fitted (or set) value of pi:\n")
	print(x$pi)
	cat(paste("\nLog-likelihood:",round(x$logLik,digits),"\n"))
	cat(paste("\nOptimization method used was \"",x$method,"\"\n\n",sep=""))
}

## summary method for objects of class "fitmultiMk"
summary.fitmultiMk<-function(object,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-6
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	Q<-list()
	for(i in 1:length(object$regimes)){
		if(!quiet) cat(paste("Fitted value of Q[",object$regimes[i],
			"]:\n",sep=""))
		Q[[i]]<-matrix(NA,length(object$states),length(object$states))
		Q[[i]][]<-c(0,object$rates)[object$index.matrix+1]
		diag(Q[[i]])<-0
		diag(Q[[i]])<--rowSums(Q[[i]])
		colnames(Q[[i]])<-rownames(Q[[i]])<-object$states
		if(!quiet) print(round(Q[[i]],digits))
		cat("\n")
	}
	names(Q)<-object$regimes
	if(!quiet) cat(paste("Log-likelihood:",round(object$logLik,digits),"\n\n"))
	invisible(list(Q=Q,logLik=object$logLik))
}

## logLik method for objects of class "fitmultiMk"
logLik.fitmultiMk<-function(object,...){ 
	lik<-object$logLik
	attr(lik,"df")<-length(object$rates)
	lik
}
