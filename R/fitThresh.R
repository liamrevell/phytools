## compute trace (used internally)
tr<-function(X) sum(diag(X))

## convert discrete character to liability matrix
thresh2bin<-function(liability,threshold,X){
	Y<-matrix(0,length(liability),nrow(X),
		dimnames=list(round(liability,3),rownames(X)))
	interval<-mean(liability[2:length(liability)]-
			liability[1:(length(liability)-1)])
	down<-liability-interval/2
	## down[1]<--Inf
	up<-liability+interval/2
	## up[length(up)]<-Inf
	xx<-rep(0,length(liability))
	for(i in 1:(length(threshold)-1)){
		xx[]<-0
		ind<-intersect(which(up>threshold[i]),
			which(down<threshold[i+1]))
		xx[ind]<-1
		ind<-setdiff(ind,
			intersect(which(down>threshold[i]),
				which(up<threshold[i+1])))
		if(length(ind)>0){
			for(j in 1:length(ind)){
				xx[ind[j]]<-if(down[ind[j]]>threshold[i]&&
						down[ind[j]]<threshold[i+1]) 
					(threshold[i+1]-down[ind[j]])/interval else
						abs(up[ind[j]]-threshold[i])/interval
			}
		}
		Y[,which(X[,i]==1)]<-xx
	}
	t(Y)
}

## fit threshold model using discrete approximation
fitThresh<-function(tree,x,sequence=NULL,...){
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-0
	if(hasArg(levs)) levs<-list(...)$levs
	else levs<-200
	if(hasArg(root)) root<-list(...)$root
	else root<-"fitzjohn"
	if(hasArg(rand_start)) rand_start<-list(...)$rand_start
	else rand_start<-TRUE
	C<-vcv(tree)
	N<-Ntip(tree)
	Ed<-tr(C)/N-sum(C)/(N^2)
	lims<-qnorm(c(0.005,0.995),sd=sqrt(Ed))
	liability<-seq(lims[1],lims[2],length.out=levs)
	if(!is.factor(x)) x<-setNames(as.factor(x),names(x))
	if(is.null(sequence)) sequence<-levels(x)
	X<-to.matrix(x,sequence)
	q<-1/(2*(diff(lims)/levs)^2)
	model<-matrix(0,levs,levs,
		dimnames=list(round(liability,3),round(liability,3)))
	ind<-cbind(1:(nrow(model)-1),2:nrow(model))
	model[ind]<-1
	model[ind[,2:1]]<-1
	Q<-model*q
	diag(Q)<--rowSums(Q)
	pw<-if(!is.null(attr(tree,"order"))&&
			attr(tree,"order")=="postorder") tree else 
				reorder(tree,"postorder")
	P<-lapply(pw$edge.length,function(e,Q) expm(e*Q),
		Q=Q)
	if(length(sequence)==2){
		Y<-thresh2bin(liability,c(-Inf,0,Inf),X)
		lnL<-lik_thresh(c(),pw,fixed_threshold=0,
			liability,X,P,trace=trace,pi=root)
		attr(lnL,"df")<-1
		threshold<-0
		opt<-list(convergence=0)
	} else if(length(sequence)>2){
		fixed_threshold<-min(liability)+
			diff(range(liability))/length(sequence)
		if(length(sequence)==3){
			opt<-optimize(lik_thresh,
				c(fixed_threshold,max(liability)),
				pw=pw,fixed_threshold=fixed_threshold,
				liability=liability,x=X,P.all=P,trace=trace,
				pi=root,maximum=TRUE)
			lnL<-opt$objective
			attr(lnL,"df")<-2
			threshold<-c(fixed_threshold,opt$maximum)
		} else if(length(sequence)>3){
			init<-seq(fixed_threshold,max(liability),
				length.out=length(sequence))[
					-c(1,length(sequence))]
			if(rand_start) init<-runif(n=length(init))*init
			opt<-optim(init,function(p) -lik_thresh(p,pw=pw,
				fixed_threshold=fixed_threshold,liability=liability,
				x=X,pi=root,P.all=P,trace=trace),method="L-BFGS",
				lower=rep(fixed_threshold,length(init)),
				upper=rep(max(liability),length(init)))
			lnL<--opt$value
			attr(lnL,"df")<-length(sequence)-1
			threshold<-c(fixed_threshold,sort(opt$par))
		}
	}
	Y<-thresh2bin(liability,c(-Inf,threshold,Inf),X)
	fit<-list(
		logLik=lnL,
		rates=q,
		index.matrix=model,
		states=colnames(Y),
		pi=lik_thresh(threshold,pw,c(),liability,X,P,
			trace=0,pi=root,return="pi"),
		opt_results=opt,
		data=Y,
		tree=tree)
	class(fit)<-"fitMk"
	object<-list(sigsq=1.0,bounds=lims,ncat=levs,
		liability=liability,threshold=threshold,
		logLik=lnL,tree=tree,data=X,mk_fit=fit)
	class(object)<-c("fitThresh","fitMk")
	object
}

## print for "fitThresh" object
print.fitThresh<-function(x, digits=6, ...){
	spacer<-if(length(x$threshold)>2) "\n        " else ""
	cat("Object of class \"fitThresh\".\n\n")
	cat("    Set value of sigsq (of the liability) = 1.0\n\n")
	cat(paste("	  Set or estimated threshold(s) =",spacer,"[",
		paste(round(x$threshold,digits),collapse=", "),"]*\n\n"))
	cat(paste("    Log-likelihood:", round(x$logLik, digits),
		"\n\n"))
	cat("(*lowermost threshold is fixed)\n\n")
}

## logLik for "fitThresh" object
logLik.fitThresh<-function(object,...) object$logLik

## marginal ancestral states of "fitThresh" object
ancr.fitThresh<-function(object,...){
	anc_mk<-ancr(object$mk_fit,type="marginal")
	anc<-matrix(0,object$tree$Nnode,ncol(object$data),
		dimnames=list(1:object$tree$Nnode+Ntip(object$tree),
			colnames(object$data)))
	tmp<-diag(rep(1,ncol(object$data)))
	colnames(tmp)<-colnames(object$data)
	xx<-thresh2bin(object$liability,
		c(-Inf,object$threshold,Inf),tmp)
	for(i in 1:nrow(anc)){
		for(j in 1:nrow(xx)){
			anc[i,j]<-sum(anc_mk$ace[i,]*xx[j,])
		}
	}
	result<-list(ace=anc,logLik=logLik(object$mk_fit))
	attr(result,"type")<-"marginal"
	attr(result,"tree")<-object$tree
	attr(result,"data")<-object$data
	class(result)<-"ancr"
	result
}

## likelihood function for thresholds
lik_thresh<-function(threshold,pw,fixed_threshold,
	liability,x,P.all,trace=1,...){
	if(hasArg(return)) return<-list(...)$return
	else return<-"likelihood"
	th<-c(-Inf,fixed_threshold,sort(threshold),Inf)
	Y<-thresh2bin(liability,th,x)
	k<-ncol(Y)
	if(hasArg(pi)) pi<-list(...)$pi
	else pi<-rep(1/k,k)
	L<-rbind(Y[pw$tip.label,],
		matrix(0,pw$Nnode,k,
			dimnames=list(1:pw$Nnode+Ntip(pw))))
	nn<-unique(pw$edge[,1])
	pp<-vector(mode="numeric",length=length(nn))
	root<-min(nn)
	for(i in 1:length(nn)){
		ee<-which(pw$edge[,1]==nn[i])
		PP<-matrix(NA,length(ee),k)
		for(j in 1:length(ee)){
			P<-P.all[[ee[j]]]
			PP[j,]<-P%*%L[pw$edge[ee[j],2],]
		}
		L[nn[i],]<-apply(PP,2,prod)
		if(nn[i]==root){
			if(pi[1]=="fitzjohn") pi<-L[nn[i],]/sum(L[nn[i],])
			else if(pi[1]=="mle") 
				pi<-as.numeric(L[nn[i],]==max(L[nn[i],]))
			L[nn[i],]<-pi*L[nn[i],]
		}
		pp[i]<-sum(L[nn[i],])
		L[nn[i],]<-L[nn[i],]/pp[i]
	}
	prob<-sum(log(pp))
	if(trace>0) print(c(fixed_threshold,sort(threshold),
		prob))
	if(return=="likelihood") 
		if(is.na(prob)||is.nan(prob)) 
			return(-Inf) else return(prob)
	else if(return=="conditional") L
	else if(return=="pi") pi
}
