## functions to do the pruning algorithm and reconstruct ancestral states

ancr<-function(object,...) UseMethod("ancr")

ancr.default<-function(object,...){
	warning(paste(
		"ancr does not know how to handle objects of class ",
		class(object),".\n"))
}

logLik.ancr<-function(object,...){
	lik<-object$logLik
	attr(lik,"df")<-attr(object$logLik,"df")
	lik
}

print.ancr<-function(x,digits=6,printlen=6,...){
	cat("Marginal ancestral state estimates:\n")
	if (is.null(printlen)) 
		print(round(x$ace,digits))
	else {
		print(round(x$ace[1:printlen,],digits))
		cat("...\n")
	}
	cat(paste("\nLog-likelihood =",round(logLik(x),
		digits),"\n\n"))
}

## marginal ancestral states for "anova.fitMk" object
ancr.anova.fitMk<-function(object,...){
	if(hasArg(weighted)) weighted<-list(...)$weighted
	else weighted<-TRUE
	if(weighted){
		w<-object$weight
		fits<-attr(object,"models")
		anc<-lapply(fits,function(x,...) ancr(x,...)$ace,...)
		ss<-sort(unique(unlist(lapply(anc,colnames))))
		if(any(sapply(attr(object,"models"),class)=="fitHRM")){
				for(i in 1:length(anc)){
						tmp<-anc[[i]]
						anc[[i]]<-matrix(0,nrow(tmp),length(ss),
							dimnames=list(rownames(tmp),ss))
						anc[[i]][rownames(tmp),colnames(tmp)]<-tmp
				}
		}
		anc<-mapply("*",anc,w,SIMPLIFY=FALSE)
		anc<-Reduce("+",anc)
		foo<-function(obj) unclass(as.Qmatrix(obj))
		Q<-lapply(fits,foo)
		if(any(sapply(attr(object,"models"),class)=="fitHRM")){
			for(i in 1:length(Q)){
				tmp<-Q[[i]]
				Q[[i]]<-matrix(0,length(ss),length(ss),
					dimnames=list(ss,ss))
				Q[[i]][rownames(tmp),colnames(tmp)]<-tmp
			}
		}
		Q<-mapply("*",Q,w,SIMPLIFY=FALSE)
		Q<-Reduce("+",Q)
		model<-matrix(0,nrow(Q),ncol(Q))
		k<-nrow(Q)*(nrow(Q)-1)
		model[col(model)!=row(model)]<-1:k
		q<-sapply(1:k,function(i,Q,model) Q[which(model==i)],
			Q=Q,model=model)
		TREE<-fits[[1]]$tree
		DATA<-fits[[2]]$data
		if(any(sapply(attr(object,"models"),class)=="fitHRM")){
			levs<-unique(gsub("*","",ss,fixed=TRUE))
			tmp<-DATA[,levs]
			DATA<-matrix(0,nrow(tmp),length(ss),
				dimnames=list(rownames(tmp),ss))
			for(i in 1:nrow(tmp)){
				for(j in 1:length(levs)){
					DATA[i,grep(levs[j],ss)]<-tmp[i,levs[j]]
				}
			}
		}
		log_lik<-pruning(q,TREE,DATA,model=model)
		attr(log_lik,"df")<-max(model)
		obj<-list(ace=anc,logLik=log_lik)
		class(obj)<-"ancr"
		return(obj)
	} else {
		best<-which(object$AIC==min(object$AIC))
		return(ancr(object$fits[[best]],...))
	}
}

## marginal ancestral states for "fitHRM" object
ancr.fitHRM<-function(object,...) ancr.fitMk(object,...)

hide.hidden<-function(object,...) UseMethod("hide.hidden")

hide.hidden.default<-function(object,...){
	warning(paste(
		"hide.hidden does not know how to handle objects of class ",
		class(object),".\n"))
}

hide.hidden.ancr<-function(object,...){
	ss<-colnames(object$ace)
	ii<-grep("*",ss,fixed=TRUE)
	if(length(ii)>0) ss<-ss[-ii]
	anc<-matrix(0,nrow(object$ace),length(ss),
		dimnames=list(rownames(object$ace),ss))
	for(i in 1:length(ss)){
		anc[,ss[i]]<-rowSums(object$ace[,grep(ss[i],
			colnames(object$ace)),drop=FALSE])
	}
	anc
}

## marginal ancestral states for "fitpolyMk" object
ancr.fitpolyMk<-function(object,...) ancr.fitMk(object,...)

## marginal ancestral states for "fitMk" object
ancr.fitMk<-function(object,...){
	x<-object$data
	tree<-object$tree
	q<-object$rates
	model<-object$index.matrix
	model[is.na(model)]<-0
	pi=object$pi
	plik<-pruning(q,tree,x,model=model,pi=pi,
		return="conditional")
	if(hasArg(tips)) tips<-list(...)$tips
	else tips<-FALSE
	ace<-marginal_asr(q,tree,plik,model,tips)
	result<-list(ace=ace,
		logLik=pruning(q,tree,x,model=model,pi=pi))
	attr(result$logLik,"df")<-max(model)
	class(result)<-"ancr"
	result
}

## marginal ancestral states for "fitMk" object
ancr.fitMk<-function(object,...){
	x<-object$data
	tree<-object$tree
	q<-object$rates
	model<-object$index.matrix
	model[is.na(model)]<-0
	pi=object$pi
	plik<-pruning(q,tree,x,model=model,pi=pi,
		return="conditional")
	if(hasArg(tips)) tips<-list(...)$tips
	else tips<-FALSE
	ace<-marginal_asr(q,tree,plik,model,tips)
	result<-list(ace=ace,
		logLik=pruning(q,tree,x,model=model,pi=pi))
	attr(result$logLik,"df")<-max(model)
	class(result)<-"ancr"
	result
}

pruning<-function(q,tree,x,model=NULL,...){
	if(hasArg(return)) return<-list(...)$return
	else return<-"likelihood"
	pw<-if(!is.null(attr(tree,"order"))&&
		attr(tree,"order")=="postorder") tree else 
		reorder(tree,"postorder")
	k<-ncol(x)
	if(is.null(model)){
		model<-matrix(1,k,k)
		diag(model)<-0
	}
	if(hasArg(pi)) pi<-list(...)$pi
	else pi<-rep(1/k,k)
	Q<-matrix(0,k,k)
	Q[]<-c(0,q)[model+1]
	diag(Q)<--rowSums(Q)
	L<-rbind(x[pw$tip.label,],
		matrix(0,tree$Nnode,k,
		dimnames=list(1:tree$Nnode+Ntip(tree))))
	nn<-unique(pw$edge[,1])
	pp<-vector(mode="numeric",length=length(nn))
	root<-min(nn)
	for(i in 1:length(nn)){
		ee<-which(pw$edge[,1]==nn[i])
		PP<-matrix(NA,length(ee),k)
		for(j in 1:length(ee)){
			P<-expm(Q*pw$edge.length[ee[j]])
			PP[j,]<-P%*%L[pw$edge[ee[j],2],]
		}
		L[nn[i],]<-apply(PP,2,prod)
		if(nn[i]==root){
			if(pi[1]=="fitzjohn") pi<-L[nn[i],]/sum(L[nn[i],])
			L[nn[i],]<-pi*L[nn[i],]
		}
		pp[i]<-sum(L[nn[i],])
		L[nn[i],]<-L[nn[i],]/pp[i]
	}
	prob<-sum(log(pp))
	if(return=="likelihood") 
		if(is.na(prob)||is.nan(prob)) 
			return(-Inf) else return(prob)
	else if(return=="conditional") L
	else if(return=="pi") pi
}

marginal_asr<-function(q,tree,L,model=NULL,tips=FALSE){
	pw<-reorder(tree,"postorder")
	k<-ncol(L)
	if(is.null(model)){
		model<-matrix(1,k,k)
		diag(model)<-0
	}
	Q<-matrix(0,k,k)
	Q[]<-c(0,q)[model+1]
	diag(Q)<--rowSums(Q)
	nn<-unique(pw$edge[,1])
	for(i in length(nn):1){
		ee<-which(pw$edge[,1]==nn[i])
		for(j in 1:length(ee)){
			P<-expm(Q*pw$edge.length[ee[j]])
			pp<-t(L[nn[i],]/(P%*%L[pw$edge[ee[j],2],]))
			pp[is.nan(pp)]<-0
			L[pw$edge[ee[j],2],]<-(pp%*%P)*
				L[pw$edge[ee[j],2],]
		}
	}
	anc<-L/rep(rowSums(L),k)
	if(tips) anc else anc[1:Nnode(tree)+Ntip(tree),]
}
