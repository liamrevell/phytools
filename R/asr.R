## functions to do the pruning algorithm and reconstruct ancestral states

anceb<-function(object,...) UseMethod("anceb")

anceb.default<-function(object,...){
	warning(paste(
		"anceb does not know how to handle objects of class ",
		class(object),".\n"))
}

logLik.anceb<-function(object,...){
	lik<-object$logLik
	attr(lik,"df")<-attr(object$logLik,"df")
	lik
}

print.anceb<-function(x,digits=6,printlen=6,...){
	cat("Marginal ancestral character estimates:\n")
	if (is.null(printlen)) 
		print(round(x$ace,digits))
	else {
		print(round(x$ace[1:printlen,],digits))
		cat("...\n")
	}
	cat(paste("\nLog-likelihood =",round(x$logLik,digits), 
		"\n\n"))
}

## object expected to be of class "fitMk"
anceb.fitMk<-function(object,...){
	x<-object$data
	tree<-object$tree
	q<-object$rates
	model<-object$index.matrix
	model[is.na(model)]<-0
	pi=object$pi
	plik<-pruning(q,tree,x,model=model,pi=pi,
		return="conditional")
	ace<-marginal_asr(q,tree,plik,model=model)
	result<-list(ace=ace,
		logLik=pruning(q,tree,x,model=model,pi=pi))
	attr(result$logLik,"df")<-max(model)
	class(result)<-"anceb"
	result
}

pruning<-function(q,tree,x,model=NULL,...){
	if(hasArg(return)) return<-list(...)$return
	else return<-"likelihood"
	pw<-reorder(tree,"postorder")
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
	L<-rbind(x,matrix(0,tree$Nnode,k,
	dimnames=list(1:tree$Nnode+Ntip(tree))))
	nn<-unique(pw$edge[,1])
	for(i in 1:length(nn)){
		ee<-which(pw$edge[,1]==nn[i])
		PP<-matrix(NA,length(ee),k)
		for(j in 1:length(ee)){
			P<-expm(Q*pw$edge.length[ee[j]])
			PP[j,]<-P%*%L[pw$edge[ee[j],2],]
		}
		L[nn[i],]<-apply(PP,2,prod)
	}
	prob<-log(sum(pi*L[nn[i],]))
	if(return=="likelihood") prob
	else if(return=="conditional") L
}

marginal_asr<-function(q,tree,L,model=NULL){
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
	anc[1:Nnode(tree)+Ntip(tree),]
}
