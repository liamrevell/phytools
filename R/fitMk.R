## function for conditional likelihoods at nodes
## written by Liam J. Revell 2015
## with input from (& structural similarity to) function ace by E. Paradis et al. 2013

fitMk<-function(tree,x,model,fixedQ=NULL,...){
	if(hasArg(output.liks)) output.liks<-list(...)$output.liks
	else output.liks<-FALSE
	N<-Ntip(tree)
	M<-tree$Nnode
	if(is.matrix(x)){
		x<-x[tree$tip.label,]
		m<-ncol(x)
		states<-colnames(x)
	} else {
		x<-to.matrix(x,sort(unique(x)))
		m<-ncol(x)
		states<-colnames(x)
	}
	if(hasArg(pi)) pi<-list(...)$pi
	else pi<-"equal"
	if(is.character(pi)&&pi=="equal") pi<-setNames(rep(1/m,m),states)
	if(is.null(fixedQ)){
		if(is.character(model)){
			rate<-matrix(NA,m,m)
			if(model=="ER") k<-rate[]<-1
			else if(model=="ARD"){
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
		Q<-matrix(0,m,m)
	} else {
		rate<-matrix(NA,m,m)
		k<-m*(m-1)
		rate[col(rate)!=row(rate)]<-1:k
		Q<-fixedQ
	}
	index.matrix<-rate
	tmp<-cbind(1:m,1:m)
	rate[tmp]<-0
	rate[rate==0]<-k+1
	liks<-rbind(x,matrix(0,M,m,dimnames=list(1:M+N,states)))
	pw<-reorder(tree,"pruningwise")
	lik<-function(pp,output.liks=FALSE,pi){
		if(any(is.nan(pp))||any(is.infinite(pp))) return(1e50)
		comp<-vector(length=N+M,mode="numeric")
		Q[]<-c(pp,0)[rate]
		diag(Q)<--rowSums(Q)
		parents<-unique(pw$edge[,1])
		root<-min(parents)
		for(i in 1:length(parents)){
			anc<-parents[i]
			ii<-which(pw$edge[,1]==parents[i])
			desc<-pw$edge[ii,2]
			el<-pw$edge.length[ii]
			v<-vector(length=length(desc),mode="list")
			for(j in 1:length(v))
				v[[j]]<-expm(Q*el[j])%*%liks[desc[j],]
			vv<-if(anc==root) Reduce('*',v)[,1]*pi else Reduce('*',v)[,1]
			comp[anc]<-sum(vv)
			liks[anc,]<-vv/comp[anc]
		}
		if(output.liks) return(liks[1:M+N,])
		logL<--2*sum(log(comp[1:M+N]))
		return(if(is.na(logL)) Inf else logL)
	}
	if(is.null(fixedQ)){
		fit<-nlminb(rep(0.1,k),function(p) lik(p,pi=pi),lower=rep(0,k),upper=rep(1e50,k))
		obj<-list(logLik=-fit$objective/2,
			rates=fit$par,
			index.matrix=index.matrix,
			states=states)
		if(output.liks) obj$lik.anc<-lik(obj$rates,TRUE,pi=pi)
	} else {
		fit<-lik(fixedQ[sapply(1:k,function(x,y) which(x==y),index.matrix)],pi=pi)
		obj<-list(logLik=fit/2,
			rates=fixedQ[sapply(1:k,function(x,y) which(x==y),index.matrix)],
			index.matrix=index.matrix,
			states=states)
		if(output.liks) obj$lik.anc<-lik(obj$rates,TRUE,pi=pi)
	}
	class(obj)<-"fitMk"
	obj
}

print.fitMk<-function(x,digits=6,...){
	cat("Object of class \"fitMk\".\n\n")
	cat("Fitted (or set) value of Q:\n")
	Q<-matrix(NA,length(obj$states),length(obj$states))
	Q[]<-obj$rates[obj$index.matrix]
	diag(Q)<-0
	diag(Q)<--rowSums(Q)
	colnames(Q)<-rownames(Q)<-obj$states
	print(round(Q,digits))
	cat(paste("\nLog-likelihood:",round(obj$logLik,digits),"\n\n"))
} 
