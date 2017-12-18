mcmcMk<-function(tree,x,model="ER",ngen=10000,...){
	log.prior<-function(x,prior) sum(dexp(x,prior,log=TRUE))
	proposal<-function(q,pv) abs(q+rnorm(n=length(q),sd=sqrt(pv)))
	makeQ<-function(m,q,index.matrix){
		Q<-matrix(0,m,m)
		Q[]<-c(0,q)[index.matrix+1]
		diag(Q)<-0
		diag(Q)<--rowSums(Q)
		Q
	}		
	if(hasArg(q)) q<-list(...)$q
	else q<-if(is.matrix(x)) nrow(x)/sum(tree$edge.length) else
		length(unique(x))/sum(tree$edge.length)
	if(hasArg(prop.var)) prop.var<-list(...)$prop.var
	else prop.var<-if(is.matrix(x)) sqrt(0.01*nrow(x)/sum(tree$edge.length)) else
		sqrt(0.01*length(unique(x))/sum(tree$edge.length))
	if(hasArg(prior.rate)) prior.rate<-list(...)$prior.rate
	else prior.rate<-if(is.matrix(x)) 1/(nrow(x)/sum(tree$edge.length)) else
		1/(length(unique(x))/sum(tree$edge.length))
	if(hasArg(print)) print<-list(...)$print
	else print<-100
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
	index.matrix<-rate
	if(length(q)!=k) q<-rep(q[1],k)
	if(length(prop.var)!=k) prop.var<-rep(prop.var[1],k)
	if(length(prior.rate)!=k) prior.rate<-rep(prior.rate[1],k)
	likQ<-logLik(fitMk(tree,x,model=model,fixedQ=makeQ(m,q,index.matrix),pi=pi))
	nn<-vector(length=k,mode="character")
	for(i in 1:k) nn[i]<-paste("[",paste(states[which(rate==i,arr.ind=TRUE)[1,]],
		collapse=","),"]",sep="")
	PS<-matrix(NA,ngen,k+2,dimnames=list(1:ngen,c("gen",nn,"logLik")))
	PS[1,]<-c(1,q,likQ)
	cat("Running MCMC....\n")
	if(print){
		cat(paste(colnames(PS),collapse=" \t"))
		cat("\n")
	}
	flush.console()
	qp<-q
	for(i in 2:ngen){
		qp[i%%k+1]<-proposal(q[i%%k+1],prop.var[i%%k+1])
		likQp<-logLik(fitMk(tree,x,model=model,fixedQ=makeQ(m,qp,index.matrix),pi=pi))
		por<-exp(likQp-likQ+log.prior(qp,prior.rate)-log.prior(q,prior.rate))
		if(por>runif(n=1)){
			q<-qp
			likQ<-likQp
		}
		PS[i,]<-c(i,q,likQ)
		if(print) if(i%%print==0){
			cat(paste(round(PS[i,1]),paste(round(PS[i,1:k+1],4),collapse="\t"),
				round(PS[i,ncol(PS)],4),sep="\t"))
			cat("\n")
			flush.console()
		}
	}
	cat("Done.\n")
	class(PS)<-"mcmcMk"
	PS
}	

print.mcmcMk<-function(x,...){
	cat("\n   Posterior sample from mcmcMk consisting of a matrix.\n\n")
}

