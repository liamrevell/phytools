# function creates a stochastic character mapped tree as a modified "phylo" object
# written by Liam Revell 2013, 2014, 2015

make.simmap<-function(tree,x,model="SYM",nsim=1,...){
	if(inherits(tree,"multiPhylo")){
		ff<-function(yy,x,model,nsim,...){
			zz<-make.simmap(yy,x,model,nsim,...)
			if(nsim>1) class(zz)<-NULL
			return(zz)
		}	
		if(nsim>1) mtrees<-unlist(sapply(tree,ff,x,model,nsim,...,simplify=FALSE),recursive=FALSE)
		else mtrees<-sapply(tree,ff,x,model,nsim,...,simplify=FALSE)
		class(mtrees)<-c("multiSimmap","multiPhylo")
	} else {
		## get optional arguments
		if(hasArg(pi)) pi<-list(...)$pi
		else pi<-"equal"
		if(hasArg(message)) pm<-list(...)$message
		else pm<-TRUE
		if(hasArg(tol)) tol<-list(...)$tol
		else tol<-0
		if(hasArg(Q)) Q<-list(...)$Q
		else Q<-"empirical"
		if(hasArg(burnin)) burnin<-list(...)$burnin
		else burnin<-1000
		if(hasArg(samplefreq)) samplefreq<-list(...)$samplefreq
		else samplefreq<-100
		if(hasArg(vQ)) vQ<-list(...)$vQ
		else vQ<-0.1
		prior<-list(alpha=1,beta=1,use.empirical=FALSE)
		if(hasArg(prior)){ 
			pr<-list(...)$prior
			prior[names(pr)]<-pr
		}
		## done optional arguments
		# check
		if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
		# if vector convert to binary matrix
		if(!is.matrix(x)) xx<-to.matrix(x,sort(unique(x)))
		else xx<-x
		xx<-xx[tree$tip.label,]
		xx<-xx/rowSums(xx)
		# reorder to cladewise
		tree<-bt<-reorder.phylo(tree,"cladewise")
		if(!is.binary.tree(bt)) bt<-multi2di(bt)
		# some preliminaries
		N<-length(tree$tip)
		m<-ncol(xx)
		root<-N+1
		# get conditional likelihoods & model
		if(is.character(Q)&&Q=="empirical"){
			XX<-getPars(bt,xx,model,Q=NULL,tree,tol,m,pi=pi,...)
			L<-XX$L
			Q<-XX$Q
			logL<-XX$loglik
			if(pi[1]=="equal") pi<-setNames(rep(1/m,m),colnames(L)) # set equal
			else if(pi[1]=="estimated") pi<-statdist(Q) # set from stationary distribution
			else pi<-pi/sum(pi) # obtain from input
			if(pm) printmessage(Q,pi,method="empirical")
			mtrees<-replicate(nsim,smap(tree,x,N,m,root,L,Q,pi,logL),simplify=FALSE)
		} else if(is.character(Q)&&Q=="mcmc"){
			if(prior$use.empirical){
				qq<-fitMk(bt,xx,model)$rates
				prior$alpha<-qq*prior$beta
			}
			XX<-mcmcQ(bt,xx,model,tree,tol,m,burnin,samplefreq,nsim,vQ,prior)
			L<-lapply(XX,function(x) x$L)
			Q<-lapply(XX,function(x) x$Q)
			logL<-lapply(XX,function(x) x$loglik)
			if(pi[1]=="equal"){
				pi<-setNames(rep(1/m,m),colnames(L)) # set equal
				pi<-lapply(1:nsim,function(x,y) y, y=pi)
			} else if(pi[1]=="estimated"){
				pi<-lapply(Q,statdist) # set from stationary distribution
			} else { 
				pi<-pi/sum(pi) # obtain from input
				pi<-lapply(1:nsim,function(x,y) y, y=pi)
			}
			if(pm) printmessage(Reduce('+',Q)/length(Q),pi,method="mcmc")
			mtrees<-mapply(smap,L=L,Q=Q,pi=pi,logL=logL,MoreArgs=list(tree=tree,x=x,N=N,m=m,root=root),SIMPLIFY=FALSE)
		} else if(is.matrix(Q)){
			XX<-getPars(bt,xx,model,Q=Q,tree,tol,m,pi=pi)
			L<-XX$L
			logL<-XX$loglik
			if(pi[1]=="equal") pi<-setNames(rep(1/m,m),colnames(L)) # set equal
			else if(pi[1]=="estimated") pi<-statdist(Q) # set from stationary distribution
			else pi<-pi/sum(pi) # obtain from input
			if(pm) printmessage(Q,pi,method="fixed")
			mtrees<-replicate(nsim,smap(tree,x,N,m,root,L,Q,pi,logL),simplify=FALSE)
		}
		if(length(mtrees)==1) mtrees<-mtrees[[1]]
		else class(mtrees)<-c("multiSimmap","multiPhylo")
	}	
	(if(hasArg(message)) list(...)$message else TRUE)
	if((if(hasArg(message)) list(...)$message else TRUE)&&inherits(tree,"phylo")) message("Done.")
	return(mtrees)
}

# print message
# written by Liam J. Revell 2013
printmessage<-function(Q,pi,method){
	if(method=="empirical"||method=="fixed")
		cat("make.simmap is sampling character histories conditioned on the transition matrix\n\nQ =\n")
	else if(method=="mcmc"){
		cat("make.simmap is simulating with a sample of Q from\nthe posterior distribution\n")
		cat("\nMean Q from the posterior is\nQ =\n")
	}
	print(Q)
	if(method=="empirical") cat("(estimated using likelihood);\n")
	else if(method=="fixed") cat("(specified by the user);\n")
	cat("and (mean) root node prior probabilities\npi =\n")
	if(is.list(pi)) pi<-Reduce("+",pi)/length(pi)
	print(pi)
	flush.console()
}

# mcmc for Q used in Q="mcmc"
# written by Liam J. Revell 2013
mcmcQ<-function(bt,xx,model,tree,tol,m,burnin,samplefreq,nsim,vQ,prior){
	update<-function(x){
		## x<-exp(log(x)+rnorm(n=np,mean=0,sd=sqrt(vQ)))
		x<-abs(x+rnorm(n=np,mean=0,sd=sqrt(vQ)))
		return(x)
	}
	# get model matrix
	if(is.character(model)){
		rate<-matrix(NA,m,m)
		if(model=="ER"){
			np<-rate[]<-1
			diag(rate)<-NA
		}
		if(model=="ARD"){
			np<-m*(m-1)
			rate[col(rate)!=row(rate)]<-1:np
		}
		if(model=="SYM") {
			np<-m*(m-1)/2
			sel<-col(rate)<row(rate)
			rate[sel]<-1:np
			rate<-t(rate)
			rate[sel]<-1:np
		}
	} else {
		if(ncol(model)!=nrow(model)) stop("the matrix given as 'model' is not square")
		if(ncol(model)!=m) stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
		rate<-model
		np<-max(rate)
	}
	# burn-in
	p<-rgamma(np,prior$alpha,prior$beta)
	Q<-matrix(c(0,p)[rate+1],m,m)
	diag(Q)<--rowSums(Q,na.rm=TRUE)
	yy<-getPars(bt,xx,model,Q,tree,tol,m,pi=pi)
	cat("Running MCMC burn-in. Please wait....\n")
	flush.console()
	for(i in 1:burnin){
		pp<-update(p)
		Qp<-matrix(c(0,pp)[rate+1],m,m)
		diag(Qp)<--rowSums(Qp,na.rm=TRUE)
		zz<-getPars(bt,xx,model,Qp,tree,tol,m,FALSE,pi=pi)
		p.odds<-exp(zz$loglik+sum(dgamma(pp,prior$alpha,prior$beta,log=TRUE))-
			yy$loglik-sum(dgamma(p,prior$alpha,prior$beta,log=TRUE)))
		if(p.odds>=runif(n=1)){
			yy<-zz
			p<-pp
		}
	}
	# now run MCMC generation, sampling at samplefreq
	cat(paste("Running",samplefreq*nsim,"generations of MCMC, sampling every",samplefreq,"generations.\nPlease wait....\n\n"))
	flush.console()
	XX<-vector("list",nsim)
	for(i in 1:(samplefreq*nsim)){
		pp<-update(p)
		Qp<-matrix(c(0,pp)[rate+1],m,m)
		diag(Qp)<--rowSums(Qp,na.rm=TRUE)
		zz<-getPars(bt,xx,model,Qp,tree,tol,m,FALSE,pi=pi)
		p.odds<-exp(zz$loglik+sum(dgamma(pp,prior$alpha,prior$beta,log=TRUE))-
			yy$loglik-sum(dgamma(p,prior$alpha,prior$beta,log=TRUE)))
		if(p.odds>=runif(n=1)){
			yy<-zz
			p<-pp
		}
		if(i%%samplefreq==0){
			Qi<-matrix(c(0,p)[rate+1],m,m)
			diag(Qi)<--rowSums(Qi,na.rm=TRUE)
			XX[[i/samplefreq]]<-getPars(bt,xx,model,Qi,tree,tol,m,TRUE,pi=pi)
		}
	}
	return(XX)
}

# get pars
# written by Liam J. Revell 2013
getPars<-function(bt,xx,model,Q,tree,tol,m,liks=TRUE,pi,...){
	obj<-fitMk(bt,xx,model,fixedQ=Q,output.liks=liks,pi=pi,...)
	N<-length(bt$tip.label)
	II<-obj$index.matrix+1
	lvls<-obj$states
	if(liks){
		L<-obj$lik.anc
		rownames(L)<-N+1:nrow(L)
		if(!is.binary.tree(tree)){
			ancNames<-matchNodes(tree,bt)
			L<-L[as.character(ancNames[,2]),]
			rownames(L)<-ancNames[,1]
		}
		L<-rbind(xx,L)
		rownames(L)[1:N]<-1:N
	} else L<-NULL	
	Q<-matrix(c(0,obj$rates)[II],m,m,dimnames=list(lvls,lvls))
	if(any(rowSums(Q,na.rm=TRUE)<tol)){
		message(paste("\nWarning: some rows of Q not numerically distinct from 0; setting to",tol,"\n"))
		ii<-which(rowSums(Q,na.rm=TRUE)<tol)
		for(i in 1:length(ii)) Q[ii[i],setdiff(1:ncol(Q),ii[i])]<-tol/(ncol(Q)-1)
	}
	diag(Q)<--rowSums(Q,na.rm=TRUE)
	return(list(Q=Q,L=L,loglik=logLik(obj)))
}

# convert vector of x to binary matrix
# written by Liam J. Revell 2012
to.matrix<-function(x,seq){
	X<-matrix(0,length(x),length(seq),dimnames=list(names(x),seq))
	for(i in 1:length(seq)) X[x==seq[i],i]<-1
	return(X)
}

# function does the stochastic mapping, conditioned on our model & given the conditional likelihoods
# written by Liam J. Revell 2013
smap<-function(tree,x,N,m,root,L,Q,pi,logL){
	# create the map tree object
	mtree<-tree; mtree$maps<-list()
	mtree$mapped.edge<-matrix(0,nrow(tree$edge),m,dimnames=list(paste(tree$edge[,1],",",tree$edge[,2],sep=""),colnames(L)))
	# now we want to simulate the node states & histories by pre-order traversal
	NN<-matrix(NA,nrow(tree$edge),2) # our node values
	NN[which(tree$edge[,1]==root),1]<-rstate(L[as.character(root),]*pi/sum(L[as.character(root),]*pi)) # assign root
	for(j in 1:nrow(tree$edge)){
		# conditioned on the start value, assign end value of node (if internal)
		p<-expm(Q*tree$edge.length[j])[NN[j,1],]*L[as.character(tree$edge[j,2]),]
		NN[j,2]<-rstate(p/sum(p))
		NN[which(tree$edge[,1]==tree$edge[j,2]),1]<-NN[j,2]
		# now simulate on the branches
		accept<-FALSE	
		while(!accept){
			map<-sch(NN[j,1],tree$edge.length[j],Q)
			if(names(map)[length(map)]==NN[j,2]) accept=TRUE
		}	
		mtree$maps[[j]]<-map
		for(k in 1:length(mtree$maps[[j]])) 
			mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]<-mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]+mtree$maps[[j]][k]
	}
	mtree$Q<-Q
	mtree$logL<-logL
	if(!inherits(mtree,"simmap")) class(mtree)<-c("simmap",setdiff(class(mtree),"simmap"))
	return(mtree)
}

# function generates a history along a branch
# written by Liam J. Revell 2013
sch<-function(start,t,Q){
	tol<-t*1e-12
	dt<-setNames(0,start)
	while(sum(dt)<(t-tol)){
		s<-names(dt)[length(dt)]
		dt[length(dt)]<-if(-Q[s,s]>0) rexp(n=1,rate=-Q[s,s]) else t-sum(dt)
		if(sum(dt)<(t-tol)){
			dt<-c(dt,0)
			if(sum(Q[s,][-match(s,colnames(Q))])>0)
				names(dt)[length(dt)]<-rstate(Q[s,][-match(s,colnames(Q))]/sum(Q[s,][-match(s,colnames(Q))]))
			else names(dt)[length(dt)]<-s
		} else dt[length(dt)]<-dt[length(dt)]-sum(dt)+t
	}
	return(dt)
}

# function uses numerical optimization to solve for the stationary distribution
# written by Liam J. Revell 2013
statdist<-function(Q){
	foo<-function(theta,Q){
		Pi<-c(theta[1:(nrow(Q)-1)],1-sum(theta[1:(nrow(Q)-1)]))
		sum((Pi%*%Q)^2)
	}
	k<-nrow(Q)
	if(nrow(Q)>2){ 
		fit<-optim(rep(1/k,k-1),foo,Q=Q,control=list(reltol=1e-16))
		return(setNames(c(fit$par[1:(k-1)],1-sum(fit$par[1:(k-1)])),rownames(Q)))
	} else {
		fit<-optimize(foo,interval=c(0,1),Q=Q)
		return(setNames(c(fit$minimum,1-fit$minimum),rownames(Q)))
	}
}

## S3 print method for objects of class "simmap" & multiSimmap
## based on print.phylo in ape
print.simmap<-function(x,printlen=6,...){
	N<-Ntip(x)
    	M<-x$Nnode
	cat(paste("\nPhylogenetic tree with",N,"tips and",M,"internal nodes.\n\n"))
	cat("Tip labels:\n")
	if(N>printlen) cat(paste("\t",paste(x$tip.label[1:printlen],collapse=", "),", ...\n",sep=""))
	else print(x$tip.label)
	ss<-sort(unique(c(getStates(x,"tips"),getStates(x,"nodes"))))
	cat(paste("\nThe tree includes a mapped, ",length(ss),"-state discrete character with states:\n",
		sep=""))
	if(length(ss)>printlen) cat(paste("\t",paste(ss[1:printlen],collapse=", "),", ...\n",sep=""))
	else cat(paste("\t",paste(ss,collapse=", "),"\n",sep=""))
    	rlab<-if(is.rooted(x)) "Rooted" else "Unrooted"
	cat("\n",rlab,"; includes branch lengths.\n",sep="")
}
print.multiSimmap<-function(x,details=FALSE,...){
	N<-length(x)
	cat(N,"phylogenetic trees with mapped discrete characters\n")
    	if(details){
		n<-sapply(x,Ntip)
		s<-sapply(x,function(x) length(unique(c(getStates(x,"tips"),getStates(x,"nodes")))))
		for(i in 1:N) cat("tree",i,":",n[i],"tips,",s[i],"mapped states\n")
	}
}

## S3 summary method for objects of class "simmap" & "multiSimmap"
summary.simmap<-function(object,...) describe.simmap(object,...)
summary.multiSimmap<-function(object,...) describe.simmap(object,...)

## for backward compatibility with any function using apeAce internally
apeAce<-function(tree,x,model,fixedQ=NULL,...){
	if(hasArg(output.liks)){ 
		output.liks<-list(...)$output.liks
		return(fitMk(tree,x,model,fixedQ,...))
	} else { 
		output.liks<-TRUE
		return(fitMk(tree,x,model,fixedQ,output.liks=TRUE,...))
	}
}

## S3 logLik methods for "simmap" & "multiSimmap"
logLik.simmap<-function(object,...) object$logL
logLik.multiSimmap<-function(object,...)
	sapply(object,function(x) x$logL)