## function to impute missing values in a multivariate data matrix
## written by Liam J. Revell 2019

phylo.impute<-function(tree,X,...){
	if(hasArg(maxit)) maxit<-list(...)$maxit
	else maxit<-5000
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-TRUE
	if(hasArg(fixed)) fixed<-list(...)$fixed
	else fixed<-FALSE
	if(hasArg(p)) p<-list(...)$p
	else p<-2.0
	if(is.data.frame(X)) X<-as.matrix(X)
	ii<-which(is.na(X),arr.ind=TRUE)
	lik<-function(theta,tree,X,ii){
		X[ii]<-theta
		object<-evol.vcv(tree,X)
		if(!quiet) print(object$logL1)
		-object$logL1
	}
	lower<-apply(X,2,expand.range,na.rm=TRUE,p=p)[1,ii[,2]]
	upper<-apply(X,2,expand.range,na.rm=TRUE,p=p)[2,ii[,2]]
	if(hasArg(x.init)) x.init<-list(...)$x.init
	else x.init<-"ace"
	if(x.init[1]=="random") x.init<-runif(n=length(ii),lower,upper)
	else if(x.init[1]=="ace") { 
		x.init<-vector()
		for (i in 1:nrow(ii)){
			tip<-rownames(X)[ii[i]]
			tt<-drop.tip(root(tree,outgroup=tip),
				names(ii[which(ii[,2]==ii[i,2]),1]))
			x.init[i]<-fastAnc(tt,X[!is.na(X[,ii[i,2]]),
				ii[i,2]])[1]
		}
	}
	if(fixed==FALSE){
		fit<-optim(x.init,lik,tree=tree,X=X,ii=ii,method="L-BFGS-B",
			lower=lower,upper=upper,control=list(maxit=maxit))
	} else {
		fit<-list(par=x.init,value=lik(x.init,tree,X,ii),
		counts=c(0,0),convergence=0,
		message="Fixed value of par.")
	}
	X[ii]<-fit$par
	attr(X,"optim")<-list(logLik=-fit$value,counts=fit$counts,
		convergence=fit$convergence,message=fit$message)
	class(X)<-c("matrix","phylo.impute")
	X
}

expand.range<-function(x,na.rm=FALSE,p=2.0){
	rr<-range(x,na.rm=na.rm)
	mean(rr)+c(p*(rr[1]-mean(rr)),p*(rr[2]-mean(rr)))
}

print.phylo.impute<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-c(7,3)
	if(hasArg(n)) n<-list(...)$n
	else n<-6L
	if(n>nrow(x)) n<-nrow(x)
	cat("Results from phylogenetic imputation:\n")
	print(head(unclass(x),n),digits=digits[1])
	if(n<nrow(x)) cat("...\n\n") else cat("\n")
	cat(paste("log(L) :",round(attr(x,"optim")$logLik,digits[2]),"\n"))
	if(attr(x,"optim")$convergence==0)
		cat("R thinks it has converged on the correct solution.\n\n")
	else
		cat("R thinks it may not have converged.\n\n")
}