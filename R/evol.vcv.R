# this function fits the model of Revell & Collar (2009; Evolution)
# written by Liam J. Revell 2010, 2011, 2013, 2014, 2015, 2016, 2020, 2021

evol.vcv<-function(tree,X,maxit=2000,vars=FALSE,...){

	## start checks
	
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")

	## end checks

	## likelihood functions (to be used internally)
	
	# compute the log-likelihood (from the cholesky matrices)
	lik.chol<-function(theta,y,C,D,E){
		m<-length(y)/dim(C[[1]])[1]
		n<-length(y)/m
		p<-length(C)
		cholR<-array(data=0,dim=c(m,m,p))
		l<-1
		for(i in 1:p) for(j in 1:m) for(k in j:m){ 
			cholR[j,k,i]<-theta[l]
			l<-l+1
		}
		V<-matrix(0,nrow(D),nrow(D))
		for(i in 1:p) V<-V+kronecker(t(cholR[,,i])%*%cholR[,,i],C[[i]])
		V<-V+E
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-
			determinant(V)$modulus[1]/2
		return(-logL)
	}
	
	# compute the log-likelihood (from the original matrices)
	lik.R<-function(theta,y,C,D,E){
		m<-length(y)/nrow(C[[1]])
		n<-length(y)/m
		p<-length(C)
		R<-array(data=0,dim=c(m,m,p)); l<-1
		for(i in 1:p) for(j in 1:m) for(k in j:m){ 
			R[j,k,i]<-theta[l]
			R[k,j,i]<-theta[l]
			l<-l+1
		}
		V<-matrix(0,nrow(D),nrow(D))
		for(i in 1:p) V<-V+kronecker(R[,,i],C[[i]])
		V<-V+E
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}
	
	## end likelihood functions
	
	## start preliminaries
	
	if(is.data.frame(X)) X<-as.matrix(X)
	n<-nrow(X) # number of species
	m<-ncol(X) # number of traits
	if(hasArg(err_vcv)){ 
		err_vcv<-list(...)$err_vcv
		err_vcv<-err_vcv[tree$tip.label]
	} else { 
		err_vcv<-replicate(n,matrix(0,m,m),simplify=FALSE)
		names(err_vcv)<-tree$tip.label
	}
	E<-matrix(0,n*m,n*m)
	for(i in 1:n){
		ii<-0:(m-1)*n+i
		E[ii,ii]<-err_vcv[[i]]
	}
	tree1rate<-paintSubTree(tree,Ntip(tree)+1,"fitted")
	tree<-if(inherits(tree,"simmap")) tree else tree1rate
	p<-ncol(tree$mapped.edge) # number of states
	D<-matrix(0,n*m,m)
	for(i in 1:(n*m)) for(j in 1:m) if((j-1)*n<i&&i<=j*n) D[i,j]=1.0
	X<-X[rownames(C<-vcv(tree)),,drop=FALSE]
	y<-as.matrix(as.vector(X))	
	
	## end preliminaries 
	
	## get starting parameter values
	
	a<-colSums(solve(C))%*%X/sum(solve(C)) # ancestral states
	R<-t(X-rep(1,nrow(X))%*%a)%*%solve(C)%*%(X-rep(1,nrow(X))%*%a)/n
	
	## end starting parameter values
	
	## start fit single matrix model
	
	C<-multiC(tree1rate)
		
	# populate the starting parameter vector
	ss<-vector(mode="numeric")
	for(i in 1:p) ss<-c(ss,to.upper(chol(R)))
	# optimize using generic optimizer
	r=optim(ss,fn=lik.chol,y=y,C=C,D=D,E=E,control=list(maxit=maxit))
	# update R
	R<-matrix(data=t(upper.diag(r$par[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))]))%*%
		upper.diag(r$par[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))]),m,m,dimnames=list(colnames(X),
		colnames(X)))
	dimnames(R)<-list(rownames(R),colnames(R))
	# log-likelihood for the single model
	logL1<--r$value
	# convert R to a vector
	Rv<-to.upper(R)
	H<-hessian(lik.R,Rv,y=y,C=C,D=D,E=E)
	# convert Hessian diagonal to matrices
	if(vars){ 
		print(H)
		Vh<-diag(solve(H))
		vars.single<-matrix(data=t(to.symmetric(Vh[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))])),m,m)
		dimnames(vars.single)<-list(rownames(R),colnames(R))
	}
	# report convergence
	convergence<-r$convergence
	
	## end fitting single matrix model
	
	## start fit multi-matrix model

	if(p>1){

		C<-multiC(tree) # compute separate C for each state
		
		# compute the starting parameter values
		ss<-vector(mode="numeric")
		for(i in 1:p) ss<-c(ss,to.upper(chol(R)))
		# optimize using generic optimizer
		r=optim(ss,fn=lik.chol,y=y,C=C,D=D,E=E,control=list(maxit=maxit))
		# convert parameter estimates to matrices
		R.i<-array(dim=c(m,m,p))
		ss<-colnames(tree$mapped.edge)
		for(i in 1:p) 
			R.i[,,i]<-matrix(data=t(upper.diag(r$par[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))]))%*%
				upper.diag(r$par[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))]),m,m,dimnames=list(colnames(X),
				colnames(X)))
		dimnames(R.i)<-list(rownames(R),colnames(R),ss)
		# log-likelihood for the multi-matrix model
		logL2<--r$value

		# convert R.i to a vector
		Rv<-vector(mode="numeric")
		for(i in 1:p) Rv<-c(Rv,to.upper(R.i[,,i]))
		H<-hessian(lik.R,Rv,y=y,C=C,D=D,E=E)
		# convert Hessian diagonal to matrices
		if(vars){ 
			Vars<-array(dim=c(m,m,p))
			Vh<-diag(solve(H))
			for(i in 1:p) Vars[,,i]<-matrix(data=t(to.symmetric(Vh[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))])),m,m)
			dimnames(Vars)<-list(rownames(R),colnames(R),ss)
		}
		# report convergence
		convergence[2]==r$convergence
	}
	if(all(convergence==0)) converged<-"Optimization has converged."
	else converged<-"Optimization may not have converged. Consider increasing maxit."
	
	if(p>1){
		if(vars) obj<-list(R.single=R,vars.single=vars.single,logL1=as.numeric(logL1),
			k1=m*(m+1)/2+m,R.multiple=R.i,vars.multiple=Vars,logL.multiple=logL2,
			k2=p*m*(m+1)/2+m,P.chisq=pchisq(2*(logL2-as.numeric(logL1)),(p-1)*m*(m+1)/2,
			lower.tail=FALSE),convergence=converged)
		else obj<-list(R.single=R,logL1=as.numeric(logL1),k1=m*(m+1)/2+m,R.multiple=R.i,
			logL.multiple=logL2,k2=p*m*(m+1)/2+m,P.chisq=pchisq(2*(logL2-as.numeric(logL1)),
			(p-1)*m*(m+1)/2,lower.tail=FALSE),convergence=converged)
	} else {
		if(vars) obj<-list(R.single=R,vars.single=vars.single,logL1=as.numeric(logL1),
			k1=m*(m+1)/2+m,convergence="Optimization has converged.")
		else obj<-list(R.single=R,logL1=as.numeric(logL1),k1=m*(m+1)/2+m,
			convergence="Optimization has converged.")
	}
	class(obj)<-"evol.vcv"
	obj
}

# function puts the upper triangle of a square matrix in a vector, by row
# written by Liam J. Revell 2013
to.upper<-function(X) t(X)[lower.tri(X,diag=TRUE)]

# function puts a vector into the upper triangle of a square matrix, by row
# written by Liam J. Revell 2013
upper.diag<-function(x){
	m<-(-1+sqrt(1+8*length(x)))/2
	X<-lower.tri(matrix(NA,m,m),diag=TRUE)
	X[X==TRUE]<-x
	t(X)
}

# function converts vector to symmetric matrix
# written by Liam J. Revell 2013
to.symmetric<-function(x){
	if(length(x)==1) X<-matrix(x,1,1)
	else {
		X<-upper.diag(x)
		for(i in 2:nrow(X)) for(j in 1:(i-1)) X[i,j]<-X[j,i]
	}
	X
}

## S3 print method for object of class "evol.vcv"
## written by Liam J. Revell 2013
print.evol.vcv<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	x<-lapply(x,function(a,b) if(is.numeric(a)) round(a,b) else a,b=digits)
	cat("ML single-matrix model:\n")
	nn<-paste("R[",t(sapply(1:ncol(x$R.single),paste,
		1:ncol(x$R.single),sep=","))[upper.tri(x$R.single,diag=TRUE)],"]",
		sep="")
	cat(paste("\t",paste(nn,collapse="\t"),"\tk\tlog(L)","\n",sep=""))
	cat(paste("fitted",paste(x$R.single[upper.tri(x$R.single,diag=TRUE)],
		collapse="\t"),x$k1,x$logL1,"\n",sep="\t"))
	cat("\n")
	if(!is.null(x$R.multiple)){
		cat("ML multi-matrix model:\n")
		cat(paste("\t",paste(nn,collapse="\t"),"\tk\tlog(L)","\n",sep=""))
		for(i in 1:dim(x$R.multiple)[3]){
			if(i==1) cat(paste(dimnames(x$R.multiple)[[3]][i],
				paste(x$R.multiple[,,i][upper.tri(x$R.single,diag=TRUE)],
				collapse="\t"),x$k2,x$logL.multiple,"\n",sep="\t"))
			else cat(paste(dimnames(x$R.multiple)[[3]][i],paste(x$R.multiple[,,i][upper.tri(x$R.single,
				diag=TRUE)],collapse="\t"),"\n",sep="\t"))
		}
		cat("\n")
		cat(paste("P-value (based on X^2):",x$P.chisq,"\n\n"))
	}
	if(x$convergence[1]=="Optimization has converged.") cat("R thinks it has found the ML solution.\n\n")
	else cat("Optimization may not have converged.\n\n")
}
