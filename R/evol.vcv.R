# this function fits the model of Revell & Collar (2009; Evolution)
# written by Liam J. Revell 2010, 2011, 2013, 2014, 2015, 2016, 2020

evol.vcv<-function(tree,X,maxit=2000,vars=FALSE,...){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	if(is.data.frame(X)) X<-as.matrix(X)
	n<-nrow(X) # number of species
	m<-ncol(X) # number of traits
	if(hasArg(se)){ 
		se<-list(...)$se
		se<-se[tree$tip.label]
	} else { 
		se<-replicate(n,matrix(0,m,m),simplify=FALSE)
		names(se)<-tree$tip.label
	}
	SE<-matrix(0,n*m,n*m)
	for(i in 1:n){
		ii<-0:(m-1)*n+i
		SE[ii,ii]<-se[[i]]
	}
	p<-ncol(tree$mapped.edge) # number of states
	D<-matrix(0,n*m,m)
	for(i in 1:(n*m)) for(j in 1:m) if((j-1)*n<i&&i<=j*n) D[i,j]=1.0
	C<-vcv(tree)
	X<-X[rownames(C),,drop=FALSE]
	y<-as.matrix(as.vector(X))
	a<-colSums(solve(C))%*%X/sum(solve(C)) # ancestral states
	R<-t(X-rep(1,nrow(X))%*%a)%*%solve(C)%*%(X-rep(1,nrow(X))%*%a)/n
	# likelihood for a single matrix
	lik1<-function(rr,y,D,a,C,n,mm){
		R<-to.symmetric(rr)
		return(-t(y-D%*%t(a))%*%solve(kronecker(R,C))%*%(y-D%*%t(a))/2-n*mm*log(2*pi)/2-
			determinant(kronecker(R,C))$modulus[1]/2)
	}
	logL1<-lik1(to.upper(R),y,D,a,C,n,m)
	if(vars){
		H<-hessian(lik1,to.upper(R),y=y,D=D,a=a,C=C,n=n,mm=m) # compute Hessian
		vars.single<-to.symmetric(diag(solve(-H))) # compute variances
		dimnames(vars.single)<-list(colnames(X),colnames(X))
	}
	if(inherits(tree,"simmap")){
		C<-multiC(tree) # compute separate C for each state
		# compute the log-likelihood (from the cholesky matrices)
		lik.cholR<-function(theta,y,C,D){
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
			a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
			logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-
				determinant(V)$modulus[1]/2
			return(-logL)
		}
		# compute the starting parameter values
		ss<-vector(mode="numeric")
		for(i in 1:p) ss<-c(ss,to.upper(chol(R)))
		# optimize using generic optimizer
		r=optim(ss,fn=lik.cholR,y=y,C=C,D=D,control=list(maxit=maxit))
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
		# compute the log-likelihood (from the original matrices)
		lik2<-function(theta,y,C,D){
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
			a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
			logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
			return(-logL)
		}
		# convert R.i to a vector
		Rv<-vector(mode="numeric")
		for(i in 1:p) Rv<-c(Rv,to.upper(R.i[,,i]))
		H<-hessian(lik2,Rv,y=y,C=C,D=D)
		# convert Hessian diagonal to matrices
		if(vars){ 
			Vars<-array(dim=c(m,m,p))
			Vh<-diag(solve(H))
			for(i in 1:p) Vars[,,i]<-matrix(data=t(to.symmetric(Vh[((i-1)*m*(m+1)/2+1):(i*(m*(m+1)/2))])),m,m)
			dimnames(Vars)<-list(rownames(R),colnames(R),ss)
		}
		# report convergence
		if(r$convergence==0) converged<-"Optimization has converged."
		else converged<-"Optimization may not have converged. Consider increasing maxit."
	}
	if(inherits(tree,"simmap")){
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
