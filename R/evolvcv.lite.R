# function is simplified version of evol.vcv
# written by Liam J. Revell 2011, 2012, 2013, 2017

evolvcv.lite<-function(tree,X,maxit=2000,tol=1e-10){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")

	# model 1: common variances & correlation
	lik1<-function(theta,C,D,y){
		v<-theta[1:2]; r<-theta[3]
		R<-matrix(c(v[1],r*sqrt(v[1]*v[2]),r*sqrt(v[1]*v[2]),v[2]),2,2)
		V<-kronecker(R,C)
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}

	# model 2: different variances, same correlation
	lik2<-function(theta,C,D,y){
		p<-(length(theta)-1)/2
		v<-matrix(theta[1:(2*p)],p,2,byrow=T); r<-theta[length(theta)]
		R<-list()
		for(i in 1:p) R[[i]]<-matrix(c(v[i,1],r*sqrt(v[i,1]*v[i,2]),r*sqrt(v[i,1]*v[i,2]),v[i,2]),2,2)
		V<-matrix(0,length(y),length(y))
		for(i in 1:p) V<-V+kronecker(R[[i]],C[[i]])
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}

	# model 3: same variances, different correlation
	lik3<-function(theta,C,D,y){
		p<-length(theta)-2
		v<-theta[1:2]; r<-theta[3:length(theta)]
		R<-list()
		for(i in 1:p) R[[i]]<-matrix(c(v[1],r[i]*sqrt(v[1]*v[2]),r[i]*sqrt(v[1]*v[2]),v[2]),2,2)
		V<-matrix(0,length(y),length(y))
		for(i in 1:p) V<-V+kronecker(R[[i]],C[[i]])
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}

	# model 4: everything different
	lik4<-function(theta,C,D,y){
		p<-length(theta)/3
		v<-matrix(theta[1:(2*p)],p,2,byrow=T); r<-theta[(2*p+1):length(theta)]
		R<-list()
		for(i in 1:p) R[[i]]<-matrix(c(v[i,1],r[i]*sqrt(v[i,1]*v[i,2]),r[i]*sqrt(v[i,1]*v[i,2]),v[i,2]),2,2)
		V<-matrix(0,length(y),length(y))
		for(i in 1:p) V<-V+kronecker(R[[i]],C[[i]])
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}
		
	# done internal functions

	# bookkeeping
	X<-as.matrix(X)
	X<-X[tree$tip.label,]
	n<-nrow(X) # number of species
	m<-ncol(X) # number of traits
	if(m!=2) stop("number of traits must equal 2")
	p<-ncol(tree$mapped.edge) # number of states
	D<-matrix(0,n*m,m)
	for(i in 1:(n*m)) for(j in 1:m) if((j-1)*n<i&&i<=j*n) D[i,j]=1.0
	y<-as.matrix(as.vector(X))
	mC<-multiC(tree)
	C<-Reduce("+",mC)
	sv1<-mean(pic(X[,1],tree)^2)
	sv2<-mean(pic(X[,2],tree)^2)
	sr<-mean(pic(X[,1],tree)*pic(X[,2],tree))/sqrt(sv1*sv2)

	# now optimize models
	res1<-optim(runif(3)*c(sv1,sv2,sr),lik1,C=C,D=D,y=y,method="L-BFGS-B",lower=tol+c(0,0,-1),upper=c(Inf,Inf,1)-tol)
	res2<-optim(runif(2*p+1)*c(rep(c(sv1,sv2),p),sr),lik2,C=mC,D=D,y=y,method="L-BFGS-B",lower=tol+c(rep(0,2*p),-1),upper=c(rep(Inf,2*p),1)-tol)
	res3<-optim(runif(2+p)*c(sv1,sv2,rep(sr,p)),lik3,C=mC,D=D,y=y,method="L-BFGS-B",lower=tol+c(0,0,rep(-1,p)),upper=c(Inf,Inf,rep(1,p))-tol)
	res4<-optim(runif(3*p)*c(rep(c(sv1,sv2),p),rep(sr,p)),lik4,C=mC,D=D,y=y,method="L-BFGS-B",lower=tol+c(rep(0,2*p),rep(-1,p)),upper=c(rep(Inf,2*p),rep(1,p))-tol)

	m1<-list(description="common rates, common correlation",
				R=matrix(c(res1$par[1],res1$par[3]*sqrt(res1$par[1]*res1$par[2]),res1$par[3]*sqrt(res1$par[1]*res1$par[2]),res1$par[2]),2,2),
				logLik=-res1$value,
				convergence=res1$convergence,
				k=length(res1$par)+2,
				AIC=2*(length(res1$par)+2)+2*res1$value)
	R<-list(); for(i in 1:p) R[[i]]<-matrix(c(res2$par[2*(i-1)+1],rep(res2$par[2*p+1]*sqrt(res2$par[2*(i-1)+1]*res2$par[2*(i-1)+2]),2),res2$par[2*(i-1)+2]),2,2)
	names(R)<-colnames(tree$mapped.edge)
	m2<-list(description="different rates, common correlation",
				R=R,
				logLik=-res2$value,
				convergence=res2$convergence,
				k=length(res2$par)+2,
				AIC=2*(length(res2$par)+2)+2*res2$value)
	R<-list(); for(i in 1:p) R[[i]]<-matrix(c(res3$par[1],rep(res3$par[2+i]*sqrt(res3$par[1]*res3$par[2]),2),res3$par[2]),2,2)
	names(R)<-colnames(tree$mapped.edge)
	m3<-list(description="common rates, different correlation",
				R=R,
				logLik=-res3$value,
				convergence=res3$convergence,
				k=length(res3$par)+2,
				AIC=2*(length(res3$par)+2)+2*res3$value)
	R<-list(); for(i in 1:p) R[[i]]<-matrix(c(res4$par[2*(i-1)+1],rep(res4$par[2*p+i]*sqrt(res4$par[2*(i-1)+1]*res4$par[2*(i-1)+2]),2),res4$par[2*(i-1)+2]),2,2)
	names(R)<-colnames(tree$mapped.edge)
	m4<-list(description="no common structure",
				R=R,
				logLik=-res4$value,
				convergence=res4$convergence,
				k=length(res4$par)+2,
				AIC=2*(length(res4$par)+2)+2*res4$value)
	
	obj<-list(model1=m1,model2=m2,model3=m3,model4=m4)
	class(obj)<-"evolvcv.lite"
	obj
}

## S3 print method for object of class "evolvcv.lite"
## written by Liam J. Revell 2017
print.evolvcv.lite<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	## MODEL 1	
	m1<-lapply(x$model1,function(a,b) if(is.numeric(a)) round(a,b) else a,b=digits)
	cat(paste("Model 1:",x$model1$description,"\n"))
	nn<-paste("R[",t(sapply(1:ncol(m1$R),paste,1:ncol(m1$R),sep=","))[upper.tri(m1$R,
		diag=TRUE)],"]",sep="")
	cat(paste("\t",paste(nn,collapse="\t"),"\tk\tlog(L)\tAIC","\n",sep=""))
	cat(paste("fitted",paste(m1$R[upper.tri(m1$R,diag=TRUE)],collapse="\t"),m1$k,
		m1$logLik,m1$AIC,"\n",sep="\t"))
	if(m1$convergence==0) cat("\n(R thinks it has found the ML solution for model 1.)\n\n")
	else cat("\n(Model 1 optimization may not have converged.)\n\n")
	## MODELS 2-4
	for(i in 2:4){
		m<-lapply(x[[i]],function(a,b) if(is.numeric(a)) round(a,b) else a,b=digits)
		m$R<-lapply(m$R,function(a,b) if(is.numeric(a)) round(a,b) else a,b=digits)
		cat(paste("Model ",i,": ",m$description,"\n",sep=""))
		cat(paste("\t",paste(nn,collapse="\t"),"\tk\tlog(L)\tAIC","\n",sep=""))
		for(j in 1:length(m$R)){
			if(j==1) cat(paste(names(m$R)[j],paste(m$R[[j]][upper.tri(m$R[[j]],diag=TRUE)],
				collapse="\t"),m$k,m$logLik,m$AIC,"\n",sep="\t"))
			else cat(paste(names(m$R)[j],paste(m$R[[j]][upper.tri(m$R[[j]],diag=TRUE)],
				collapse="\t"),"\n",sep="\t"))
		}
		if(m$convergence==0) 
			cat(paste("\n(R thinks it has found the ML solution for model ",i,".)\n\n",sep=""))
		else cat(paste("\n(Model ",i,"optimization may not have converged.)\n\n",sep=""))
	}
}

