## function is simplified version of evol.vcv
## written by Liam J. Revell 2011, 2012, 2013, 2017, 2019, 2020, 2021, 2022

anova.evolvcv.lite<-function(object,...){
	models<-paste("model",sapply(strsplit(names(object),"model"),
		function(x) x[2]))
	logL<-sapply(object,function(x) x$logLik)
	df<-sapply(object,function(x) x$k)
	AIC<-sapply(object,function(x) x$AIC)
	value<-data.frame(logL=logL,df=df,AIC=AIC,weight=unclass(aic.w(AIC)))
	rownames(value)<-models
	colnames(value)<-c("log(L)","d.f.","AIC","weight")
	print(value)
	invisible(value)
}

evolvcv.lite<-function(tree,X,maxit=2000,tol=1e-10,...){
	## check 'phylo' object
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	## get some optional arguments
	if(hasArg(models)) models<-list(...)$models
	else models<-as.character(1:4)
	if(models[1]=="all models") models<-c("1","2","2b","2c","3","3b","3c","4")
	models<-as.character(models)
	if(hasArg(try.iter)) try.iter<-list(...)$try.iter
	else try.iter<-10
	if(hasArg(lower)){
		lower<-list(...)$lower
		if(length(lower)==1) lower<-c(rep(lower,2),-1)
		if(length(lower)>3) lower<-lower[1:3]
	} else lower<-NULL
	if(hasArg(upper)){
		upper<-list(...)$upper
		if(length(upper)==1) upper<-c(rep(upper,2),1)
		if(length(upper)>3) upper<-upper[1:3]
	} else upper<-NULL
	
	if(!inherits(tree,"simmap")) models<-intersect("1",models)
	
	if(length(models)==0) stop("for models!=\"1\" tree must be an object of class \"simmap\".")

	# model 1: common variances & correlation
	lik1<-function(theta,C,D,y,E){
		v<-exp(theta[1:2])
		r<-theta[3]
		R<-matrix(c(v[1],r*sqrt(v[1]*v[2]),r*sqrt(v[1]*v[2]),
			v[2]),2,2)
		V<-kronecker(R,C)+E
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-
			n*m*log(2*pi)/2-determinant(V)$modulus[1]/2
		return(-logL)
	}
	
	# model 2: different variances, same correlation
	lik2<-function(theta,C,D,y,E){
		p<-(length(theta)-1)/2
		v<-matrix(exp(theta[1:(2*p)]),p,2,byrow=T)
		r<-theta[length(theta)]
		R<-list()
		for(i in 1:p) R[[i]]<-matrix(c(v[i,1],r*sqrt(v[i,1]*v[i,2]),
			r*sqrt(v[i,1]*v[i,2]),v[i,2]),2,2)
		V<-matrix(0,length(y),length(y))
		for(i in 1:p) V<-V+kronecker(R[[i]],C[[i]])
		V<-V+E
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-
			determinant(V)$modulus[1]/2
		return(-logL)
	}
	
	## model 2b: different variances for only trait 1, same correlation
	lik2b<-function(theta,C,D,y,E){
		p<-length(theta)-2
		v<-matrix(exp(c(theta[1:p],rep(theta[p+1],p))),p,2,byrow=FALSE)
		r<-theta[length(theta)]
		R<-list()
		for(i in 1:p) R[[i]]<-matrix(c(v[i,1],r*sqrt(v[i,1]*v[i,2]),
			r*sqrt(v[i,1]*v[i,2]),v[i,2]),2,2)
		V<-matrix(0,length(y),length(y))
		for(i in 1:p) V<-V+kronecker(R[[i]],C[[i]])
		V<-V+E
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-
			determinant(V)$modulus[1]/2
		return(-logL)
	}
	
	## model 2c: different variances for only trait 2, same correlation
	lik2c<-function(theta,C,D,y,E){
		p<-length(theta)-2
		v<-matrix(exp(c(rep(theta[1],p),theta[1:p+1])),p,2,byrow=FALSE)
		r<-theta[length(theta)]
		R<-list()
		for(i in 1:p) R[[i]]<-matrix(c(v[i,1],r*sqrt(v[i,1]*v[i,2]),
			r*sqrt(v[i,1]*v[i,2]),v[i,2]),2,2)
		V<-matrix(0,length(y),length(y))
		for(i in 1:p) V<-V+kronecker(R[[i]],C[[i]])
		V<-V+E
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-
			determinant(V)$modulus[1]/2
		return(-logL)
	}	

	# model 3: same variances, different correlations
	lik3<-function(theta,C,D,y,E){
		p<-length(theta)-2
		v<-exp(theta[1:2])
		r<-theta[3:length(theta)]
		R<-list()
		for(i in 1:p) R[[i]]<-matrix(c(v[1],r[i]*sqrt(v[1]*v[2]),
			r[i]*sqrt(v[1]*v[2]),v[2]),2,2)
		V<-matrix(0,length(y),length(y))
		for(i in 1:p) V<-V+kronecker(R[[i]],C[[i]])
		V<-V+E
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-
			determinant(V)$modulus[1]/2
		return(-logL)
	}
	
	# model 3b: different variances for only trait 1, different correlation
	lik3b<-function(theta,C,D,y,E){
		p<-(length(theta)-1)/2
		v<-matrix(exp(c(theta[1:p],rep(theta[p+1],p))),p,2,byrow=FALSE)
		r<-theta[1:p+(p+1)]
		R<-list()
		for(i in 1:p) R[[i]]<-matrix(c(v[i,1],r[i]*sqrt(v[i,1]*v[i,2]),
			r[i]*sqrt(v[i,1]*v[i,2]),v[i,2]),2,2)
		V<-matrix(0,length(y),length(y))
		for(i in 1:p) V<-V+kronecker(R[[i]],C[[i]])
		V<-V+E
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-
			determinant(V)$modulus[1]/2
		return(-logL)
	}
	
	# model 3c: different variances for only trait 2, different correlations
	lik3c<-function(theta,C,D,y,E){
		p<-(length(theta)-1)/2
		v<-matrix(exp(c(rep(theta[1],p),theta[1:p+1])),p,2,byrow=FALSE)
		r<-theta[1:p+(p+1)]
		R<-list()
		for(i in 1:p) R[[i]]<-matrix(c(v[i,1],r[i]*sqrt(v[i,1]*v[i,2]),
			r[i]*sqrt(v[i,1]*v[i,2]),v[i,2]),2,2)
		V<-matrix(0,length(y),length(y))
		for(i in 1:p) V<-V+kronecker(R[[i]],C[[i]])
		V<-V+E
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-
			determinant(V)$modulus[1]/2
		return(-logL)
	}

	# model 4: everything different
	lik4<-function(theta,C,D,y,E){
		p<-length(theta)/3
		v<-matrix(exp(theta[1:(2*p)]),p,2,byrow=TRUE)
		r<-theta[(2*p+1):length(theta)]
		R<-list()
		for(i in 1:p) R[[i]]<-matrix(c(v[i,1],r[i]*sqrt(v[i,1]*v[i,2]),
			r[i]*sqrt(v[i,1]*v[i,2]),v[i,2]),2,2)
		V<-matrix(0,length(y),length(y))
		for(i in 1:p) V<-V+kronecker(R[[i]],C[[i]])
		V<-V+E
		a<-solve(t(D)%*%solve(V)%*%D)%*%(t(D)%*%solve(V)%*%y)
		logL<--t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*m*log(2*pi)/2-
			determinant(V)$modulus[1]/2
		return(-logL)
	}
		
	# done internal functions

	# bookkeeping
	X<-as.matrix(X)
	X<-X[tree$tip.label,]
	n<-nrow(X) # number of species
	m<-ncol(X) # number of traits
	if(m!=2) stop("number of traits must equal 2")
	
	## get error covariance matrices
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
	## end get error covariances
	
	## more bookkeeping
	p<-if(inherits(tree,"simmap")) ncol(tree$mapped.edge) else 1 # number of states
	D<-matrix(0,n*m,m)
	for(i in 1:(n*m)) for(j in 1:m) 
		if((j-1)*n<i&&i<=j*n) D[i,j]=1.0
	y<-as.matrix(as.vector(X))
	if(inherits(tree,"simmap")){ 
		mC<-multiC(tree)
		C<-Reduce("+",mC)
	} else C<-vcv.phylo(tree)
	sv1<-mean(pic(X[,1],multi2di(as.phylo(tree)))^2)
	sv2<-mean(pic(X[,2],multi2di(as.phylo(tree)))^2)
	sr<-mean(pic(X[,1],tree)*pic(X[,2],tree))/sqrt(sv1*sv2)
	if(is.null(lower)) lower<-c(log(0.1*c(sv1,sv2)),-1)
	if(is.null(upper)) upper<-c(log(10*c(sv1,sv2)),1)
	
	object<-list()
	class(object)<-"evolvcv.lite"
	
	loglik_model1<-if(hasArg(err_vcv)) 
		evol.vcv(as.phylo(tree),X,err_vcv=err_vcv)$logL1 else 
		evol.vcv(as.phylo(tree),X)$logL1

	# now optimize models
	msg<-function(text){
		cat(text)
		flush.console()
	}
	if("1"%in%models){
		msg("Fitting model 1: common rates, common correlation...\n")
		res1<-list()
		res1$convergence<-99
		class(res1)<-"try-error"
		iter<-0
		best<--loglik_model1+tol
		LOWER<-lower
		UPPER<-upper
		while((inherits(res1,"try-error")||res1$convergence!=0)&&iter<try.iter){
			init<-c(rnorm(n=2)*log(c(sv1,sv2)),runif(n=1,-1,1))
			res1<-try(optim(init,lik1,C=C,D=D,y=y,E=E,
				method="L-BFGS-B",lower=LOWER+tol,
				upper=UPPER-tol))
			if(inherits(res1,"try-error")){
				res1<-list()
				res1$convergence<-99
				class(res1)<-"try-error"
				res1$value<-best+tol
			}
			if(res1$value>best){
				if(res1$convergence==0) res1$convergence<-99
			} else best<-res1$value
			iter<-iter+1
		}
		if(inherits(res1,"try-error")){
			m1<-list(description="common rates, common correlation",
				R=matrix(NA,p,p),logLik=NA,k=5,AIC=NA)
		} else {
			if(any(res1$par<(LOWER+2*tol))||any(res1$par>(UPPER-2*tol)))
				res1$convergence<-77
			res1$par[1:2]<-exp(res1$par[1:2])
			m1<-list(description="common rates, common correlation",
				R=matrix(c(res1$par[1],res1$par[3]*sqrt(res1$par[1]*res1$par[2]),
				res1$par[3]*sqrt(res1$par[1]*res1$par[2]),res1$par[2]),2,2),
				logLik=-res1$value,
				convergence=res1$convergence,
				k=length(res1$par)+2,
				AIC=2*(length(res1$par)+2)+2*res1$value)
		}
		object$model1=m1
		msg(paste("Best log(L) from model 1: ",round(m1$logLik,4),".\n",sep=""))
	}  else m1<-NULL
	makeError<-function(description,names,p,k){
		list(description=description,
			R=setNames(replicate(p,matrix(NA,2,2),simplify=FALSE),names),
			logLik=NA,
			convergence=99,
			k=k,
			AIC=NA)
	}
	if("2"%in%models){
		msg("Fitting model 2: different rates, common correlation...\n")
		res2<-list()
		res2$convergence<-99
		class(res2)<-"try-error"
		iter<-0
		best<-if(is.null(m1)) -loglik_model1+tol else 
			-m1$logLik
		LOWER<-c(rep(lower[1:2],p),lower[3])
		UPPER<-c(rep(upper[1:2],p),upper[3])
		while((inherits(res2,"try-error")||res2$convergence!=0)&&iter<try.iter){
			init<-c(rnorm(n=2*p)*log(c(sv1,sv2)),runif(n=1,-1,1))
			res2<-try(optim(init,lik2,C=mC,
				D=D,y=y,E=E,method="L-BFGS-B",
				lower=tol+LOWER,upper=UPPER-tol))
			if(inherits(res2,"try-error")){
				res2<-list()
				res2$convergence<-99
				class(res2)<-"try-error"
				res2$value<-best+tol
			}
			if(res2$value>best){
				if(res2$convergence==0) res2$convergence<-99
			} else best<-res2$value
			iter<-iter+1
		}
		if(inherits(res2,"try-error")){
			m2<-makeError("different rates, common correlation",
				colnames(tree$mapped.edge),
				ncol(tree$mapped.edge),
				2*ncol(tree$mapped.edge)+1)
		} else {
			if(any(res2$par<(LOWER+2*tol))||any(res2$par>(UPPER-2*tol)))
				res2$convergence<-77
			R<-list()
			res2$par[1:(2*p)]<-exp(res2$par[1:(2*p)])
			for(i in 1:p) R[[i]]<-matrix(c(res2$par[2*(i-1)+1],
				rep(res2$par[2*p+1]*sqrt(res2$par[2*(i-1)+1]*res2$par[2*(i-1)+2]),2),
				res2$par[2*(i-1)+2]),2,2)
			names(R)<-colnames(tree$mapped.edge)
			m2<-list(description="different rates, common correlation",
				R=R,
				logLik=-res2$value,
				convergence=res2$convergence,
				k=length(res2$par)+2,
				AIC=2*(length(res2$par)+2)+2*res2$value)
		}
		object$model2<-m2
		msg(paste("Best log(L) from model 2: ",round(m2$logLik,4),".\n",sep=""))
	}
	if("2b"%in%models){
		msg("Fitting model 2b: different rates (trait 1), common correlation...\n")
		res2b<-list()
		res2b$convergence<-99
		class(res2b)<-"try-error"
		iter<-0
		best<-if(is.null(m1)) -loglik_model1+tol else 
			(-m1$logLik)
		LOWER<-c(rep(lower[1],p),lower[2],lower[3])
		UPPER<-c(rep(upper[1],p),upper[2],upper[3])
		while((inherits(res2b,"try-error")||res2b$convergence!=0)&&iter<try.iter){
			init<-c(rnorm(n=p+1)*log(c(rep(sv1,p),sv2)),runif(n=1,-1,1))
			res2b<-try(optim(init,lik2b,C=mC,
				D=D,y=y,E=E,method="L-BFGS-B",
				lower=LOWER+tol,upper=UPPER-tol))
			if(inherits(res2b,"try-error")){
				res2b<-list()
				res2b$convergence<-99
				class(res2b)<-"try-error"
				res2b$value<-best+tol
			}
			if(res2b$value>best){
				if(res2b$convergence==0) res2b$convergence<-99
			} else best<-res2b$value
			iter<-iter+1
		}
		if(inherits(res2b,"try-error")){
			m2b<-makeError("different rates (trait 1), common correlation",
				colnames(tree$mapped.edge),
				ncol(tree$mapped.edge),
				2*ncol(tree$mapped.edge)+2)
		} else {
			if(any(res2b$par<(LOWER+2*tol))||any(res2b$par>(UPPER-2*tol)))
				res2b$convergence<-77
			R<-list()
			res2b$par[1:(p+1)]<-exp(res2b$par[1:(p+1)])
			for(i in 1:p) R[[i]]<-matrix(c(res2b$par[i],
				rep(res2b$par[p+2]*sqrt(res2b$par[i]*res2b$par[p+1]),2),
				res2b$par[p+1]),2,2)
			names(R)<-colnames(tree$mapped.edge)
			m2b<-list(description="different rates (trait 1), common correlation",
				R=R,
				logLik=-res2b$value,
				convergence=res2b$convergence,
				k=length(res2b$par)+2,
				AIC=2*(length(res2b$par)+2)+2*res2b$value)
		}
		object$model2b<-m2b
		msg(paste("Best log(L) from model 2b: ",round(m2b$logLik,4),".\n",sep=""))
	}
	if("2c"%in%models){
		msg("Fitting model 2c: different rates (trait 2), common correlation...\n")
		res2c<-list()
		res2c$convergence<-99
		class(res2c)<-"try-error"
		iter<-0
		best<-if(is.null(m1)) -loglik_model1+tol else 
			(-m1$logLik)
		LOWER<-c(lower[1],rep(lower[2],p),lower[3])
		UPPER<-c(upper[1],rep(upper[2],p),upper[3])
		while((inherits(res2c,"try-error")||res2c$convergence!=0)&&iter<try.iter){
			init<-c(rnorm(n=p+1)*log(c(sv1,rep(sv2,p))),runif(n=1,-1,1))
			res2c<-try(optim(init,lik2c,C=mC,
				D=D,y=y,E=E,method="L-BFGS-B",
				lower=LOWER+tol,upper=UPPER-tol))
			if(inherits(res2c,"try-error")){
				res2c<-list()
				res2c$convergence<-99
				class(res2c)<-"try-error"
				res2c$value<-best+tol
			}
			if(res2c$value>best){
				if(res2c$convergence==0) res2c$convergence<-99
			} else best<-res2c$value
			iter<-iter+1
		}
		if(inherits(res2c,"try-error")){
			m2c<-makeError("different rates (trait 2), common correlation",
				colnames(tree$mapped.edge),
				ncol(tree$mapped.edge),
				2*ncol(tree$mapped.edge)+2)
		} else {
			if(any(res2c$par<(LOWER+2*tol))||any(res2c$par>(UPPER-2*tol)))
				res2c$convergence<-77
			R<-list()
			res2c$par[1:(p+1)]<-exp(res2c$par[1:(p+1)])
			for(i in 1:p) R[[i]]<-matrix(c(res2c$par[1],
				rep(res2c$par[p+2]*sqrt(res2c$par[1]*res2c$par[i+1]),2),
				res2c$par[i+1]),2,2)
			names(R)<-colnames(tree$mapped.edge)
			m2c<-list(description="different rates (trait 2), common correlation",
				R=R,
				logLik=-res2c$value,
				convergence=res2c$convergence,
				k=length(res2c$par)+2,
				AIC=2*(length(res2c$par)+2)+2*res2c$value)
		}
		object$model2c<-m2c
		msg(paste("Best log(L) from model 2c: ",round(m2c$logLik,4),".\n",sep=""))
	}
	if("3"%in%models){
		msg("Fitting model 3: common rates, different correlation...\n")
		res3<-list()
		res3$convergence<-99
		class(res3)<-"try-error"
		iter<-0
		best<-if(is.null(m1)) -loglik_model1+tol else 
			(-m1$logLik)
		LOWER<-c(lower[1],lower[2],rep(lower[3],p))
		UPPER<-c(upper[1],upper[2],rep(upper[3],p))
		while((inherits(res3,"try-error")||res3$convergence!=0)&&iter<try.iter){
			init<-c(rnorm(n=2)*log(c(sv1,sv2)),runif(n=p,-1,1))
			res3<-try(optim(init,lik3,C=mC,D=D,
				y=y,E=E,method="L-BFGS-B",lower=tol+LOWER,
				upper=UPPER-tol))
			if(inherits(res3,"try-error")){
				res3<-list()
				res3$convergence<-99
				class(res3)<-"try-error"
				res3$value<-best+tol
			}
			if(res3$value>best){
				if(res3$convergence==0) res3$convergence<-99
			} else best<-res3$value
			iter<-iter+1
		} 
		if(inherits(res3,"try-error")){
			m3<-makeError("common rates, different correlation",
				colnames(tree$mapped.edge),
				ncol(tree$mapped.edge),
				ncol(tree$mapped.edge)+2)
		} else {
			if(any(res3$par<(LOWER+2*tol))||any(res3$par>(UPPER-2*tol)))
				res3$convergence<-77
			R<-list()
			res3$par[1:2]<-exp(res3$par[1:2])
			for(i in 1:p) R[[i]]<-matrix(c(res3$par[1],
				rep(res3$par[2+i]*sqrt(res3$par[1]*res3$par[2]),2),res3$par[2]),2,2)
			names(R)<-colnames(tree$mapped.edge)
			m3<-list(description="common rates, different correlation",
				R=R,
				logLik=-res3$value,
				convergence=res3$convergence,
				k=length(res3$par)+2,
				AIC=2*(length(res3$par)+2)+2*res3$value)
		}
		object$model3<-m3
		msg(paste("Best log(L) from model 3: ",round(m3$logLik,4),".\n",sep=""))
	}
	if("3b"%in%models){
		msg("Fitting model 3b: different rates (trait 1), different correlation...\n")
		res3b<-list()
		res3b$convergence<-99
		class(res3b)<-"try-error"
		iter<-0
		best<-if(is.null(m1)) -loglik_model1+tol else 
			(-m1$logLik)
		LOWER<-c(rep(lower[1],p),lower[2],rep(lower[3],p))
		UPPER<-c(rep(upper[1],p),upper[2],rep(upper[3],p))
		while((inherits(res3b,"try-error")||res3b$convergence!=0)&&iter<try.iter){
			init<-c(rnorm(n=p+1)*log(c(rep(sv1,p),sv2)),runif(n=p,-1,1))
			res3b<-try(optim(init,lik3b,C=mC,
				D=D,y=y,E=E,method="L-BFGS-B",
				lower=LOWER+tol,upper=UPPER-tol))
			if(inherits(res3b,"try-error")){
				res3b<-list()
				res3b$convergence<-99
				class(res3b)<-"try-error"
				res3b$value<-best+tol
			}
			if(res3b$value>best){
				if(res3b$convergence==0) res3b$convergence<-99
			} else best<-res3b$value
			iter<-iter+1
		}
		if(inherits(res3b,"try-error")){
			m3b<-makeError("different rates (trait 1), different correlation",
				colnames(tree$mapped.edge),
				ncol(tree$mapped.edge),
				2*ncol(tree$mapped.edge)+1)
		} else {
			if(any(res3b$par<(LOWER+2*tol))||any(res3b$par>(UPPER-2*tol)))
				res3b$convergence<-77
			R<-list()
			res3b$par[1:(p+1)]<-exp(res3b$par[1:(p+1)])
			for(i in 1:p) R[[i]]<-matrix(c(res3b$par[i],
				rep(res3b$par[p+1+i]*sqrt(res3b$par[i]*res3b$par[p+1]),2),
				res3b$par[p+1]),2,2)
			names(R)<-colnames(tree$mapped.edge)
			m3b<-list(description="different rates (trait 1), different correlation",
				R=R,
				logLik=-res3b$value,
				convergence=res3b$convergence,
				k=length(res3b$par)+2,
				AIC=2*(length(res3b$par)+2)+2*res3b$value)
		}
		object$model3b<-m3b
		msg(paste("Best log(L) from model 3b: ",round(m3b$logLik,4),".\n",sep=""))
	}
	if("3c"%in%models){
		msg("Fitting model 3c: different rates (trait 2), different correlation...\n")
		res3c<-list()
		res3c$convergence<-99
		class(res3c)<-"try-error"
		iter<-0
		best<-if(is.null(m1)) -loglik_model1+tol else 
			(-m1$logLik)
		LOWER<-c(lower[1],rep(lower[2],p),rep(lower[3],p))
		UPPER<-c(upper[1],rep(upper[2],p),rep(upper[3],p))
		while((inherits(res3c,"try-error")||res3c$convergence!=0)&&iter<try.iter){
			init<-c(rnorm(n=p+1)*log(c(sv1,rep(sv2,p))),runif(n=p,-1,1))
			res3c<-try(optim(init,lik3c,C=mC,
				D=D,y=y,E=E,method="L-BFGS-B",
				lower=LOWER+tol,upper=UPPER-tol))
			if(inherits(res3c,"try-error")){
				res3c<-list()
				res3c$convergence<-99
				class(res3c)<-"try-error"
				res3c$value<-best+tol
			}
			if(res3c$value>best){
				if(res3c$convergence==0) res3c$convergence<-99
			} else best<-res3c$value
			iter<-iter+1
		}
		if(inherits(res3c,"try-error")){
			m3c<-makeError("different rates (trait 2), different correlation",
				colnames(tree$mapped.edge),
				ncol(tree$mapped.edge),
				2*ncol(tree$mapped.edge)+1)
		} else {
			if(any(res3c$par<(LOWER+2*tol))||any(res3c$par>(UPPER-2*tol)))
				res3c$convergence<-77
			R<-list()
			res3c$par[1:(p+1)]<-exp(res3c$par[1:(p+1)])
			for(i in 1:p) R[[i]]<-matrix(c(res3c$par[1],
				rep(res3c$par[p+1+i]*sqrt(res3c$par[1]*res3c$par[1+i]),2),
				res3c$par[1+i]),2,2)
			names(R)<-colnames(tree$mapped.edge)
			m3c<-list(description="different rates (trait 2), different correlation",
				R=R,
				logLik=-res3c$value,
				convergence=res3c$convergence,
				k=length(res3c$par)+2,
				AIC=2*(length(res3c$par)+2)+2*res3c$value)
		}
		object$model3c<-m3c
		msg(paste("Best log(L) from model 3c: ",round(m3c$logLik,4),".\n",sep=""))
	}
	if("4"%in%models){
		msg("Fitting model 4: no common structure...\n")
		res4<-list()
		res4$convergence<-99
		class(res4)<-"try-error"
		iter<-0
		best<-if(is.null(m1)) -loglik_model1+tol else 
			(-m1$logLik)
		LOWER<-c(rep(lower[1:2],p),rep(lower[3],p))
		UPPER<-c(rep(upper[1:2],p),rep(upper[3],p))
		while((inherits(res4,"try-error")||res4$convergence!=0)&&iter<try.iter){
			init<-c(rnorm(n=2*p)*log(rep(c(sv1,sv2),p)),runif(n=p,-1,1))
			res4<-try(optim(init,lik4,C=mC,
				D=D,y=y,E=E,method="L-BFGS-B",lower=LOWER+tol,
				upper=UPPER-tol))
			if(inherits(res4,"try-error")){
				res4<-list()
				res4$convergence<-99
				class(res4)<-"try-error"
				res4$value<-best+tol
			}
			if(res4$value>best){
				if(res4$convergence==0) res4$convergence<-99
			} else best<-res4$value
			iter<-iter+1
		}
		if(inherits(res4,"try-error")){
			m4<-makeError("no common structure",
				colnames(tree$mapped.edge),
				ncol(tree$mapped.edge),
				3*ncol(tree$mapped.edge)+2)
		} else {
			if(any(res4$par<(LOWER+2*tol))||any(res4$par>(UPPER-2*tol)))
				res4$convergence<-77
			R<-list()
			res4$par[1:(2*p)]<-exp(res4$par[1:(2*p)])
			for(i in 1:p) R[[i]]<-matrix(c(res4$par[2*(i-1)+1],
				rep(res4$par[2*p+i]*sqrt(res4$par[2*(i-1)+1]*res4$par[2*(i-1)+2]),
				2),res4$par[2*(i-1)+2]),2,2)
			names(R)<-colnames(tree$mapped.edge)
			m4<-list(description="no common structure",
				R=R,
				logLik=-res4$value,
				convergence=res4$convergence,
				k=length(res4$par)+2,
				AIC=2*(length(res4$par)+2)+2*res4$value)
		}
		object$model4<-m4
		msg(paste("Best log(L) from model 4: ",round(m4$logLik,4),".\n",sep=""))
	}
	object
}

## S3 print method for object of class "evolvcv.lite"
## written by Liam J. Revell 2017, 2019, 2021
print.evolvcv.lite<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	
	if(!is.null(x$model1)){
		nn<-paste("R[",t(sapply(1:ncol(x$model1$R),paste,
			1:ncol(x$model1$R),sep=","))[upper.tri(x$model1$R,
			diag=TRUE)],"]",sep="")
		## MODEL 1	
		m1<-lapply(x$model1,function(a,b) if(is.numeric(a)) round(a,b) else a,b=digits)
		cat(paste("Model 1:",x$model1$description,"\n"))
		cat(paste("\t",paste(nn,collapse="\t"),"\tk\tlog(L)\tAIC","\n",sep=""))
		cat(paste("fitted",paste(m1$R[upper.tri(m1$R,diag=TRUE)],collapse="\t"),m1$k,
			m1$logLik,m1$AIC,"\n",sep="\t"))
		if(m1$convergence==0) cat("\n(R thinks it has found the ML solution for model 1.)\n\n")
		else if(m1$convergence==77) cat("\n(Model 1 optimization may be at bounds.)\n\n")
		else cat("\n(Model 1 optimization may not have converged.)\n\n")
		ii<-which(names(x)=="model1")
		x<-x[-ii]
	}
	
	if(length(x)>0){
		nn<-paste("R[",t(sapply(1:ncol(x[[1]]$R[[1]]),paste,
			1:ncol(x[[1]]$R[[1]]),sep=","))[upper.tri(x[[1]]$R[[1]],
			diag=TRUE)],"]",sep="")
		for(i in 1:length(x)){
			m<-lapply(x[[i]],function(a,b) if(is.numeric(a)) round(a,b) else a,b=digits)
			m$R<-lapply(m$R,function(a,b) if(is.numeric(a)) round(a,b) else a,b=digits)
			model<-strsplit(names(x)[i],"model")[[1]][2]
			cat(paste("Model ",model,": ",m$description,"\n",sep=""))
			cat(paste("\t",paste(nn,collapse="\t"),"\tk\tlog(L)\tAIC","\n",sep=""))
			for(j in 1:length(m$R)){
				if(j==1) cat(paste(names(m$R)[j],paste(m$R[[j]][upper.tri(m$R[[j]],diag=TRUE)],
					collapse="\t"),m$k,m$logLik,m$AIC,"\n",sep="\t"))
				else cat(paste(names(m$R)[j],paste(m$R[[j]][upper.tri(m$R[[j]],diag=TRUE)],
					collapse="\t"),"\n",sep="\t"))
			}
			if(m$convergence==0) 
				cat(paste("\n(R thinks it has found the ML solution for model ",model,".)\n\n",sep=""))
			else if(m$convergence==77)
				cat(paste("\n(Model ",model," optimization may be at bounds.)\n\n",sep=""))
			else cat(paste("\n(Model ",model," optimization may not have converged.)\n\n",sep=""))
		}
	}
}
