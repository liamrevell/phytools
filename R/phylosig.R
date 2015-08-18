# function for computing phylogenetic signal by the lambda (Pagel 1999) of K (Blomberg et al. 2003) methods
# written by Liam J. Revell 2011/2012

phylosig<-function(tree,x,method="K",test=FALSE,nsim=1000,se=NULL,start=NULL,control=list()){
	# some minor error checking
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	x<-matchDatatoTree(tree,x,"x")
	tree<-matchTreetoData(tree,x,"x")
	if(!is.null(se)){
		se<-matchDatatoTree(tree,se,"se")
		tree<-matchTreetoData(tree,se,"se")
		me=TRUE
		M<-diag(se^2)
		rownames(M)<-colnames(M)<-names(se)
	} else me=FALSE
	if(!is.null(start)&&!is.null(se)){
		if(start[1]<=0||start[2]<0||start[2]>maxLambda(tree)){
			message("some of the elements of 'start' are invalid, resetting to random")
			start<-NULL
		}
	}
	# done error handling
	if(method=="K"){
		C<-vcv.phylo(tree)
		x<-x[rownames(C)]
		n<-nrow(C)
		if(!me){
			invC<-solve(C)
			a<-sum(invC%*%x)/sum(invC)
			K<-(t(x-a)%*%(x-a)/(t(x-a)%*%invC%*%(x-a)))/((sum(diag(C))-n/sum(invC))/(n-1)) # calculate K
			if(!test)
				return(as.numeric(K)) # return K
			else {
				P=0.0
				simX<-x
				for(i in 1:nsim){
					a<-sum(invC%*%simX)/sum(invC)
					simK<-(t(simX-a)%*%(simX-a)/(t(simX-a)%*%invC%*%(simX-a)))/((sum(diag(C))-n/sum(invC))/(n-1))
					if(simK>=K) P<-P+1/nsim # calculate P-value for randomization test
					simX<-sample(simX) # randomize x
				}
				return(list(K=as.numeric(K),P=P)) # return K & P
			}
		} else {
			likelihoodK<-function(theta,C,M,y){
				Ce<-theta*C+M
				invCe<-solve(Ce)
				a<-as.numeric(sum(invCe%*%y)/sum(invCe))
				logL<--t(y-a)%*%invCe%*%(y-a)/2-n*log(2*pi)/2-determinant(Ce,logarithm=TRUE)$modulus/2
				return(logL)
			}
			M<-M[rownames(C),colnames(C)]
			invC<-solve(C)
			maxSig2<-as.numeric(t(x-as.numeric(sum(invC%*%x)/sum(invC)))%*%invC%*%(x-as.numeric(sum(invC%*%x)/sum(invC)))/n)
			res<-optimize(f=likelihoodK,interval=c(0,maxSig2),y=x,C=C,M=M,maximum=TRUE) # optimize sig2
			sig2<-res$maximum*n/(n-1)
			Ce<-sig2*C+M
			invCe<-solve(Ce)
			a<-as.numeric(sum(invCe%*%x)/sum(invCe))
			K<-(t(x-a)%*%(x-a)/(t(x-a)%*%invCe%*%(x-a)))/((sum(diag(Ce))-n/sum(invCe))/(n-1)) # calculate K
			if(!test)
				return(list(K=as.numeric(K),sig2=as.numeric(sig2),logL=res$objective[1,1]))
			else {
				P=0.0
				simX<-x
				for(i in 1:nsim){
					maxSig2<-as.numeric(t(simX-as.numeric(sum(invC%*%simX)/sum(invC)))%*%invC%*%(simX-as.numeric(sum(invC%*%simX)/sum(invC)))/n)
					simRes<-optimize(f=likelihoodK,interval=c(0,maxSig2),y=simX,C=C,M=M,maximum=TRUE) # optimize sig2
					simSig2<-simRes$maximum*n/(n-1)
					Ce<-simSig2*C+M
					invCe<-solve(Ce)
					a<-as.numeric(sum(invCe%*%simX)/sum(invCe))
					simK<-(t(simX-a)%*%(simX-a)/(t(simX-a)%*%invCe%*%(simX-a)))/((sum(diag(Ce))-n/sum(invCe))/(n-1)) # calculate K
					if(simK>=K) P<-P+1/nsim # calculate P-value for randomization test
					o<-sample(1:n)
					simX<-x[o]; M<-diag(se[o]^2) # randomize x & errors
				}
				return(list(K=as.numeric(K),P=P,sig2=as.numeric(sig2),logL=res$objective[1,1])) # return K & P
			}
		}
	} else if(method=="lambda"){
		# function to compute C with lambda
		lambda.transform<-function(C,lambda){
			dC<-diag(diag(C))
			C<-lambda*(C-dC)+dC
			return(C)
		}
		# likelihood function
		likelihoodLambda<-function(theta,C,y){
			Cl<-lambda.transform(C,theta)
			invCl<-solve(Cl)
			n<-nrow(Cl)
			y<-y[rownames(Cl)]
			a<-as.numeric(sum(invCl%*%y)/sum(invCl))
			sig2<-as.numeric(t(y-a)%*%invCl%*%(y-a)/n)
			logL<--t(y-a)%*%(1/sig2*invCl)%*%(y-a)/2-n*log(2*pi)/2-determinant(sig2*Cl,logarithm=TRUE)$modulus/2
			return(logL)
		}
		# likelihood function with error
		likelihoodLambda.me<-function(theta,C,y,M){
			Cl<-theta[1]*lambda.transform(C,theta[2])
			V<-Cl+M
			invV<-solve(V)
			n<-nrow(Cl)
			y<-y[rownames(Cl)]
			a<-as.numeric(sum(invV%*%y)/sum(invV))
			logL<--t(y-a)%*%invV%*%(y-a)/2-n*log(2*pi)/2-determinant(V,logarithm=TRUE)$modulus/2
			return(-logL)
		}
		C<-vcv.phylo(tree)
		x<-x[rownames(C)]
		maxlam<-maxLambda(tree)
		if(!me){
			res<-optimize(f=likelihoodLambda,interval=c(0,maxlam),y=x,C=C,maximum=TRUE) # optimize lambda
			if(!test)
				return(list(lambda=res$maximum,logL=res$objective[1,1])) # return lambda and log-likelihood
			else {
				logL0<-likelihoodLambda(theta=0,C=C,y=x) # compute likelihood of lambda=0
				P<-as.numeric(pchisq(2*(res$objective[1,1]-logL0),df=1,lower.tail=FALSE)) # P-value
				return(list(lambda=res$maximum,logL=res$objective[1,1],logL0=logL0[1,1],P=P)) # return lambda, logL, and P-value
			}
		} else {
			M<-M[rownames(C),colnames(C)]
			if(is.null(start)) s<-c(0.02*runif(n=1)*mean(pic(x,multi2di(tree))^2),runif(n=1))
			else s<-start
			res<-optim(s,likelihoodLambda.me,C=C,y=x,M=M,method="L-BFGS-B",lower=c(0,0),upper=c(Inf,maxlam),control=control)
			if(!test)
				return(list(lambda=res$par[2],sig2=res$par[1],logL=-res$value,convergence=res$convergence,message=res$message))
			else {
				res0<-optim(c(s[1],0),likelihoodLambda.me,C=C,y=x,M=M,method="L-BFGS-B",lower=c(0,0),upper=c(Inf,1e-10),control=control)
				P<-as.numeric(pchisq(2*(res0$value-res$value),df=1,lower.tail=FALSE))
				return(list(lambda=res$par[2],sig2=res$par[1],logL=-res$value,convergence=res$convergence,message=res$message,logL0=-res0$value,P=P))
			}
		}
	} else
		stop(paste("do not recognize method = \"",method,"\"; methods are \"K\" and \"lambda\"",sep=""))
}

