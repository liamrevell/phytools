## function for computing phylogenetic signal by the lambda (Pagel 1999) of K (Blomberg et al. 2003) methods
## written by Liam J. Revell 2011/2012, 2019

phylosig<-function(tree,x,method="K",test=FALSE,nsim=1000,se=NULL,start=NULL,control=list()){
	# some minor error checking
	if(!inherits(tree,"phylo")) 
		stop("tree should be an object of class \"phylo\".")
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
			K<-(t(x-a)%*%(x-a)/(t(x-a)%*%invC%*%(x-a)))/((sum(diag(C))-
				n/sum(invC))/(n-1)) # calculate K
			if(!test){
				object<-as.numeric(K)
			} else {
				P=0.0
				simX<-x
				simK<-vector()
				for(i in 1:nsim){
					a<-sum(invC%*%simX)/sum(invC)
					simK[i]<-(t(simX-a)%*%(simX-a)/(t(simX-a)%*%invC%*%
						(simX-a)))/((sum(diag(C))-n/sum(invC))/(n-1))
					## calculate P-value for randomization test
					if(simK[i]>=K) P<-P+1/nsim 
					simX<-sample(simX) # randomize x
				}
				object<-list(K=as.numeric(K),P=P,sim.K=simK)
			}
		} else {
			likelihoodK<-function(theta,C,M,y){
				Ce<-theta*C+M
				invCe<-solve(Ce)
				a<-as.numeric(sum(invCe%*%y)/sum(invCe))
				logL<--t(y-a)%*%invCe%*%(y-a)/2-n*log(2*pi)/2-
					determinant(Ce,logarithm=TRUE)$modulus/2
				logL[1,1]
			}
			M<-M[rownames(C),colnames(C)]
			invC<-solve(C)
			maxSig2<-as.numeric(t(x-as.numeric(sum(invC%*%x)/
				sum(invC)))%*%invC%*%(x-as.numeric(sum(invC%*%x)/
				sum(invC)))/n)
			res<-optimize(f=likelihoodK,interval=c(0,maxSig2),y=x,
				C=C,M=M,maximum=TRUE) # optimize sig2
			sig2<-res$maximum*n/(n-1)
			Ce<-sig2*C+M
			invCe<-solve(Ce)
			a<-as.numeric(sum(invCe%*%x)/sum(invCe))
			K<-(t(x-a)%*%(x-a)/(t(x-a)%*%invCe%*%(x-a)))/((sum(diag(Ce))-
				n/sum(invCe))/(n-1)) # calculate K
			if(!test){
				object<-list(K=as.numeric(K),sig2=as.numeric(sig2),
					logL=res$objective,
					lik=function(sig2) likelihoodK(sig2,C=C,M=M,y=x))
			} else {
				P=0.0
				simX<-x
				simK<-vector()
				for(i in 1:nsim){
					maxSig2<-as.numeric(t(simX-as.numeric(sum(invC%*%
						simX)/sum(invC)))%*%invC%*%(simX-as.numeric(sum(invC%*%
						simX)/sum(invC)))/n)
					simRes<-optimize(f=likelihoodK,interval=c(0,maxSig2),y=simX,
						C=C,M=M,maximum=TRUE) # optimize sig2
					simSig2<-simRes$maximum*n/(n-1)
					Ce<-simSig2*C+M
					invCe<-solve(Ce)
					a<-as.numeric(sum(invCe%*%simX)/sum(invCe))
					# calculate K
					simK[i]<-(t(simX-a)%*%(simX-a)/(t(simX-a)%*%invCe%*%
						(simX-a)))/((sum(diag(Ce))-n/sum(invCe))/(n-1)) 
					# calculate P-value for randomization test
					if(simK[i]>=K) P<-P+1/nsim 
					o<-sample(1:n)
					simX<-x[o]
					M<-diag(se[o]^2) # randomize x & errors
				}
				object<-list(K=as.numeric(K),P=P,sim.K=simK,
					sig2=as.numeric(sig2),logL=res$objective,
					lik=function(sig2) likelihoodK(sig2,C=C,M=M,y=x))
			}
		}
	} else if(method=="lambda"){
		# function to compute C with lambda
		lambda.transform<-function(C,lambda){
			dC<-diag(diag(C))
			C<-lambda*(C-dC)+dC
			C
		}
		# likelihood function
		likelihoodLambda<-function(theta,C,y){
			Cl<-lambda.transform(C,theta)
			invCl<-solve(Cl)
			n<-nrow(Cl)
			y<-y[rownames(Cl)]
			a<-as.numeric(sum(invCl%*%y)/sum(invCl))
			sig2<-as.numeric(t(y-a)%*%invCl%*%(y-a)/n)
			logL<--t(y-a)%*%(1/sig2*invCl)%*%(y-a)/2-n*log(2*pi)/2-
				determinant(sig2*Cl,logarithm=TRUE)$modulus/2
			logL[1,1]
		}
		# likelihood function with error
		likelihoodLambda.me<-function(theta,C,y,M){
			Cl<-theta[1]*lambda.transform(C,theta[2])
			V<-Cl+M
			invV<-solve(V)
			n<-nrow(Cl)
			y<-y[rownames(Cl)]
			a<-as.numeric(sum(invV%*%y)/sum(invV))
			logL<--t(y-a)%*%invV%*%(y-a)/2-n*log(2*pi)/2-
				determinant(V,logarithm=TRUE)$modulus/2
			logL[1,1]
		}
		C<-vcv.phylo(tree)
		x<-x[rownames(C)]
		maxlam<-maxLambda(tree)
		if(!me){
			res<-optimize(f=likelihoodLambda,interval=c(0,maxlam),
				y=x,C=C,maximum=TRUE) # optimize lambda
			if(!test){
				object<-list(lambda=res$maximum,logL=res$objective,
					lik=function(lambda) likelihoodLambda(lambda,
					C=C,y=x))
			} else {
				# compute likelihood of lambda=0
				logL0<-likelihoodLambda(theta=0,C=C,y=x) 
				P<-as.numeric(pchisq(2*(res$objective-logL0),df=1,
					lower.tail=FALSE)) # P-value
				object<-list(lambda=res$maximum,logL=res$objective,
					logL0=logL0,P=P,lik=function(lambda) 
					likelihoodLambda(lambda,C=C,y=x))
			}
		} else {
			control$fnscale=-1
			M<-M[rownames(C),colnames(C)]
			if(is.null(start)) s<-c(0.02*runif(n=1)*mean(pic(x,
				multi2di(tree))^2),runif(n=1))
			else s<-start
			res<-optim(s,likelihoodLambda.me,C=C,y=x,M=M,
				method="L-BFGS-B",lower=c(0,0),upper=c(Inf,maxlam),
				control=control)
			if(!test){
				object<-list(lambda=res$par[2],sig2=res$par[1],
					logL=res$value,convergence=res$convergence,
					message=res$message,lik=function(lambda) 
					likelihoodLambda(lambda,C=res$par[1]*C+M,y=x))
			} else {
				res0<-optim(c(s[1],0),likelihoodLambda.me,C=C,
					y=x,M=M,method="L-BFGS-B",lower=c(0,0),
					upper=c(Inf,1e-10),control=control)
				P<-as.numeric(pchisq(2*(res$value-res0$value),df=1,
					lower.tail=FALSE))
				object<-list(lambda=res$par[2],sig2=res$par[1],
					logL=res$value,convergence=res$convergence,
					message=res$message,logL0=res0$value,P=P,
					lik=function(lambda) likelihoodLambda.me(c(res$par[1],
					lambda),C=C,M=M,y=x))
			}
		}
	} else
		stop(paste("do not recognize method = \"",method,
			"\"; methods are \"K\" and \"lambda\"",sep=""))
	attr(object,"class")<-"phylosig"
	attr(object,"method")<-method
	attr(object,"test")<-test
	attr(object,"se")<-!is.null(se)
	object
}

print.phylosig<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-6
	cat("\n")
	if(attr(x,"method")=="K"){
		if(attr(x,"test")||attr(x,"se")){
			cat(paste("Phylogenetic signal K :",signif(x$K,
				digits),"\n"))
			if(attr(x,"se")){
				cat(paste("MLE(sig2) :",signif(x$sig2,digits),
					"\n"))
				cat(paste("logL(sig2) :",signif(x$logL,digits),
					"\n"))
			}
			if(attr(x,"test"))
				cat(paste("P-value (based on",length(x$sim.K),
					"randomizations) :",signif(x$P,digits),"\n"))
		} else cat(paste("Phylogenetic signal K :",signif(x[1],
			digits),"\n"))
	} else if(attr(x,"method")=="lambda"){
		cat(paste("Phylogenetic signal lambda :",signif(x$lambda,
			digits),"\n"))
		cat(paste("logL(lambda) :",signif(x$logL,digits),"\n"))
		if(attr(x,"se")) cat(paste("MLE(sig2) :",signif(x$sig2,
			digits),"\n"))
		if(attr(x,"test")){
			cat(paste("LR(lambda=0) :",signif(2*(x$logL-x$logL0),
				digits),"\n"))
			cat(paste("P-value (based on LR test) :",signif(x$P,
				digits),"\n"))
		}
	}
	cat("\n")			
}

plot.phylosig<-function(x,...){
	if(hasArg(what)) what<-list(...)$what
	else what<-if(attr(x,"method")=="lambda") "lambda" else 
		if(attr(x,"method")=="K"&&attr(x,"test")) "K" else "sig2"
	if(hasArg(res)) res<-list(...)$res
	else res<-100
	if(attr(x,"method")=="lambda"){
		lambda<-seq(0,max(c(1,x$lambda)),length.out=res)
		logL<-sapply(lambda,x$lik)
		plot(lambda,logL,xlab=expression(lambda),ylab="log(L)",
			type="l",bty="l")
		lines(rep(x$lambda,2),c(par()$usr[3],x$logL),lty="dotted")
		text(x=x$lambda+0.01*diff(par()$usr[1:2]),
			par()$usr[3]+0.5*diff(par()$usr[3:4]),
			expression(paste("MLE(",lambda,")")),srt=90,adj=c(0.5,1))
		if(attr(x,"test")){
			lines(c(0,x$lambda),rep(x$logL0,2),lty="dotted")
			text(x=0.5*x$lambda,x$logL0-0.01*diff(par()$usr[3:4]),
				expression(paste("logL(",lambda,"=0)")),adj=c(0.5,1))
		}
	} else if(attr(x,"method")=="K"){
		if(what=="K"){
			if(attr(x,"test")==FALSE)
				cat("Sorry. This is not a valid plotting option for your object.\n\n")
			else {
				hist(x$sim.K,breaks=min(c(max(12,round(length(x$sim.K)/10)),
					20)),bty="l",col="lightgrey",border="lightgrey",
					main="",xlab="K",ylab="null distribution of K")
				arrows(x0=x$K,y0=par()$usr[4],y1=0,length=0.12,
					col=make.transparent("blue",0.5),lwd=2)
				text(x$K,0.98*par()$usr[4],"observed value of K",
					pos=if(x$K>mean(x$sim.K)) 2 else 4)
			}
		} else if(what=="sig2"){
			if(attr(x,"se")==FALSE)
				cat("Sorry. This is not a valid plotting option for your object.\n\n")
			else {
				sig2<-seq(0.25*x$sig2,1.75*x$sig2,length.out=res)
				logL<-sapply(sig2,x$lik)
				plot(sig2,logL,xlab=expression(sigma^2),ylab="log(L)",
					type="l",bty="l")
				lines(rep(x$sig2,2),c(par()$usr[3],x$logL),lty="dotted")
				text(x=x$sig2+0.01*diff(par()$usr[1:2]),
					par()$usr[3]+0.5*diff(par()$usr[3:4]),
					expression(paste("MLE(",sigma^2,")")),srt=90,
					adj=c(0.5,1))
			}
		}	
	}
}