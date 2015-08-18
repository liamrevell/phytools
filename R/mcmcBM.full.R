# function
# written by Liam J. Revell 2011

mcmcBM.full<-function(tree,x,ngen=10000,control=list()){

	# starting values (for now)
	n<-length(tree$tip)
	temp<-aggregate(x,list(species=as.factor(names(x))),mean)
	xbar<-temp[,2]; names(xbar)<-temp[,1]; xbar<-xbar[tree$tip.label]
	sig2<-mean(pic(xbar,tree)^2)
	a<-mean(xbar)
	v<-rep(mean(aggregate(x,list(species=as.factor(names(x))),var)[,2],na.rm=T),n)
	names(v)<-names(xbar)
	prop<-c(0.01*sig2,0.01*sig2,rep(0.01*sig2*max(vcv(tree)),n),0.01*v)
	pr.mean<-c(1000,rep(0,n+1),rep(1000,n))
	pr.var<-c(pr.mean[1]^2,rep(1000,n+1),pr.mean[n+2+1:n]^2)

	# populate control list
	con=list(sig2=sig2,a=a,xbar=xbar,v=v,pr.mean=pr.mean,pr.var=pr.var,prop=prop,sample=100)
	names(con$v)<-gsub("\\)","",gsub("var\\(","",names(con$v)))
	con[(namc<-names(control))]<-control
	con<-con[!sapply(con,is.null)]
	# print control parameters to screen
	message("Control parameters (set by user or default):"); str(con)
	
	# function returns the log-likelihood
	likelihood<-function(C,invC,detC,x,sig2,a,xbar,v){
		z<-xbar-a
		logLik<--z%*%invC%*%z/(2*sig2)-nrow(C)*log(2*pi)/2-nrow(C)*log(sig2)/2-detC/2+sum(dnorm(x,xbar[names(x)],sd=sqrt(v[names(x)]),log=T))
		return(logLik)
	}

	# function returns the log prior probability
	log.prior<-function(sig2,a,xbar,v){
		pp<-dexp(sig2,rate=1/con$pr.mean[1],log=T)+sum(dnorm(c(a,xbar),mean=con$pr.mean[1+1:(n+1)],sd=sqrt(con$pr.var[1+1:(n+1)]),log=T))+sum(dexp(v,rate=1/con$pr.mean[n+2+1:n],log=T))
		return(pp)
	}

	# compute C
	C<-vcv.phylo(tree)
	invC<-solve(C)
	detC<-determinant(C,logarithm=TRUE)$modulus[1]

	# now set starting values for MCMC
	sig2<-con$sig2; a<-con$a; xbar<-con$xbar; v<-con$v
	L<-likelihood(C,invC,detC,x,sig2,a,xbar,v)
	Pr<-log.prior(sig2,a,xbar,v)

	# store
	X<-matrix(NA,ngen/con$sample+1,2*n+4,dimnames=list(NULL,c("gen","sig2","a",tree$tip.label,paste("var(",tree$tip.label,")",sep=""),"logLik")))
	X[1,]<-c(0,sig2,a,xbar,v,L)

	message("Starting MCMC...")

	# start MCMC
	for(i in 1:ngen){
		j<-(i-1)%%(2*n+2)
		if(j==0){
			# update sig2
			sig2.prime<-sig2+rnorm(n=1,sd=sqrt(con$prop[j+1]))
			if(sig2.prime<0) sig2.prime<--sig2.prime
			L.prime<-likelihood(C,invC,detC,x,sig2.prime,a,xbar,v)
			Pr.prime<-log.prior(sig2.prime,a,xbar,v)
			post.odds<-min(1,exp(Pr.prime+L.prime-Pr-L),na.rm=T)
			if(post.odds>runif(n=1)){
				if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2.prime,a,xbar,v,L.prime)
				sig2<-sig2.prime
				L<-L.prime
				Pr<-Pr.prime
			} else if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,xbar,v,L)
		} else if(j==1){
			# update a
			a.prime<-a+rnorm(n=1,sd=sqrt(con$prop[j+1]))
			L.prime<-likelihood(C,invC,detC,x,sig2,a.prime,xbar,v)
			Pr.prime<-log.prior(sig2,a.prime,xbar,v)
			post.odds<-min(1,exp(Pr.prime+L.prime-Pr-L),na.rm=T)
			if(post.odds>runif(n=1)){ 
				if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a.prime,xbar,v,L.prime)
				a<-a.prime
				L<-L.prime
				Pr<-Pr.prime
			} else if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,xbar,v,L)
		} else if(j>1&&j<=(n+1)) {
			k<-j-1 # update tip mean k
			xbar.prime<-xbar
			xbar.prime[k]<-xbar[k]+rnorm(n=1,sd=sqrt(con$prop[j+1]))
			L.prime<-likelihood(C,invC,detC,x,sig2,a,xbar.prime,v)
			Pr.prime<-log.prior(sig2,a,xbar.prime,v)
			post.odds<-min(1,exp(Pr.prime+L.prime-Pr-L),na.rm=T)
			if(post.odds>runif(n=1)){ 
				if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,xbar.prime,v,L.prime)
				xbar<-xbar.prime
				L<-L.prime
				Pr<-Pr.prime
			} else if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,xbar,v,L)
		} else if(j>(n+1)){
			k<-j-n-1 # update var
			v.prime<-v
			v.prime[k]<-v[k]+rnorm(n=1,sd=sqrt(con$prop[j+1]))
			if(v.prime[k]<0) v.prime[k]<--v.prime[k]
			L.prime<-likelihood(C,invC,detC,x,sig2,a,xbar,v.prime)
			Pr.prime<-log.prior(sig2,a,xbar,v.prime)
			post.odds<-min(1,exp(Pr.prime+L.prime-Pr-L),na.rm=T)
			if(post.odds>runif(n=1)){ 
				if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,xbar,v.prime,L.prime)
				v<-v.prime
				L<-L.prime
				Pr<-Pr.prime
			} else if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,xbar,v,L)	
		}
	}

	# done MCMC
	message("Done MCMC.")
	return(X)
}
