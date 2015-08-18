# function
# written by Liam J. Revell 2011

mcmcBM<-function(tree,x,ngen=10000,control=list()){

	# starting values (for now)
	n<-length(tree$tip)
	temp<-aggregate(x,list(species=as.factor(names(x))),mean)
	xbar<-temp[,2]; names(xbar)<-temp[,1]; xbar<-xbar[tree$tip.label]
	sig2<-mean(pic(xbar,tree)^2)
	a<-mean(xbar)
	intV<-mean(aggregate(x,list(species=as.factor(names(x))),var)[,2],na.rm=T)
	prop<-c(0.01*sig2,0.01*sig2,rep(0.01*sig2*max(vcv(tree)),n),0.01*intV)
	pr.mean<-c(1000,rep(0,n+1),1000)
	pr.var<-c(pr.mean[1]^2,rep(1000,n+1),pr.mean[length(pr.mean)]^2)

	# populate control list
	con=list(sig2=sig2,a=a,xbar=xbar,intV=intV,pr.mean=pr.mean,pr.var=pr.var,prop=prop,sample=100)
	con[(namc<-names(control))]<-control
	con<-con[!sapply(con,is.null)]
	# print control parameters to screen
	message("Control parameters (set by user or default):"); str(con)
	
	# function returns the log-likelihood
	likelihood<-function(C,invC,detC,x,sig2,a,xbar,intV){
		z<-xbar-a
		logLik<--z%*%invC%*%z/(2*sig2)-nrow(C)*log(2*pi)/2-nrow(C)*log(sig2)/2-detC/2+sum(dnorm(x,xbar[names(x)],sd=sqrt(intV),log=T))
		return(logLik)
	}

	# function returns the log prior probability
	log.prior<-function(sig2,a,xbar,intV){
		pp<-dexp(sig2,rate=1/con$pr.mean[1],log=T)+sum(dnorm(c(a,xbar),mean=con$pr.mean[1+1:(n+1)],sd=sqrt(con$pr.var[1+1:(n+1)]),log=T))+dexp(intV,rate=1/con$pr.mean[length(con$pr.mean)],log=T)
		return(pp)
	}

	# compute C
	C<-vcv.phylo(tree)
	invC<-solve(C)
	detC<-determinant(C,logarithm=TRUE)$modulus[1]

	# now set starting values for MCMC
	sig2<-con$sig2; a<-con$a; xbar<-con$xbar; intV<-con$intV
	L<-likelihood(C,invC,detC,x,sig2,a,xbar,intV)
	Pr<-log.prior(sig2,a,xbar,intV)

	# store
	X<-matrix(NA,ngen/con$sample+1,n+5,dimnames=list(NULL,c("gen","sig2","a",tree$tip.label,"var","logLik")))
	X[1,]<-c(0,sig2,a,xbar,intV,L)

	message("Starting MCMC...")

	# start MCMC
	for(i in 1:ngen){
		j<-(i-1)%%(n+3)
		if(j==0){
			# update sig2
			sig2.prime<-sig2+rnorm(n=1,sd=sqrt(con$prop[j+1]))
			if(sig2.prime<0) sig2.prime<--sig2.prime
			L.prime<-likelihood(C,invC,detC,x,sig2.prime,a,xbar,intV)
			Pr.prime<-log.prior(sig2.prime,a,xbar,intV)
			post.odds<-min(1,exp(Pr.prime+L.prime-Pr-L),na.rm=T)
			if(post.odds>runif(n=1)){
				if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2.prime,a,xbar,intV,L.prime)
				sig2<-sig2.prime
				L<-L.prime
				Pr<-Pr.prime
			} else if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,xbar,intV,L)
		} else if(j==1){
			# update a
			a.prime<-a+rnorm(n=1,sd=sqrt(con$prop[j+1]))
			L.prime<-likelihood(C,invC,detC,x,sig2,a.prime,xbar,intV)
			Pr.prime<-log.prior(sig2,a.prime,xbar,intV)
			post.odds<-min(1,exp(Pr.prime+L.prime-Pr-L),na.rm=T)
			if(post.odds>runif(n=1)){ 
				if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a.prime,xbar,intV,L.prime)
				a<-a.prime
				L<-L.prime
				Pr<-Pr.prime
			} else if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,xbar,intV,L)
		} else if(j>1&&j<=(n+1)) {
			k<-j-1 # update tip mean k
			xbar.prime<-xbar
			xbar.prime[k]<-xbar[k]+rnorm(n=1,sd=sqrt(con$prop[j+1]))
			L.prime<-likelihood(C,invC,detC,x,sig2,a,xbar.prime,intV)
			Pr.prime<-log.prior(sig2,a,xbar.prime,intV)
			post.odds<-min(1,exp(Pr.prime+L.prime-Pr-L),na.rm=T)
			if(post.odds>runif(n=1)){ 
				if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,xbar.prime,intV,L.prime)
				xbar<-xbar.prime
				L<-L.prime
				Pr<-Pr.prime
			} else if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,xbar,intV,L)
		} else if(j>(n+1)){
			# update var
			intV.prime<-intV+rnorm(n=1,sd=sqrt(con$prop[j+1]))
			if(intV.prime<0) intV.prime<--intV.prime
			L.prime<-likelihood(C,invC,detC,x,sig2,a,xbar,intV.prime)
			Pr.prime<-log.prior(sig2,a,xbar,intV.prime)
			post.odds<-min(1,exp(Pr.prime+L.prime-Pr-L),na.rm=T)
			if(post.odds>runif(n=1)){ 
				if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,xbar,intV.prime,L.prime)
				intV<-intV.prime
				L<-L.prime
				Pr<-Pr.prime
			} else if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,xbar,intV,L)	
		}
	}

	# done MCMC
	message("Done MCMC.")
	return(X)
}
