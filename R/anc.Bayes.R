## function does Bayes ancestral character estimation
## written by Liam J. Revell 2011, 2013, 2015, 2017, 2020

anc.Bayes<-function(tree,x,ngen=10000,control=list()){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	# give the function some defaults (in case none are provided)
	temp<-phyl.vcv(as.matrix(x),vcv(tree),1)
	sig2<-temp$R[1,1]
	a<-temp$alpha
	y<-rep(a,tree$Nnode-1)
	pr.mean<-c(1000,rep(0,tree$Nnode))
	pr.var<-c(pr.mean[1]^2,rep(1000,tree$Nnode))
	prop<-rep(0.01*max(temp$C)*sig2,tree$Nnode+1)

	# populate control list
	con=list(sig2=sig2,a=a,y=y,pr.mean=pr.mean,pr.var=pr.var,prop=prop,sample=100)
	con[(namc<-names(control))]<-control
	con<-con[!sapply(con,is.null)]
	# print control parameters to screen
	message("Control parameters (set by user or default):"); str(con)

	# function returns the log-likelihood
	likelihood<-function(C,invC,detC,x,sig2,a,y){
		z<-c(x,y)-rep(a,nrow(C))
		logLik<--z%*%invC%*%z/(2*sig2)-nrow(C)*log(2*pi)/2-nrow(C)*log(sig2)/2-detC/2
		return(logLik)
	}
	
	# function returns the prior
	log.prior<-function(pr.mean,pr.var,sig2,a,y)
		pp<-dexp(sig2,rate=1/pr.mean[1],log=TRUE)+sum(dnorm(c(a,y),mean=pr.mean[2:length(pr.mean)],sd=sqrt(pr.var[1+1:tree$Nnode]),log=TRUE))

	# compute C
	C<-vcvPhylo(tree)
	# check to make sure that C will be non-singular
	if(any(tree$edge.length<=(10*.Machine$double.eps)))
		stop("some branch lengths are 0 or nearly zero")
	invC<-solve(C)
	detC<-determinant(C,logarithm=TRUE)$modulus[1]

	# now set starting values for MCMC
	sig2<-con$sig2; a<-con$a; y<-con$y
	x<-x[tree$tip.label]
	if(is.null(names(y))) names(y)<-length(tree$tip)+2:tree$Nnode
	else y[as.character(length(tree$tip)+2:tree$Nnode)]
	L<-likelihood(C,invC,detC,x,sig2,a,y)
	Pr<-log.prior(con$pr.mean,con$pr.var,sig2,a,y)

	# store
	X<-matrix(NA,ngen/con$sample+1,tree$Nnode+3,dimnames=list(NULL,c("gen","sig2",length(tree$tip)+1:tree$Nnode,"logLik")))
	X[1,]<-c(0,sig2,a,y,L)

	message("Starting MCMC...")

	# start MCMC
	for(i in 1:ngen){
		j<-(i-1)%%(tree$Nnode+1)
		if(j==0){
			# update sig2
			sig2.prime<-sig2+rnorm(n=1,sd=sqrt(con$prop[j+1]))
			if(sig2.prime<0) sig2.prime<--sig2.prime
			L.prime<-likelihood(C,invC,detC,x,sig2.prime,a,y)
			Pr.prime<-log.prior(con$pr.mean,con$pr.var,sig2.prime,a,y)
			post.odds<-min(1,exp(Pr.prime+L.prime-Pr-L))
			if(post.odds>runif(n=1)){ 
				if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2.prime,a,y,L.prime)
				sig2<-sig2.prime
				L<-L.prime
				Pr<-Pr.prime
			} else if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,y,L)
		} else if(j==1){
			# update a
			a.prime<-a+rnorm(n=1,sd=sqrt(con$prop[j+1]))
			L.prime<-likelihood(C,invC,detC,x,sig2,a.prime,y)
			Pr.prime<-log.prior(con$pr.mean,con$pr.var,sig2,a.prime,y)
			post.odds<-min(1,exp(Pr.prime+L.prime-Pr-L))
			if(post.odds>runif(n=1)){ 
				if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a.prime,y,L.prime)
				a<-a.prime
				L<-L.prime
				Pr<-Pr.prime
			} else if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,y,L)
		} else {
			k<-j-1 # update node k
			y.prime<-y
			y.prime[k]<-y[k]+rnorm(n=1,sd=sqrt(con$prop[j+1]))
			L.prime<-likelihood(C,invC,detC,x,sig2,a,y.prime)
			Pr.prime<-log.prior(con$pr.mean,con$pr.var,sig2,a,y.prime)
			post.odds<-min(1,exp(Pr.prime+L.prime-Pr-L))
			if(post.odds>runif(n=1)){ 
				if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,y.prime,L.prime)
				y<-y.prime
				L<-L.prime
				Pr<-Pr.prime
			} else if(i%%con$sample==0) X[i/con$sample+1,]<-c(i,sig2,a,y,L)
		}
	}

	# done MCMC
	message("Done MCMC.")
	obj<-list(mcmc=as.data.frame(X),tree=tree)
	class(obj)<-"anc.Bayes"
	obj
}

## S3 methods
print.anc.Bayes<-function(x,digits=6,printlen=NULL,...){
	cat("\nObject of class \"anc.Bayes\" consisting of a posterior")
	cat("\n   sample from a Bayesian ancestral state analysis:\n")
	if(hasArg(burnin)) burnin<-list(...)$burnin
	else burnin<-0.2*max(x$mcmc$gen)
	ii<-which(((x$mcmc$gen-burnin)^2)==min((x$mcmc$gen-burnin)^2))
	Nnode<-x$tree$Nnode
	cat("\nMean ancestral states from posterior distribution:\n")
	ace<-colMeans(x$mcmc[ii:nrow(x$mcmc),as.character(1:Nnode+Ntip(x$tree))])
	if(is.null(printlen)||printlen>=Nnode) print(round(ace,digits))
	else printDotDot(ace,digits,printlen)
	cat(paste("\nBased on a burn-in of ",burnin," generations.\n\n",sep=""))
	invisible(ace)
}			

summary.anc.Bayes<-function(object,...) print(object,...)

plot.anc.Bayes<-function(x,...){
	args<-list(...)
	if(is.null(args$what)) what<-"logLik"
	else {
		what<-args$what
		args$what<-NULL
	}
	if(is.null(args$burnin)) burnin<-0.2*max(x$mcmc$gen)
	else {
		burnin<-args$burnin
		args$burnin<-NULL
	}
	if(what=="logLik"){
		args$x<-x$mcmc$gen
		args$y<-x$mcmc$logLik
		if(is.null(args$xlab)) args$xlab<-"generation"
		if(is.null(args$ylab)) args$ylab<-"log(L)"
		if(is.null(args$type)) args$type<-"l"
		if(is.null(args$col)) args$col<-make.transparent("blue",0.5)
		do.call(plot,args)
	} else {
		args$x<-x$mcmc$gen
		args$y<-x$mcmc[,as.character(what)]
		if(is.null(args$xlab)) args$xlab<-"generation"
		if(is.null(args$ylab)) args$ylab<-paste("state for node",what)
		if(is.null(args$type)) args$type<-"l"
		if(is.null(args$col)) args$col<-make.transparent("blue",0.5)
		if(is.null(args$ylim)) args$ylim<-range(x$mcmc[,2:(ncol(x$mcmc)-1)])
		do.call(plot,args)
	}
}

density.anc.Bayes<-function(x,...){
	if(hasArg(what)) what<-list(...)$what
	else what<-Ntip(x$tree)+1
	if(hasArg(burnin)) burnin<-list(...)$burnin
	else burnin<-0.2*max(x$mcmc$gen)
	ii<-which(abs(x$mcmc$gen-burnin)==min(abs(x$mcmc$gen-burnin)))
	args<-list()
	if(hasArg(bw)) args$bw<-list(...)$bw
	else args$bw<-"nrd0"
	args$x<-x$mcmc[ii:nrow(x$mcmc),as.character(what)]
	d<-do.call(density,args)
	d$call<-match.call()
	d$data.name<-if(what=="logLik") what else paste("node",what)
	d
}