## Function fits the threshold model for two characters using Bayesian MCMC.
## All characters should be provided either in the form of a numeric matrix, X, or as
## a data frame in which the discrete character is coded as a factor
## Row names of X should match the species names of the tree.
## types=c("discrete","discrete"), c("discrete","continuous"), c("cont","disc") etc. should 
## be used to indicate the data type of each column in X
## written by Liam J. Revell 2012, 2014, 2017, 2018, 2020, 2021

threshBayes<-function(tree,X,types=NULL,ngen=10000,control=list(),...){

	if(hasArg(burnin)) burnin<-list(...)$burnin
	else burnin<-round(0.2*ngen)
	
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	
	if(hasArg(auto.tune)) auto.tune<-list(...)$auto.tune
	else auto.tune<-TRUE
	if(is.logical(auto.tune)) if(auto.tune==TRUE) auto.tune<-0.234
	else if(is.numeric(auto.tune)) if(auto.tune>=1.0||auto.tune<=0.0){
		cat("value for auto.tune outside allowable range. Resetting....\n")
		auto.tune<-0.234
	}

	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")

	# likelihood function for the liabilities
	lik<-function(a,V,invV,detV,D,Y){
		y<-as.vector(Y)
		logL<--t(y-D%*%a)%*%invV%*%(y-D%*%a)/2-n*m*log(2*pi)/2-detV/2
		return(logL[1,1])
	}

	# function for the log-prior (presently returns 0, i.e. a flat prior)
	logPrior<-function(sig2,a,r,Yp){
		pp<-dexp(sig2[1],rate=1/con$pr.mean[1],log=T)+
			dexp(sig2[2],rate=1/con$pr.mean[2],log=T)+
			dnorm(a[1],mean=con$pr.mean[3],sd=sqrt(con$pr.var[3]),log=T)+
			dnorm(a[2],mean=con$pr.mean[4],sd=sqrt(con$pr.var[4]),log=T)+
			dunif(r,min=con$pr.mean[5]-0.5*sqrt(12*con$pr.var[5]),
				max=con$pr.mean[5]+0.5*sqrt(12*con$pr.var[5]),log=T)
		return(pp)		
	}

	# function to crop the first 4 characters from a vector of strings
	crop<-function(x){
		for(i in 1:length(x)) x[i]<-paste(strsplit(x[i],"")[[1]][1:4],collapse="")
		return(x)
	}
	
	## bookkeeping
	C<-vcv.phylo(tree)
	n<-nrow(X)
	m<-ncol(X)
	if(is.null(types)){
		types<-rep("cont",m)
		if(is.data.frame(X)){
			ii<-which(sapply(X,is.factor))
			if(length(ii)>0) types[ii]<-"disc"
		} else if(is.matrix(X)){
			ii<-which(apply(X,2,function(x) all(x%in%c(0,1))))
			if(length(ii)>0) types[ii]<-"disc"
		}
	}		
	X<-X[tree$tip.label,]
	levels<-as.list(rep(NA,m))
	if(is.data.frame(X)){
		ii<-sapply(X,is.character)
		if(any(ii)) X[,ii]<-sapply(X[,ii],as.factor)
		ii<-sapply(X,is.factor)
		if(any(ii)){
			d<-which(ii)
			for(i in 1:length(d)) levels[[d[i]]]<-levels(X[,d[i]])
			X[,ii]<-sapply(X[,ii],as.numeric)-1
		}
		X<-as.matrix(X)
	} else {
		d<-which(crop(types)=="disc")
		if(length(d)>0) for(i in 1:length(d)) 
			levels[[d[i]]]<-as.character(0:1)
	}

	if(m!=2) stop("number of traits must equal 2")
	npar<-n*m+5 # or npar<-n*m+1
	
	# design matrix
	D<-matrix(0,n*m,m)
	for(i in 1:(n*m)) for(j in 1:m) if((j-1)*n<i&&i<=j*n) D[i,j]<-1
	
	# set starting values
	if(is.null(types)){
		Y<-matrix(-1,Ntip(tree),2,dimnames=list(rownames(X),NULL))
		Y[X==1]<-1 # for liability
		disc<-1:m
		cont<-vector()
	} else {
		Y<-X
		types<-crop(types)
		disc<-which(types=="disc")
		cont<-which(types=="cont")
		for(i in 1:length(disc)) Y[Y[,disc[i]]==0,disc[i]]<--1
	}
	Y[,disc]<-Y[,disc]*abs(rnorm(n=length(as.vector(Y[,disc])),
		sd=sqrt(max(nodeHeights(tree)))))

	npar<-npar-n*(m-length(disc))-length(disc)
	if(is.null(colnames(X))) colnames(Y)<-paste("V",1:m,sep="")
	else colnames(Y)<-colnames(X)
	r<-0 # for the correlation
	sig2<-colMeans(apply(Y,2,pic,phy=multi2di(tree,random=FALSE))^2)
	sig2[disc]<-1
	a<-colMeans(Y)
	a[disc]<-0
	R<-matrix(c(sig2[1],r*sqrt(sig2[1]*sig2[2]),r*sqrt(sig2[1]*sig2[2]),
		sig2[2]),2,2)
	V<-kronecker(R,C)
	invV<-solve(V)
	detV<-determinant(V)$modulus[1]

	# populate control list
	con=list(sample=100,
			propvar=c(0.3*sig2,0.3*sqrt(sig2),0.2),
			propliab=0.3,
			pr.mean=c(rep(1000,m),rep(0,m),0),
			pr.liab=0,
			pr.var=c(rep(1000,m)^2,rep(1000,m),4/12),
			pr.vliab=1000,
			quiet=FALSE,
			print.interval=1000)
	con[(namc<-names(control))]<-control
	con<-con[!sapply(con,is.null)]

	## permit different proposal variances for each tip/trait
	if(length(con$propliab)!=(n*(m-length(cont)))) con$propliab<-rep(con$propliab,n*(m-length(cont)))
		
	# for updating
	P<-Pp<-matrix(NA,n,m)
	P[Y<0]<-0
	P[Y>0]<-1
	Vp<-V
	invVp<-invV
	detVp<-detV
	
	# for storing the posterior sample
	Z<-matrix(NA,ngen/con$sample+1,8,dimnames=list(NULL,c("gen","sig1","sig2",
		"a1","a2","r","logL","accept_rate")))
	Z[1,]<-c(0,sig2,a,r,lik1<-lik(a,V,invV,detV,D,Y),0)
	L<-matrix(NA,ngen/con$sample+1,m*Ntip(tree)+1,dimnames=list(NULL,c("gen",
		as.vector(apply(matrix(colnames(Y)),1,paste,".",tree$tip.label,sep="")))))
	L[1,]<-c(0,as.vector(Y))

	logL.trace<-r.trace<-rep(0,ngen)
	running.accept<-0
	
	## maybe help with memory
	Yp<-Y
	
	# start MCMC
	cat("Starting MCMC....\n")
	flush.console()
	for(i in 1:ngen){
		if(i%%10000==1){
			## reset acceptance rates
			accept.rate<-accept<-hit<-rep(0,npar)
		}
		d<-i%%npar
		if(ngen>=con$print.interval) if(i%%con$print.interval==0) if(!con$quiet){ 
			cat(paste("generation: ",i,"; mean acceptance rate: ",round(mean(accept.rate),2),"\n",sep=""))
			flush.console()
		}
		Yp[]<-Y
		sig2p<-sig2
		ap<-a
		rp<-r
		if(d<=length(Y[,disc])&&d>0){
			# update liabilities
			Yp[,disc][d]<-Y[,disc][d]+rnorm(n=1,sd=sqrt(con$propliab[d]))
		} else {
			key<-c("sig1","sig2","a1","a2","r")[c(cont,3,4,5)]
			key<-key[if(d>0) d-length(Y[,disc]) else length(key)]
			if(key=="sig1"){
				sig2p[1]<-sig2[1]+rnorm(n=1,sd=sqrt(con$propvar[1]))
				if(sig2p[1]<0) sig2p[1]=-sig2p[1]
				R<-matrix(c(sig2p[1],rp*sqrt(sig2p[1]*sig2p[2]),
					rp*sqrt(sig2p[1]*sig2p[2]),sig2p[2]),2,2)
				Vp<-kronecker(R,C)
				invVp<-solve(Vp)
				detVp<-determinant(Vp)$modulus[1]
			} else if(key=="sig2"){
				sig2p[2]<-sig2[2]+rnorm(n=1,sd=sqrt(con$propvar[2]))
				if(sig2p[2]<0) sig2p[2]=-sig2p[2]
				R<-matrix(c(sig2p[1],rp*sqrt(sig2p[1]*sig2p[2]),
					rp*sqrt(sig2p[1]*sig2p[2]),sig2p[2]),2,2)
				Vp<-kronecker(R,C)
				invVp<-solve(Vp)
				detVp<-determinant(Vp)$modulus[1]
			} else if(key=="a1") {
				ap[1]<-a[1]+rnorm(n=1,sd=sqrt(con$propvar[3]))
			} else if(key=="a2") {
				ap[2]<-a[2]+rnorm(n=1,sd=sqrt(con$propvar[4]))
			} else if(key=="r"){
				rp<-r+rnorm(n=1,sd=sqrt(con$propvar[5]))
				while(rp>1||rp< -1){
					if(rp>1) rp<-2-rp
					if(rp< -1) rp<--2-rp
				}
				R<-matrix(c(sig2p[1],rp*sqrt(sig2p[1]*sig2p[2]),
					rp*sqrt(sig2p[1]*sig2p[2]),sig2p[2]),2,2)
				Vp<-kronecker(R,C)
				invVp<-solve(Vp)
				detVp<-determinant(Vp)$modulus[1]
			}
		}
		Pp[Yp<0]<-0
		Pp[Yp>0]<-1
		lik2<-lik(ap,Vp,invVp,detVp,D,Yp)+log(all(Pp[,disc]==X[,disc]))
		p.odds<-min(c(1,exp(lik2+logPrior(sig2p,ap,rp,Yp)-lik1-
			logPrior(sig2,a,r,Yp))))
		hit[if(d>0) d else npar]<-hit[if(d>0) d else npar]+1
		if(p.odds>runif(n=1)){
			accept[if(d>0) d else npar]<-accept[if(d>0) d else npar]+1
			running.accept<-running.accept+1/con$sample
			Y[]<-Yp
			sig2<-sig2p
			a<-ap
			r<-rp
			V<-Vp
			invV<-invVp
			detV<-detVp
			logL<-lik2
			lik1<-lik2
		} else logL<-lik1
		accept.rate[if(d>0) d else npar]<-accept[if(d>0) d else npar]/hit[if(d>0) d else npar]
		if(is.numeric(auto.tune)){
			if(accept.rate[if(d>0) d else npar]>auto.tune){
				if(d>=1&&d<=(length(Y[,disc]))) con$propliab[d]<-1.1*con$propliab[d]
				else {
					ii<-which(c("sig1","sig2","a1","a2","r")==key)
					con$propvar[ii]<-1.1*con$propvar[ii]
					if(key=="r"&&con$propvar[ii]>1) con$propvar[ii]<-1 
				}
			} else if(accept.rate[if(d>0) d else npar]<auto.tune){
				if(d>=1&&d<=(length(Y[,disc]))) con$propliab[d]<-con$propliab[d]/1.1
				else {
					ii<-which(c("sig1","sig2","a1","a2","r")==key)
					con$propvar[ii]<-con$propvar[ii]/1.1
				}
			}
		}
		logL.trace[i]<-logL
		r.trace[i]<-r
		if(plot&&(i%%100==1)){
			dev.hold()
			par(mfrow=c(3,1))
			par(mar=c(5.1,4.1,2.1,1.1))
			plot(1:i,logL.trace[1:i],type="l",bty="l",col=make.transparent("grey",0.5),
				xlab="generation",ylab="log(L)",xlim=c(0,ngen))
			mtext("a) log-likelihood trace",side=3,line=0,cex=1,at=0,outer=FALSE,
				adj=0)
			par(mar=c(5.1,4.1,2.1,1.1))
			h<-barplot(accept.rate,ylim=c(0,1),col=make.transparent("blue",0.25),
				border=NA)
			if(is.numeric(auto.tune)) lines(range(h),rep(auto.tune,2),lty="dotted")
			mtext("b) acceptance rates (resets every 10,000 generations)",side=3,line=0,
				cex=1,at=0,outer=FALSE,adj=0)
			plot(1:i,r.trace[1:i],type="l",bty="l",col=make.transparent("blue",0.5),
				xlab="generation",ylab="r",xlim=c(0,ngen))
			mtext("c) trace of the correlation coefficient, r",side=3,line=0,cex=1,at=0,
				outer=FALSE,adj=0)
			lines(c(par()$usr[1],ngen),rep(mean(r.trace[1:i]),2),lty="dotted")
			dev.flush()
		}
		if(i%%con$sample==0){ 
			Z[i/con$sample+1,]<-c(i,sig2,a,r,logL,running.accept)
			L[i/con$sample+1,]<-c(i,Y[,1],Y[,2])
			running.accept<-0 ## reset running acceptance rate
		}
	}
	cat("Done MCMC.\n")
	obj<-list(par=as.data.frame(Z),liab=as.data.frame(L),
		burnin=burnin,levels=levels,types=types)
	attr(obj,"auto.tune")<-auto.tune
	class(obj)<-"threshBayes"
	obj
}

## S3 methods for the object class

print.threshBayes<-function(x,...){
	cat("\nObject of class \"threshBayes\" consisting of a matrix (L) of\n")
	cat("sampled liabilities for the tips of the tree & a second matrix\n")
	cat("(par) with the sample model parameters & correlation.\n")
	cat(paste("\nMean correlation (r) from the posterior sample is: ",
		round(mean(x$par[,"r"]),5),".\n",sep=""))
	if(any(x$types=="disc")){
		cat("\nOrdination of discrete traits:\n\n")
		for(i in 1:length(x$types)){		
			if(x$types[i]=="disc") cat(paste("\tTrait ",i,": ",x$levels[[i]][1]," <-> ",
				x$levels[[i]][2],"\n",sep=""))
		}
	}
	cat("\n")
}

plot.density.threshBayes<-function(x,...){
	d<-x$density
	r<-x$mean.r
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-c(-1,1)
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-c(0,1.2*max(d$y))
	if(hasArg(bty)) bty<-list(...)$bty
	else bty<-"n"
	if(hasArg(cex.lab)) cex.lab<-list(...)$cex.lab
	else cex.lab<-1
	if(hasArg(cex.axis)) cex.axis<-list(...)$cex.axis
	else cex.axis<-1
	plot(d,xlim=xlim,ylim=ylim,col="blue",xlab="Posterior sample of r",
		ylab="Density",main="",bty=bty,cex.lab=cex.lab,cex.axis=cex.axis)
	polygon(x=c(min(d$x),d$x,max(d$x)),y=c(0,d$y,0),
		col=make.transparent("blue",0.2))
	lines(rep(r,2),c(0,max(d$y)),col="blue",lty="dashed",
		lwd=2)
	text(r,max(d$y),"mean post-burnin\nvalue of r",cex=0.7,
		pos=if(r>0) 2 else 4,font=3)
}

density.threshBayes<-function(x,...){
	if(hasArg(burnin)) burnin<-list(...)$burnin
	else burnin<-x$burnin
	if(hasArg(bw)) bw<-list(...)$bw
	else bw<-0.05
	ii<-which(((x$par$gen-burnin)^2)==min((x$par$gen-burnin)^2))[1]+1
	d<-density(x$par$r[ii:nrow(x$par)],bw=bw)
	r<-mean(x$par$r[ii:nrow(x$par)])
	d$data.name<-deparse(substitute(x))
	d$call<-paste("density.threshBayes(x =",d$data.name," bw = bw)")
	d$bw<-bw
	object<-list(density=d,mean.r=r)
	class(object)<-"density.threshBayes"
	object
}

print.density.threshBayes<-function(x,...){
	print(x$density)
	cat("\n")
	cat("To plot enter plot(\'object_name\') at the command line interface.\n\n")
}

plot.threshBayes<-function(x,...){
	if(hasArg(bw)) bw<-list(...)$bw
	else bw<-floor(length(x$par$gen)/100)
	if(hasArg(bty)) bty<-list(...)$bty
	else bty<-"n"
	if(hasArg(las)) las<-list(...)$las
	else las<-1
	if(hasArg(cex.main)) cex.main<-list(...)$cex.main
	else cex.main<-1
	par(mfrow=c(3,1))
	par(mar=c(5.1,4.1,2.1,1.1))
	plot(x$par$gen,x$par$logL,type="l",bty=bty,col=make.transparent("grey",0.5),
		xlab="generation",ylab="log(L)",las=las)
	mtext("a) log-likelihood trace",side=3,line=0.5,cex=cex.main,at=0,
		outer=FALSE,adj=0)
	par(mar=c(5.1,4.1,2.1,1.1))
	accept<-vector()
	for(i in 1:length(x$par$gen))
		accept[i]<-mean(x$par$accept_rate[max(c(1,i-bw)):i])
	plot(x$par$gen,accept,type="l",bty=bty,col=make.transparent("red",0.5),
		xlab="generation",ylab="mean acceptance rate",las=las)
	if(is.numeric(attr(x,"auto.tune"))) lines(c(par()$usr[1],max(x$par$gen)),
		rep(attr(x,"auto.tune"),2),lty="dotted")
	mtext(paste("b) mean acceptance rate (sliding window: bw=",bw,")",sep=""),
		side=3,line=0.5,cex=cex.main,at=0,outer=FALSE,adj=0)
	plot(x$par$gen,x$par$r,type="l",bty=bty,col=make.transparent("blue",0.5),
		xlab="generation",ylab="r",las=las)
	mtext("c) trace of the correlation coefficient, r",side=3,line=0.5,
		cex=cex.main,at=0,outer=FALSE,adj=0)
}
