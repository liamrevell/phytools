## Function fits the threshold model for two characters using Bayesian MCMC.
## All characters should be provided either in the form of a numeric matrix, X, or as
## a data frame in which the discrete character is coded as a factor
## Row names of X should match the species names of the tree.
## types=c("discrete","discrete"), c("discrete","continuous"), c("cont","disc") etc. should 
## be used to indicate the data type of each column in X
## written by Liam J. Revell 2012, 2014, 2017, 2018

threshBayes<-function(tree,X,types=NULL,ngen=10000,control=list(),...){

	if(hasArg(burnin)) burnin<-list(...)$burnin
	else burnin<-round(0.2*ngen)

	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")

	# likelihood function for the liabilities
	lik<-function(a,V,invV,detV,D,Y){
		y<-as.vector(Y)
		logL<--t(y-D%*%a)%*%invV%*%(y-D%*%a)/2-n*m*log(2*pi)/2-detV/2
		return(logL)
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
		if(length(d)>0) for(i in 1:length(d)) levels[[d[i]]]<-as.character(0:1)
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
	} else {
		Y<-X
		types<-crop(types)
		disc<-which(types=="disc")
		for(i in 1:length(disc)) Y[Y[,disc[i]]==0,disc[i]]<--1
	}
	npar<-npar-n*(m-length(disc))
	if(is.null(colnames(X))) colnames(Y)<-paste("V",1:m,sep="")
	else colnames(Y)<-colnames(X)
	r<-0 # for the correlation
	sig2<-colMeans(apply(Y,2,pic,phy=multi2di(tree))^2); sig2[disc]<-1
	a<-colMeans(Y); a[disc]<-0
	R<-matrix(c(sig2[1],r*sqrt(sig2[1]*sig2[2]),r*sqrt(sig2[1]*sig2[2]),sig2[2]),2,2)
	V<-kronecker(R,C); invV<-solve(V); detV<-determinant(V)$modulus[1]

	# populate control list
	con=list(sample=100,
			propvar=c(0.3*sig2,0.3*sqrt(sig2),0.2),
			propliab=0.3,
			pr.mean=c(rep(1000,m),rep(0,m),0),
			pr.liab=0,
			pr.var=c(rep(1000,m)^2,rep(1000,m),4/12),
			pr.vliab=1000,
			quiet=FALSE)
	con[(namc<-names(control))]<-control
	con<-con[!sapply(con,is.null)]
	
	# for updating
	P<-Pp<-matrix(NA,n,m)
	P[Y<0]<-0
	P[Y>0]<-1
	Vp<-V
	invVp<-invV
	detVp<-detV
	
	# for storing the posterior sample
	Z<-matrix(NA,ngen/con$sample+1,7,dimnames=list(NULL,c("gen","sig1","sig2",
		"a1","a2","r","logL")))
	Z[1,]<-c(0,sig2,a,r,lik(a,V,invV,detV,D,Y))
	L<-matrix(NA,ngen/con$sample+1,m*Ntip(tree)+1,dimnames=list(NULL,c("gen",
		as.vector(apply(matrix(colnames(Y)),1,paste,".",tree$tip.label,sep="")))))
	L[1,]<-c(0,as.vector(Y))

	# start MCMC
	cat("Starting MCMC....\n")
	flush.console()
	for(i in 1:ngen){
		lik1<-lik(a,V,invV,detV,D,Y)+log(all(P[,disc]==X[,disc]))
		d<-i%%npar
		if(ngen>=1000) if(i%%1000==0) if(!con$quiet){ 
			cat(paste("gen ",i,"\n",sep=""))
			flush.console()
		}
		Yp<-Y
		sig2p<-sig2
		ap<-a
		rp<-r
		if(d<=length(Y[,disc])&&d>0){
			# update liabilities
			ind<-c(d%%n,disc[ceiling(d/n)])
			if(ind[1]==0) ind[1]<-n
			Yp[ind[1],ind[2]]<-Y[ind[1],ind[2]]+rnorm(n=1,sd=sqrt(con$propliab))
		} else {
			if((d-length(Y[,disc]))==1||(d-length(Y[,disc]))==2){
				# update sig2
				if(!((d-length(Y[,disc]))%in%disc)){
					j<-d-length(Y[,disc])
					sig2p[j]<-sig2[j]+rnorm(n=1,sd=sqrt(con$propvar[j]))
					if(sig2p[j]<0) sig2p[j]=-sig2p[j]
					R<-matrix(c(sig2p[1],rp*sqrt(sig2p[1]*sig2p[2]),
						rp*sqrt(sig2p[1]*sig2p[2]),sig2p[2]),2,2)
					Vp<-kronecker(R,C)
					invVp<-solve(Vp)
					detVp<-determinant(Vp)$modulus[1]
				}
			} else if((d-length(Y[,disc]))==3||(d-length(Y[,disc]))==4){
				# update a
				j<-d-length(Y[,disc])
				ap[j-2]<-a[j-2]+rnorm(n=1,sd=sqrt(con$propvar[j]))
			} else if(d==0) {
				# update r
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
		if(p.odds>runif(n=1)){
			Y<-Yp
			sig2<-sig2p
			a<-ap
			r<-rp
			V<-Vp
			invV<-invVp
			detV<-detVp
			logL<-lik2
		} else logL<-lik1
		if(i%%con$sample==0){ 
			Z[i/con$sample+1,]<-c(i,sig2,a,r,logL)
			L[i/con$sample+1,]<-c(i,Y[,1],Y[,2])
		}
	}
	cat("Done MCMC.\n")
	obj<-list(par=as.data.frame(Z),liab=as.data.frame(L),
		burnin=burnin,levels=levels,types=types)
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

plot.threshBayes<-function(x,...){
	if(hasArg(burnin)) burnin<-list(...)$burnin
	else burnin<-x$burnin
	if(hasArg(bw)) bw<-list(...)$bw
	else bw<-0.05
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-c(-1,1)
	ii<-which(((x$par$gen-burnin)^2)==min((x$par$gen-burnin)^2))[1]+1
	d<-density(x$par$r[ii:nrow(x$par)],bw=bw)
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-c(0,1.2*max(d$y))
	plot(d,xlim=xlim,ylim=ylim,col="blue",xlab="Posterior sample of r",
		ylab="Density",main="")
	polygon(x=c(min(d$x),d$x,max(d$x)),y=c(0,d$y,0),
		col=make.transparent("blue",0.2))
	r<-mean(x$par$r[ii:nrow(x$par)])
	lines(rep(r,2),c(0,par()$usr[4]),col="blue",lty="dashed",
		lwd=2)
	text(r,0.95*par()$usr[4],"mean post-burnin\nvalue of r",cex=0.7,
		pos=if(r>0) 2 else 4)
}