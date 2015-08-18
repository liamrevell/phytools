# function fits the threshold model for two characters using Bayesian MCMC
# all characters should be provided in the numeric matrix X: discrete characters as 0,1
# row names of X should match the species names of the tree
# types=c("discrete","discrete"), c("discrete","continuous"), c("cont","disc") etc. should 
# be used to indicate the data type of each column in X
# written by Liam J. Revell 2012, 2014

threshBayes<-function(tree,X,types=NULL,ngen=1000,control=list()){

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
			dunif(r,min=con$pr.mean[5]-0.5*sqrt(12*con$pr.var[5]),max=con$pr.mean[5]+0.5*sqrt(12*con$pr.var[5]),log=T)
		return(pp)		
	}

	# function to crop the first 4 characters from a vector of strings
	crop<-function(x){
		for(i in 1:length(x)) x[i]<-paste(strsplit(x[i],"")[[1]][1:4],collapse="")
		return(x)
	}

	# bookkeeping
	X<-X[tree$tip.label,]
	C<-vcv.phylo(tree)
	n<-nrow(X);	m<-ncol(X)
	if(m!=2) stop("number of traits must equal 2")
	npar<-n*m+5 # or npar<-n*m+1
	
	# design matrix
	D<-matrix(0,n*m,m)
	for(i in 1:(n*m)) for(j in 1:m) if((j-1)*n<i&&i<=j*n) D[i,j]<-1
	
	# set starting values
	if(is.null(types)){
		Y<-matrix(-1,length(tree$tip),2,dimnames=list(rownames(X),NULL))
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
	P[Y<0]<-0; P[Y>0]<-1
	Vp<-V; invVp<-invV; detVp<-detV
	
	# for storing the posterior sample
	Z<-matrix(NA,ngen/con$sample+1,7,dimnames=list(NULL,c("gen","sig1","sig2","a1","a2","r","logL")))
	Z[1,]<-c(0,sig2,a,r,lik(a,V,invV,detV,D,Y))
	L<-matrix(NA,ngen/con$sample+1,m*length(tree$tip)+1,dimnames=list(NULL,c("gen",as.vector(apply(matrix(colnames(Y)),1,paste,".",tree$tip.label,sep="")))))
	L[1,]<-c(0,as.vector(Y))

	# start MCMC
	for(i in 1:ngen){
		lik1<-lik(a,V,invV,detV,D,Y)+log(all(P[,disc]==X[,disc]))
		d<-i%%npar
		if(ngen>=1000) if(i%%1000==0) if(!con$quiet) cat(paste("gen ",i,"\n",sep=""))
		Yp<-Y; sig2p<-sig2; ap<-a; rp<-r
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
					R<-matrix(c(sig2p[1],rp*sqrt(sig2p[1]*sig2p[2]),rp*sqrt(sig2p[1]*sig2p[2]),sig2p[2]),2,2)
					Vp<-kronecker(R,C); invVp<-solve(Vp); detVp<-determinant(Vp)$modulus[1]
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
				R<-matrix(c(sig2p[1],rp*sqrt(sig2p[1]*sig2p[2]),rp*sqrt(sig2p[1]*sig2p[2]),sig2p[2]),2,2)
				Vp<-kronecker(R,C); invVp<-solve(Vp); detVp<-determinant(Vp)$modulus[1]
			}
		}
		Pp[Yp<0]<-0; Pp[Yp>0]<-1
		lik2<-lik(ap,Vp,invVp,detVp,D,Yp)+log(all(Pp[,disc]==X[,disc]))
		p.odds<-min(c(1,exp(lik2+logPrior(sig2p,ap,rp,Yp)-lik1-logPrior(sig2,a,r,Yp))))
		if(p.odds>runif(n=1)){
			Y<-Yp; sig2<-sig2p; a<-ap; r<-rp
			V<-Vp; invV<-invVp; detV<-detVp
			logL<-lik2
		} else logL<-lik1
		if(i%%con$sample==0){ 
			Z[i/con$sample+1,]<-c(i,sig2,a,r,logL)
			L[i/con$sample+1,]<-c(i,Y[,1],Y[,2])
		}
	}
	return(list(par=Z,liab=L))
}
