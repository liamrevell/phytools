# function for phylogenetic paired t-test (Lindenfors et al. 2010)
# written by Liam Revell 2011, 2013, 2015

phyl.pairedttest<-function(tree,x1,x2=NULL,se1=NULL,se2=NULL,lambda=1.0,h0=0.0,fixed=FALSE){
	# check tree
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	# convert x1 to a matrix if necessary
	if(is.data.frame(x1)) x1<-as.matrix(x1)
	# if x2 is NULL
	if(is.null(x2)){
		# check to see that x1 has two columns
		if(dim(x1)[2]!=2) stop("user must provide two data vectors or matrix with two variables")
		else {
			x2<-x1[,2]
			x1<-x1[,1]
		}
	}
	if(is.data.frame(x2)) x2<-as.matrix(x2) # to be safe
	if(is.matrix(x1)) x1<-x1[,1]
	if(is.matrix(x2)) x2<-x2[,1]
	if(is.null(se1)){ 
		v1<-rep(0,length(tree$tip))
		names(v1)<-tree$tip.label
	} else v1<-se1^2
	if(is.null(se2)){
		v2<-rep(0,length(tree$tip))
		names(v2)<-tree$tip.label
	} else v2<-se2^2
	# compute C and sort tree$tip.label
	C<-vcv.phylo(tree)
	x1<-x1[tree$tip.label]
	x2<-x2[tree$tip.label]
	v1<-v1[tree$tip.label]
	v2<-v2[tree$tip.label]
	V.diff<-diag(v1+v2); dimnames(V.diff)<-list(tree$tip.label,tree$tip.label)
	# compute difference
	d<-x1-x2
	# lambda transformation
	lambda.transform<-function(C,lambda) lambda*(C-diag(diag(C)))+diag(diag(C))
	# likelihood function
	likelihood<-function(theta,d,C,V.diff){
		theta[1]->sig2
		theta[2]->lambda
		theta[3]->dbar
		V<-sig2*lambda.transform(C,lambda)+V.diff
		logL<-as.numeric(-t(d-dbar)%*%solve(V,d-dbar)/2-determinant(V)$modulus[1]/2-length(d)*log(2*pi)/2)
		return(logL)
	}
	# maximize the likelihood
	if(!fixed) res=optim(c(mean(pic(d,multi2di(tree))^2),lambda,h0),likelihood,d=d,C=C,V.diff=V.diff,method="L-BFGS-B",lower=c(1e-8,0,-Inf),upper=c(Inf,1,Inf),hessian=TRUE,control=list(fnscale=-1))
	else res=optim(c(mean(pic(d,multi2di(tree))^2),lambda,h0),likelihood,d=d,C=C,V.diff=V.diff,method="L-BFGS-B",lower=c(1e-8,lambda-1e-8,-Inf),upper=c(Inf,lambda,Inf),hessian=TRUE,control=list(fnscale=-1))
	# test
	se.dbar<-sqrt(-1/res$hessian[3,3])
	t<-(res$par[3]-h0)/se.dbar
	P<-2*pt(abs(t),df=length(tree$tip)-3,lower.tail=F)
	return(list(dbar=res$par[3],se=se.dbar,sig2=res$par[1],lambda=round(res$par[2],7),logL=res$value,t.dbar=t,P.dbar=P))
}
