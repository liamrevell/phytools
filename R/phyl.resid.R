# function to fit phylogenetic regression and compute residuals
# multiple morphological traits in Y, size in x
# written by Liam Revell 2011, 2015 ref. Revell (2009; Evolution)

phyl.resid<-function(tree,x,Y,method="BM"){
	# check tree
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	# check and sort data
	# X
	X<-cbind(1,x)
	if(is.null(rownames(X))){
		print("x has no names. function will assume that the order of x matches tree$tip.label")
		rownames(X)<-tree$tip.label
	} else X<-X[tree$tip.label,] # sort
	# Y
	Y<-as.matrix(Y)
	if(is.null(rownames(Y))){
		print("y has no names. function will assume that the order of y matches tree$tip.label")
		rownames(Y)<-tree$tip.label
	} else Y<-as.matrix(Y[tree$tip.label,]) # sort

	# analyze
	if(method=="BM"){ 
		C<-vcv.phylo(tree)
		beta<-solve(t(X)%*%solve(C)%*%X)%*%t(X)%*%solve(C)%*%Y
		resid<-Y-X%*%beta
		return(list(beta=beta,resid=resid))
	} else if(method=="lambda"){
		C<-vcv.phylo(tree)
		maxLambda<-max(C)/max(C[upper.tri(C)])
		lambda<-vector()
		logL<-vector()
		beta<-matrix(NA,ncol(X),ncol(Y),dimnames=list(colnames(X),colnames(Y)))
		for(i in 1:ncol(Y)){
			res<-optimize(f=likelihood.lambda,interval=c(0,maxLambda),y=Y[,i],X=X,C=C,maximum=TRUE)
			lambda[i]<-res$maximum
			logL[i]<-as.numeric(res$objective)
			C.l<-lambda.transform(lambda[i],C)
			beta[,i]<-solve(t(X)%*%solve(C.l)%*%X)%*%t(X)%*%solve(C.l)%*%Y[,i]
		}
		resid<-Y-X%*%beta
		return(list(beta=beta,lambda=lambda,logL=logL,resid=resid))
	}
}

# likelihood function for the regression model with lambda
likelihood.lambda<-function(lambda,y,X,C){
	n<-nrow(C)
	C.lambda<-lambda.transform(lambda,C)
	beta<-solve(t(X)%*%solve(C.lambda)%*%X)%*%(t(X)%*%solve(C.lambda)%*%y)
	sig2e<-as.double((1/n)*(t(y-X%*%beta)%*%solve(C.lambda)%*%(y-X%*%beta)))
	logL<--(1/2)*t(y-X%*%beta)%*%solve(sig2e*C.lambda)%*%(y-X%*%beta)-(1/2)*determinant(sig2e*C.lambda,logarithm=TRUE)$modulus-(n/2)*log(2*pi)
	return(logL)
}



