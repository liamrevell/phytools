# function does phylogenetic canonical correlation analysis (Revell & Harrison 2008)
# written by Liam Revell 2011, 2012, 2013, 2015

phyl.cca<-function(tree,X,Y,lambda=1.0,fixed=TRUE){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	# misc
	n<-length(tree$tip)
	mX<-ncol(X)
	mY<-ncol(Y)
	# compute C
	C<-vcv.phylo(tree)
	# if X & Y are data frames, convert to matrix
	if(is.data.frame(X)) X<-as.matrix(X)
	if(is.data.frame(Y)) Y<-as.matrix(Y)
	# check to see if X & Y have rownames
	if(is.null(rownames(X))){
		message("X is missing rownames; assuming same order as tree$tip.label")
		if(nrow(X)!=length(tree$tip)) warning("X does not have the correct number of rows")
		else rownames(X)<-tree$tip.label
	}
	if(is.null(rownames(Y))){
		message("Y is missing rownames; assuming same order as tree$tip.label")
		if(nrow(Y)!=length(tree$tip)) warning("Y does not have the correct number of rows")
		else rownames(Y)<-tree$tip.label
	}
	# reorder Y & C by tree$tip.label
	C<-C[tree$tip.label,tree$tip.label] # I think this is superfluous
	X<-as.matrix(X[tree$tip.label,])
	Y<-as.matrix(Y[tree$tip.label,])
	# set or optimize lambda
	if(fixed) logL<-likMlambda(lambda,cbind(X,Y),C)
	else {
		temp<-optimize(f=likMlambda,interval=c(0,maxLambda(tree)),X=cbind(X,Y),C=C,maximum=TRUE)
		logL<-temp$objective
		lambda<-temp$maximum
	}
	C<-lambda.transform(lambda,C)
	# invert C
	invC<-solve(C)
	# compute means
	aX<-colSums(invC%*%X)/sum(invC)
	aY<-colSums(invC%*%Y)/sum(invC)
	# compute cov & cross-cov matrices
	one<-as.matrix(rep(1,n))
	SigXX<-t(X-one%*%aX)%*%invC%*%(X-one%*%aX)/(n-1)
	SigXY<-t(X-one%*%aX)%*%invC%*%(Y-one%*%aY)/(n-1)
	SigYX<-t(SigXY)
	SigYY<-t(Y-one%*%aY)%*%invC%*%(Y-one%*%aY)/(n-1)
	# compute canonical coefficients
	A<-eigen(solve(SigXX)%*%SigXY%*%solve(SigYY)%*%SigYX)
	B<-eigen(solve(SigYY)%*%SigYX%*%solve(SigXX)%*%SigXY)
	# compute canonical variables, rescale
	U<-X%*%A$vectors[,1:min(mX,mY)]
	aU<-colSums(invC%*%U)/sum(invC)
	vcvU<-t(U-one%*%aU)%*%invC%*%(U-one%*%aU)/(n-1)
	U<-(U-one%*%aU)%*%Diag(sqrt(Diag(1/vcvU)))
	V<-Y%*%B$vectors[,1:min(mX,mY)]
	aV<-colSums(invC%*%V)/sum(invC)
	vcvV<-t(V-one%*%aV)%*%invC%*%(V-one%*%aV)/(n-1)
	V<-(V-one%*%aV)%*%Diag(sqrt(Diag(1/vcvV)))
	# compute canonical correlations
	aU<-colSums(invC%*%U)/sum(invC)
	aV<-colSums(invC%*%V)/sum(invC)
	Ccv<-round(t(cbind(U,V))%*%invC%*%cbind(U,V)/(n-1),10)
	ccs<-Diag(Ccv[1:min(mX,mY),(1+min(mX,mY)):(2*min(mX,mY))])
	if(all(Im(ccs)==0)) ccs<-Re(ccs)
	pos<-2*(as.numeric(ccs>0)-0.5)
	ccs<-ccs*pos
	# reorient variables, reorient & rescale coefficents
	U<-U*one%*%pos
	xcoef<-A$vectors[,1:min(mX,mY)]*matrix(1,mX,1)%*%pos%*%Diag(sqrt(Diag(1/vcvU)))
	ycoef<-B$vectors[,1:min(mX,mY)]%*%Diag(sqrt(Diag(1/vcvV)))
	# conduct hypothesis tests
	W_lh<-rep(1,min(mX,mY))
	chiSq<-vector()
	df<-vector()
	for(i in 1:min(mX,mY)){ 
		for(j in i:min(mX,mY)) W_lh[i]<-W_lh[i]*(1-ccs[j]^2)
		chiSq[i]<--((n-1)-(mX+mY+1)/2)*log(W_lh[i])
		df[i]<-(mX+1-i)*(mY+1-i)
	}
	pvalues<-pchisq(chiSq,df=df,lower.tail=F)
	# add row & column names
	if(!is.null(colnames(X))) rownames(xcoef)<-colnames(X)
	if(!is.null(colnames(Y))) rownames(ycoef)<-colnames(Y)
	temp<-vector()
	for(i in 1:min(mX,mY)) temp[i]<-paste("CA",i,sep="")
	colnames(xcoef)<-temp
	colnames(ycoef)<-temp
	colnames(U)<-temp
	colnames(V)<-temp
	# return as list
	lambda<-c(lambda,logL); names(lambda)<-c("lambda","logL")
	obj<-list(cor=ccs,xcoef=xcoef,ycoef=ycoef,xscores=U,yscores=V,lambda=lambda,chisq=chiSq,p=pvalues)
	class(obj)<-"phyl.cca"
	obj
}

## internal function replace Diag
## modified from a suggestion by Travis Ingram 2013
Diag<-function(X){
	if(length(X)==1) X
	else diag(X)
}

## S3 print method
print.phyl.cca<-function(x,digits=6,...){
	cat("\nObject of class \"phyl.cca\" from a phylogenetic canonical")
	cat("\n   correlation analysis.\n\n")
	object<-data.frame(round(x$cor,digits),round(x$chisq,digits),round(x$p,digits))
	colnames(object)<-c("correlation","X^2","P-value")
	rownames(object)<-paste("CC",1:length(x$cor),sep=":")
	cat("Summary of results:\n")
	print(object)
	cat("\nAssumed or estimated value of lambda:\n")
	print(round(x$lambda,digits))
	cat("\nCanonical x coefficients:\n")
	print(as.data.frame(round(x$xcoef,digits)))
	cat("\nCanonical y coefficients:\n")
	print(as.data.frame(round(x$ycoef,digits)))
	cat("\n")
}
