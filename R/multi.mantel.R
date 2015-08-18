# function for multiple matrix regression with P-values computed by Mantel permutation of the dependent matrix
# written by Liam J. Revell 2012

multi.mantel<-function(Y,X,nperm=1000){
	y<-unfoldLower(Y)
	if(!is.list(X)) X<-list(X)
	Xv<-sapply(X,unfoldLower)
	colnames(Xv)<-paste("X",1:ncol(Xv),sep="")
	fit<-lm(y~Xv)
	coefficients<-fit$coefficients
	smmry<-summary(fit)
	r.squared<-smmry$r.squared
	fstatistic<-smmry$fstatistic[1]
	pF<-0
	tstatistic<-smmry$coefficients[,"t value"]
	pT<-rep(0,length(tstatistic))
	# begin Mantel permutations
	Y<-Yp<-as.matrix(Y)
	for(i in 1:nperm){
		y<-unfoldLower(Yp)
		fitp<-lm(y~Xv)
		smmryp<-summary(fitp)
		pF<-pF+as.numeric(smmryp$fstatistic[1]>=fstatistic)/nperm
		pT<-pT+as.numeric(abs(smmryp$coefficients[,"t value"])>=abs(tstatistic))/nperm
		rndm<-sample(1:nrow(Yp))
		Yp<-Y[rndm,rndm]
	}
	names(coefficients)<-names(tstatistic)<-names(pT)<-c("(intercept)",paste("X",1:ncol(Xv),sep=""))
	names(fstatistic)<-NULL
	residuals<-foldtoLower(fit$residuals); attr(residuals,"Labels")<-rownames(Y)
	fitted.values<-foldtoLower(fit$fitted.values); attr(fitted.values,"Labels")<-rownames(Y)
	return(list(r.squared=r.squared,
		coefficients=coefficients,
		tstatistic=tstatistic,
		fstatistic=fstatistic,
		probt=pT,
		probF=pF,
		residuals=residuals,
		fitted.values=fitted.values))
}

# function unfolds the sub-diagonal of a "dist" object or symmetric matrix into a vector
# written by Liam J. Revell 2012

unfoldLower<-function(X){
	if(class(X)=="dist") X<-as.matrix(X)
	x<-vector()
	for(i in 2:ncol(X)-1) x<-c(x,X[(i+1):nrow(X),i])
	names(x)<-NULL
	return(x)
}

# function folds vector into lower matrix "dist" object

foldtoLower<-function(x){
	n<-(1+sqrt(1+8*length(x)))/2
	X<-matrix(0,n,n)
	j<-0
	for(i in 2:n-1){
		X[(i+1):n,i]<-x[i:(n-1)+j]
		j<-j+(n-i-1)
	}
	return(as.dist(X))
}
