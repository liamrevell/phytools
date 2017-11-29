## function for multiple matrix regression with P-values computed by Mantel permutation of the dependent matrix
## written by Liam J. Revell 2012, 2017

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
	object<-list(r.squared=r.squared,
		coefficients=coefficients,
		tstatistic=tstatistic,
		fstatistic=fstatistic,
		probt=pT,
		probF=pF,
		residuals=residuals,
		fitted.values=fitted.values,
		nperm=nperm)
	class(object)<-"multi.mantel"
	object
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

## S3 methods (added 2017)

print.multi.mantel<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-6
	star<-function(p){
		obj<-if(p>0.1) "" else if(p<=0.1&&p>0.05) "." else
			if(p<=0.05&&p>0.01) "*" else if(p<=0.01&&p>0.001) "**" else
			if(p<=0.001) "***"
		obj
	}
	cat("\nResults from a (multiple) Mantel regression using \"multi.mantel\":\n\n")
	cat("Coefficients:\n")
	object<-data.frame(x$coefficients,
		x$tstatistic,x$probt,
		sapply(x$probt,star))
	rownames(object)<-names(x$coefficients)
	colnames(object)<-c("Estimate","t value","Pr(>|t|)","")
	print(object)
	cat("---\n")
	cat(paste("Signif. codes:  0 \u2018***\u2019 0.001 \u2018**\u2019 0.01", 
		"\u2018*\u2019 0.05 \u2018.\u2019 0.1 \u2018 \u2019 1\n"))
	cat(paste("Pr(>|t|) based on",x$nperm,
		"(Mantel) permutations of rows & columns together in Y.\n\n"))
	cat(paste("Multiple R-squared:",round(x$r.squared,digits),"\n"))
	cat(paste("F-statistic: ",round(x$fstatistic,digits),
		", p-value (based on ",x$nperm," permutations): ",
		round(x$probF,ceiling(log10(x$nperm))),"\n\n",sep=""))
}

residuals.multi.mantel<-function(object,...) object$residuals

fitted.multi.mantel<-function(object,...) object$fitted.values
