## this function computes a phylogenetic reduced major axis (RMA) regression
## written by Liam Revell 2010, 2011, 2012, 2015, 2016

phyl.RMA<-function(x,y,tree,method="BM",lambda=NULL,fixed=FALSE,h0=1.0){
	if(!inherits(tree,"phylo")) 
		stop("tree should be an object of class \"phylo\".")
	x<-x[tree$tip.label]; y<-y[tree$tip.label]
	# bind the x & y into columns
	X<-cbind(x,y)
	if(method=="lambda")
		if(fixed==FALSE)
			result<-optimize(f=likMlambda,interval=c(0,1),X=X,C=vcv(tree),
				maximum=TRUE)
		else
			result<-list(objective=likMlambda(lambda,X,vcv(tree)),maximum=lambda)
	else if(method=="BM")
		result<-list(objective=likMlambda(1.0,X,vcv(tree)),maximum=1.0)
	else
		stop("do not recognize method")	
	est.lambda<-result$maximum # estimated lambda
	C<-vcv(tree)
	C<-lambda.transform(est.lambda,C)
	temp<-phyl.vcv(X,vcv(tree),lambda=est.lambda)
	beta1<-sign(temp$R[1,2])*sqrt(temp$R[2,2]/temp$R[1,1])
	beta0<-temp$a[2]-beta1*temp$a[1]
	r<-y-(beta0+beta1*x)
	r2<-temp$R[1,2]^2/(temp$R[1,1]*temp$R[2,2])
	if(sign(beta1)!=sign(h0)){
		warning("b & h0 have different signs; hypothesis test invalid")
		T<-0
	} else
		T<-abs(log(abs(beta1))-log(abs(h0)))/sqrt((1-r2)/(Ntip(tree)-2))
	df<-2+(Ntip(tree)-2)/(1+0.5*r2)
	P<-2*pt(T,df=df,lower.tail=FALSE)
	test<-c(r2,T,df,P); names(test)<-c("r2","T","df","P")
	object<-list(RMA.beta=c(beta0,beta1),V=temp$R,lambda=est.lambda,
		logL=as.numeric(result$objective),test=test,h0=h0,
		resid=as.matrix(r))
	class(object)<-"phyl.RMA"
	object
}

print.phyl.RMA<-function(x,...){
	cat("\nCoefficients:\n")
	print(coef(x))
	cat("\nVCV matrix:\n")
	print(x$V)
	cat("\n")
	if(is.null(x$lambda))
		cat("Model for the covariance structure of the error is \"BM\"\n\n")
	else {
		cat("Model for the covariance structure of the error is \"lambda\"\n")
		cat("\nEstimate(s):\n")
		print(setNames(c(x$lambda,x$logL),c("lambda","log(L)")))
		cat("\n")
	}
	cat("Hypothesis test based on Clarke (1980; Biometrika):\n")
	print(round(x$test,6))
	cat(paste("\nNote that the null hypothesis test is h0 =",x$h0,
		"\n\n"))
}

coef.phyl.RMA<-function(object,...){
	val<-setNames(object$RMA.beta,c("(Intercept)","x"))
	val
}

residuals.phyl.RMA<-function(object,...) object$resid[,1]
