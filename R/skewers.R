## this function does random skewers following Cheverud & Marroig (2007)
## for is.null(method)==FALSE uses clusterGeneration::genPositiveDefMat
## written by Liam J. Revell 2013

skewers<-function(X,Y,nsim=100,method=NULL){
	m<-nrow(X)
	if(!all(sapply(c(dim(X),dim(Y)),"==",m))) 
		stop("X & Y should be square matrices of equal dimension")
	S<-matrix(runif(nsim*m,min=-1,max=1),nsim,m)
	S<-S/matrix(sqrt(rowSums(S^2)),nsim,m)
	Sx<-apply(S,1,"%*%",X)
	Sy<-apply(S,1,"%*%",Y)
	R<-colMeans(Sx*Sy)/sqrt(colMeans(Sx^2)*colMeans(Sy^2))
	r<-mean(R)
	## get null distribution
	foo<-function(m,method,rangeVar){
		if(is.null(method)){
			x<-runif(m,min=-1,max=1)
			x<-x/sqrt(sum(x^2))
			y<-runif(m)
			y<-y/sqrt(sum(y^2)) 
			r<-mean(x*y)/sqrt(mean(x^2)*mean(y^2))
		} else {
			X<-genPositiveDefMat(m,covMethod=method,rangeVar=rangeVar)$Sigma
			Y<-genPositiveDefMat(m,covMethod=method,rangeVar=rangeVar)$Sigma
			S<-matrix(runif(nsim*m,min=-1,max=1),nsim,m)
			S<-S/matrix(sqrt(rowSums(S^2)),nsim,m)
			Sx<-apply(S,1,"%*%",X)
			Sy<-apply(S,1,"%*%",Y)
			R<-colMeans(Sx*Sy)/sqrt(colMeans(Sx^2)*colMeans(Sy^2))
			r<-mean(R)
		}
		r
	}
	Rnull<-replicate(nsim,foo(m,method=method,rangeVar=range(c(diag(X),diag(Y)))))
	return(list(r=r,p=mean(Rnull>=r)))
}