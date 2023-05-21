## function to perform phylogenetic principal components analysis
## multiple morphological traits in Y
## also can use lambda transformation in which lambda is optimized by ML or REML
## written by Liam Revell 2010, 2011, 2013, 2015, 2016, 2017, 2019, 2020, 2022, 2023
## ref. Revell (2009; Evolution)

phyl.pca<-function(tree,Y,method="BM",mode="cov",...){
	## get optional argument
	if(hasArg(opt)) opt<-list(...)$opt
	else opt<-"ML"
	# check tree
	if(!inherits(tree,"phylo")) 
		stop("tree should be an object of class \"phylo\".")
	# check mode
	if(length(strsplit(mode,split="")[[1]])<=2){ 
		message(paste("mode = \"",mode,
			"\" not a valid option; setting mode = \"cov\"",sep=""))
		mode<-"cov"
	}
	if(all(strsplit(mode,split="")[[1]]==strsplit("correlation",
		split="")[[1]][1:length(strsplit(mode,split="")[[1]])])) mode<-"corr"
	else if(all(strsplit(mode,split="")[[1]]==strsplit("covariance",
		split="")[[1]][1:length(strsplit(mode,split="")[[1]])])) mode<-"cov"
	else {
		message(paste("mode = \"",mode,
			"\" not a valid option; setting mode = \"cov\"",sep=""))
		mode="cov"
	}
	if(opt=="REML") object<-reml_phyl.pca(tree,Y,method,mode,...)
	else object<-ml_phyl.pca(tree,Y,method,mode,...)
	object
}

reml_phyl.pca<-function(tree,X,method="BM",mode="cov",...){
	if(!is.binary(tree)) tree<-multi2di(tree)
	lik<-function(lambda,tree,X){
		tt<-lambdaTree(tree,lambda)
		pics<-lapply(X,pic,tt,scaled=FALSE,var.contrasts=TRUE)
		pX<-sapply(pics,function(x) x[,1]/sqrt(x[,2]))
		vcv<-t(pX)%*%pX/(Ntip(tt)-1)
		vars<-pics[[1]][,2]
		logL<-0
		for(i in 1:nrow(pX)){
			x<-sapply(pics,function(x,i) x[i,1],i=i)
			logL<-logL+dmnorm(x,varcov=vars[i]*vcv,
				log=TRUE)
		}
		logL
	}
	X<-X[tree$tip.label,]
	if(method=="lambda"){
		reml_fit<-optimize(lik,c(0,maxLambda(tree)),tree=tree,
			X=as.data.frame(X),maximum=TRUE)
		logL.lambda<-reml_fit$objective
		lambda<-reml_fit$maximum
		tree<-lambdaTree(tree,lambda)
	} else {
		logL.lambda<-lik(1,tree,as.data.frame(X))
		lambda<-1
	}
	pX<-apply(X,2,pic,phy=tree)
	vcv<-t(pX)%*%pX/(Ntip(tree)-1)
	n<-nrow(X)
	m<-ncol(X)
	if(mode=="corr"){
		X=X/matrix(rep(sqrt(diag(vcv)),n),n,m,byrow=TRUE)
		vcv=vcv/(sqrt(diag(vcv))%*%t(sqrt(diag(vcv))))
	}
	a<-apply(X,2,function(x,tree) rep(ace(x,tree,
		method="pic")$ace[1],Ntip(tree)),tree=tree)
	eig<-eigen(vcv)
	Eval<-diag(eig$values)
	colnames(Eval)<-rownames(Eval)<-paste("PC",1:nrow(Eval),sep="")
	Evec<-eig$vectors
	rownames(Evec)<-colnames(X)
	colnames(Evec)<-colnames(Eval)
	S<-as.matrix(X-a)%*%Evec
	pic_corr<-function(x,y,tree){ 
		px<-pic(x,tree)
		py<-pic(y,tree)
		mean(px*py)/sqrt(mean(px^2)*mean(py^2))
	}
	L<-apply(S,2,function(x,y,tree) apply(y,2,pic_corr,
		x=x,tree=tree),y=X,tree=tree)
	dimnames(L)<-dimnames(Evec)
	object<-list(Eval=Eval,Evec=Evec,
		S=S,L=L,lambda=lambda,
		logL.lambda=logL.lambda,
		V=vcv,a=a[1,,drop=FALSE],
		mode=mode,call=match.call())
	class(object)<-"phyl.pca"
	object
}

ml_phyl.pca<-function(tree,Y,method="BM",mode="cov",...){
	## get optional argument
	if(hasArg(opt)) opt<-list(...)$opt
	else opt<-"ML"
	# preliminaries
	n<-nrow(Y)
	m<-ncol(Y)
	# check and sort data
	if(n>Ntip(tree)) 
		stop("number of rows in Y cannot be greater than number of taxa in your tree")
	Y<-as.matrix(Y)
	if(is.null(rownames(Y))){
		if(nrow(Y)==n){ 
			print("Y has no names. function will assume that the row order of Y matches tree$tip.label")
			rownames(Y)<-tree$tip.label
		} else 
			stop("Y has no names and does not have the same number of rows as tips in tree")
	} else if(length(setdiff(rownames(Y),tree$tip.label))!=0) 
		stop("Y has rownames, but some rownames of Y not found in tree")
	# analyze
	C<-vcv.phylo(tree)[rownames(Y),rownames(Y)]
	if(method=="BM"){ 
		temp<-phyl.vcv(Y,C,1.0)
		V<-temp$R
		a<-t(temp$alpha)
		C<-temp$C
	} else if(method=="lambda"){
		if(opt=="ML") temp<-optimize(f=likMlambda,interval=c(0,maxLambda(tree)),X=Y,
			C=C,maximum=TRUE)
		else if(opt=="REML") temp<-optimize(f=remlMlambda,interval=c(0,maxLambda(tree)),
			tree=tree,X=Y,maximum=TRUE)
		else if(opt=="fixed"){
			if(hasArg(lambda)) lambda<-list(...)$lambda
			else {
				cat("  opt=\"fixed\" requires the user to specify lambda.\n")
				cat("  setting lambda to 1.0.\n")
				lambda<-1.0
			}
			temp<-list(maximum=lambda,objective=likMlambda(lambda,X=Y,C=C))
		}	
		lambda<-temp$maximum
		logL<-as.numeric(temp$objective)
		temp<-phyl.vcv(Y,C,lambda)
		V<-temp$R
		a<-t(temp$alpha)
		C<-temp$C
	}
	invC<-solve(C) ## get inverse of C
	# if correlation matrix
	if(mode=="corr"){
		Y=Y/matrix(rep(sqrt(diag(V)),n),n,m,byrow=T) # standardize Y
		V=V/(sqrt(diag(V))%*%t(sqrt(diag(V)))) # change V to correlation matrix
		a<-matrix(colSums(invC%*%Y)/sum(invC),m,1) # recalculate a
	}
	es=eigen(V) # eigenanalyze
	obj<-list()
	obj$Eval<-diag(es$values[1:min(n-1,m)])
	obj$Evec<-es$vectors[,1:min(n-1,m)]
	dimnames(obj$Eval)<-list(paste("PC",1:min(n-1,m),sep=""),
		paste("PC",1:min(n-1,m),sep=""))
	dimnames(obj$Evec)<-list(colnames(Y),paste("PC",1:min(n-1,m),sep=""))
	A<-matrix(rep(a,n),n,m,byrow=T)
	obj$S<-(Y-A)%*%obj$Evec # compute scores in the species space
	Ccv<-t(Y-A)%*%invC%*%obj$S/(n-1) # compute cross covariance matrix and loadings
	obj$L<-matrix(,m,min(n-1,m),dimnames=list(colnames(Y),paste("PC",1:min(n-1,m),sep="")))
	for(i in 1:m) for(j in 1:min(n-1,m)) obj$L[i,j]<-Ccv[i,j]/sqrt(V[i,i]*obj$Eval[j,j])
	if(method=="lambda"){ 
		obj$lambda<-lambda
		obj$logL.lambda<-logL
	}
	obj$V<-temp$R
	obj$a<-a
	obj$mode<-mode
	obj$call<-match.call()
	## assign class attribute (for S3 methods)
	class(obj)<-"phyl.pca"
	# return obj
	obj
}

## S3 method for object of class "phyl.pca
## modified from code provided by Joan Maspons

## S3 print method for "phyl.pca"
## modified from code provided by Joan Maspons
print.phyl.pca<-function(x, ...){
	cat("Phylogenetic pca\n")
	cat("Standard deviations:\n")
	print(sqrt(diag(x$Eval)))
	cat("Loads:\n")
	print(x$L)
	if("lambda" %in% names(x)){
    		cat("lambda:\n")
    		print(x$lambda)
	}
}

## S3 summary method for "phyl.pca"
## modified from code provided by Joan Maspons
summary.phyl.pca<-function(object, ...){
	cat("Importance of components:\n")
	sd<-sqrt(diag(object$Eval))
	varProp<- diag(object$Eval)/sum(object$Eval)
	impp<-rbind("Standard deviation"=sd,"Proportion of Variance"=varProp,
		"Cumulative Proportion"=cumsum(varProp))
	print(impp)
	xx<-list(sdev=sd,importance=impp)
	class(xx)<-"summary.phyl.pca"
	invisible(xx)
}

## S3 biplot method for "phyl.pca"
## modified from code provided by Joan Maspons
## written by Liam J. Revell 2015, 2017
biplot.phyl.pca<-function(x,...){
	to.do<-list(...)
	if(hasArg(choices)){ 
		choices<-list(...)$choices
		to.do$choices<-NULL
	} else choices<-c(1,2)
	to.do$x<-x$S[,choices]
	to.do$y<-x$Evec[,choices]
	do.call(biplot,to.do)
}

## lambdaTree for lambda="REML"
## written by Liam J. Revell 2013
lambdaTree<-function(tree,lambda){
	n<-length(tree$tip.label)
	h1<-nodeHeights(tree)
	ii<-which(tree$edge[,2]>n)
	tree$edge.length[ii]<-lambda*tree$edge.length[ii]
	h2<-nodeHeights(tree)
	tree$edge.length[-ii]<-tree$edge.length[-ii]+h1[-ii,2]-h2[-ii,2]
	tree
}
 
## REML function
## written by Liam J. Revell 2013
remlMlambda<-function(lambda,tree,X){
	tt<-lambdaTree(tree,lambda)
	Y<-apply(X,2,pic,phy=tt)
	V<-t(Y)%*%Y/nrow(Y)
	logL<-sum(dmnorm(Y,mean=rep(0,ncol(Y)),varcov=V,log=TRUE))
	## kronRC<-kronecker(V,diag(rep(1,nrow(Y))))
	## y<-as.vector(Y)
	## logL<--y%*%solve(kronRC,y)/2-length(y)*log(2*pi)/2-determinant(kronRC,logarithm=TRUE)$modulus/2
	## print(V)
	print(c(lambda,logL))
	logL
}

## S3 plot method (does screeplot)
plot.phyl.pca<- function(x,...){
	if(hasArg(main)) main<-list(...)$main
	else main="screeplot"
	x$sdev<-sqrt(diag(x$Eval))
	screeplot(x,main=main)
}

## S3 scores method to extract or compute scores

scores<-function(object,...) UseMethod("scores")

scores.default<-function(object,...){
	warning(paste(
		"scores does not know how to handle objects of class ",
		class(object),"."))
}

scores.phyl.pca<-function(object,...){
	if(hasArg(newdata))newdata<-list(...)$newdata
	else newdata<-NULL
	if(hasArg(dim)) dim<-list(...)$dim
	else dim<-NULL
	if(!is.null(newdata)){
		if(!is.matrix(newdata)) newdata<-as.matrix(newdata)
		if(ncol(newdata)!=nrow(object$Evec))
			stop("Dimensions of newdata incorrect.")
		n<-nrow(newdata)
		m<-ncol(newdata)
		A<-matrix(rep(object$a,n),n,m,byrow=TRUE)
		V<-object$V
		if(object$mode=="corr"){ 
			Y<-newdata/matrix(rep(sqrt(diag(V)),n),n,m,byrow=TRUE)-A
		} else Y<-newdata-A 
		Scores<-Y%*%object$Evec
		if(!is.null(dim)) Scores<-Scores[,dim,drop=FALSE]
	} else { 
		Scores<-if(!is.null(dim)) object$S[,dim,drop=FALSE] else object$S
	}
	Scores
}

## S3 as.princomp method to convert to "princomp" object class

as.princomp<-function(x,...) UseMethod("as.princomp")

as.princomp.default<-function(x,...){
	warning(paste(
		"as.princomp does not know how to handle objects of class ",
		class(x),"."))
}

as.princomp.phyl.pca<-function(x,...){
	nn<-paste("Comp.",1:ncol(x$Evec),sep="")
	obj<-list()
	obj$sdev<-setNames(sqrt(diag(x$Eval)),nn)
	obj$loadings<-x$L
	colnames(obj$loadings)<-nn
	obj$center<-setNames(x$a[1,],rownames(x$Evec))
	obj$scale<-setNames(rep(1,length(obj$center)),names(obj$center))
	obj$scores<-x$S
	colnames(obj$scores)<-nn
	obj$call<-x$call
	class(obj)<-"princomp"
	obj
}