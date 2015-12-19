## function fits Pagel '94 model of correlated evolution of two binary characters
## uses fitMk, ape::ace, or geiger::fitDiscrete internally
## written by Liam J. Revell 2014, 2015

fitPagel<-function(tree,x,y,...){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	if(hasArg(method)) method<-list(...)$method
	else method<-"fitMk"
	if(method=="fitDiscrete"){
		chk<-.check.pkg("geiger")
		if(!chk){
			cat("  method = \"fitDiscrete\" requires the package \"geiger\"\n")
			cat("  Defaulting to method = \"fitMk\"\n\n")
			method<-"fitMk"
			fitDiscrete<-function(...) NULL
		}
	}
	if(!is.factor(x)) x<-as.factor(x)
	levels.x<-levels(x)
	if(!is.factor(y)) y<-as.factor(y)
	levels.y<-levels(y)
	y<-y[names(x)]
	if(length(levels.x)!=2||length(levels.y)!=2)
		stop("Only binary characters for x & y currently permitted.")
	xy<-setNames(factor(paste(x,y,sep="|"),
		levels=sapply(levels.x,paste,levels.y,sep="|")),
		names(x))
	## fit independent model
	iQ<-matrix(c(0,1,2,0,3,0,0,2,4,0,0,1,0,4,3,0),4,4,byrow=TRUE)
	rownames(iQ)<-colnames(iQ)<-levels(xy)
	fit.iQ<-if(method=="fitDiscrete") fitDiscrete(tree,xy,model=iQ) 
		else if(method=="ace") ace(xy,tree,type="discrete",model=iQ)
		else fitMk(tree,xy,model=iQ)
	## fit dependendent model
	dQ<-matrix(c(0,1,2,0,3,0,0,4,5,0,0,6,0,7,8,0),4,4,byrow=TRUE)
	rownames(dQ)<-colnames(dQ)<-levels(xy)
	fit.dQ<-if(method=="fitDiscrete") fitDiscrete(tree,xy,model=dQ) 
		else if(method=="ace") ace(xy,tree,type="discrete",model=dQ)
		else fitMk(tree,xy,model=dQ)
	## back translate independent model
	if(method=="fitDiscrete") iQ<-.Qmatrix.from.gfit(fit.iQ)
	else {
		I<-fit.iQ$index.matrix
		I[I==0]<-NA
		iQ<-apply(I,2,function(i,x) x[i],x=fit.iQ$rates)
		iQ[is.na(iQ)]<-0
		diag(iQ)<--rowSums(iQ)
		rownames(iQ)<-colnames(iQ)
	}
	## dependent model
	if(method=="fitDiscrete") dQ<-.Qmatrix.from.gfit(fit.dQ)
	else {
		I<-fit.dQ$index.matrix
		I[I==0]<-NA
		dQ<-apply(I,2,function(i,x) x[i],x=fit.dQ$rates)
		dQ[is.na(dQ)]<-0
		diag(dQ)<--rowSums(dQ)
		rownames(dQ)<-colnames(dQ)
	}
	## assemble object to return
	obj<-list(independent.Q=iQ,
		dependent.Q=dQ,
		independent.logL=logLik(fit.iQ),
		dependent.logL=logLik(fit.dQ),
		lik.ratio=2*(logLik(fit.dQ)-logLik(fit.iQ)),
		P=pchisq(2*(logLik(fit.dQ)-logLik(fit.iQ)),
		df=length(levels(x))+length(levels(y)),
		lower.tail=FALSE),
		method=method)
	class(obj)<-"fitPagel"
	obj
}

## print method for objects of class "fitPagel"
## written by Liam J. Revell 2014

print.fitPagel<-function(x,...){
	cat("\n  Pagel's binary character correlation test:\n")
	cat("\nIndependent model rate matrix:\n")
	print(x$independent.Q)
	cat("\nDependent model rate matrix:\n")
	print(x$dependent.Q)
	cat("\nModel fit:\n")
	obj<-matrix(c(x$independent.logL,x$dependent.logL),2,1)
	rownames(obj)<-c("independent","dependent")
	colnames(obj)<-"log-likelihood"
	print(obj)
	cat("\nHypothesis test result:\n")
	cat(paste("  likelihood-ratio: ",signif(x$lik.ratio,7),"\n"))
	cat(paste("  p-value: ",signif(x$P,7),"\n"))
	cat(paste("\nModel fitting method used was",x$method,"\n\n"))
}

## function borrowed from geiger to pull the Q-matrix from a fit returned by fitDiscrete

.Qmatrix.from.gfit<-function(x){
	if(!.check.pkg("geiger")) argn<-function(...) NULL
	lik=x$lik
	numberize=function(x){
		y=gsub("q","",x)
		sp=(nn<-nchar(y))/2
		as.numeric(c(substring(y,1,sp),substring(y,sp+1,
			nn)))
	}
	att=attributes(lik)
	att$k=length(att$levels)
	Qmat=matrix(0,att$k,att$k)
	nms=att$argn[att$trns]
	other=att$argn[!att$trns]
	if("constrained"%in%class(lik)){
		cpars=x$opt[argn(lik)]
		apars=names(lik(unlist(cpars),pars.only=TRUE))
		nms=apars[!apars%in%other]
	}
	trns=x$opt[nms]
	for(i in 1:length(trns)){
		nm=names(trns)[i]
		idx=numberize(nm)
		Qmat[idx[1],idx[2]]=trns[[i]]
	}
	diag(Qmat)=-rowSums(Qmat)
	rownames(Qmat)<-colnames(Qmat)<-levels(lik)
	Qmat
}
