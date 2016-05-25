## function fits Pagel '94 model of correlated evolution of two binary characters
## uses fitMk, ape::ace, or geiger::fitDiscrete internally
## written by Liam J. Revell 2014, 2015, 2016

fitPagel<-function(tree,x,y,method="fitMk",model="ARD",dep.var="xy",...){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	if(dep.var%in%c("x","y","xy")==FALSE){
		cat("  Invalid option for argument \"dep.var\".\n")
		cat("  Setting dep.var=\"xy\" (x depends on y & vice versa)\n\n")
		dep.var<-"xy"
	}
	if(model%in%c("ER","SYM","ARD")==FALSE){
		cat("  Invalid model. Setting model=\"ARD\"\n\n")
		model<-"ARD"
	}
	if(method=="fitDiscrete"){
		chk<-.check.pkg("geiger")
		if(!chk){
			cat("  method = \"fitDiscrete\" requires the package \"geiger\"\n")
			cat("  Defaulting to method = \"fitMk\"\n\n")
			method<-"fitMk"
			fitDiscrete<-function(...) NULL
		}
	}
	if(method%in%c("fitDiscrete","ace","fitMk")==FALSE){
		cat(paste("  method = \"",method,"\" not found.\n",sep=""))
		cat("  Defaulting to method = \"fitMk\"\n\n")
		method<-"fitMk"
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
	## fit independent dep.var
	iQ<-matrix(c(0,1,2,0,3,0,0,2,4,0,0,1,0,4,3,0),4,4,byrow=TRUE)
	if(model%in%c("ER","SYM")) iQ<-make.sym(iQ)
	k.iQ<-length(unique(as.vector(iQ)))-1
	rownames(iQ)<-colnames(iQ)<-levels(xy)
	fit.iQ<-if(method=="fitDiscrete") fitDiscrete(tree,xy,model=iQ,...) 
		else if(method=="ace") ace(xy,tree,type="discrete",model=iQ,...)
		else fitMk(tree,to.matrix(xy,levels(xy)),model=iQ,...)
	## fit dependendent model
	if(dep.var=="xy")
		dQ<-matrix(c(0,1,2,0,3,0,0,4,5,0,0,6,0,7,8,0),4,4,byrow=TRUE)
	else if(dep.var=="x")
		dQ<-matrix(c(0,1,2,0,3,0,0,4,5,0,0,1,0,6,3,0),4,4,byrow=TRUE)
	else if(dep.var=="y")
		dQ<-matrix(c(0,1,2,0,3,0,0,2,4,0,0,5,0,4,6,0),4,4,byrow=TRUE)
	if(model%in%c("ER","SYM")) dQ<-make.sym(dQ)
	k.dQ<-length(unique(as.vector(dQ)))-1
	rownames(dQ)<-colnames(dQ)<-levels(xy)
	fit.dQ<-if(method=="fitDiscrete") fitDiscrete(tree,xy,model=dQ,...) 
		else if(method=="ace") ace(xy,tree,type="discrete",model=dQ,...)
		else fitMk(tree,to.matrix(xy,levels(xy)),model=dQ,...)
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
		independent.AIC=2*k.iQ-2*logLik(fit.iQ),
		dependent.AIC=2*k.dQ-2*logLik(fit.dQ),
		lik.ratio=2*(logLik(fit.dQ)-logLik(fit.iQ)),
		P=pchisq(2*(logLik(fit.dQ)-logLik(fit.iQ)),
		df=k.dQ-k.iQ,
		lower.tail=FALSE),
		method=method,
		dep.var=dep.var,
		model=model)
	class(obj)<-"fitPagel"
	obj
}

## print method for objects of class "fitPagel"
## written by Liam J. Revell 2014, 2016
print.fitPagel<-function(x,...){
	cat("\nPagel's binary character correlation test:\n")
	cat(paste("\nAssumes \"",x$model,
		"\" substitution model for both characters\n",sep=""))
	cat("\nIndependent model rate matrix:\n")
	print(x$independent.Q)
	tmp<-if(x$dep.var=="xy") "x & y" 
		else if(x$dep.var=="x") "x only" 
		else if(x$dep.var=="y") "y only"
	cat(paste("\nDependent (",tmp,") model rate matrix:\n",sep=""))
	print(x$dependent.Q)
	cat("\nModel fit:\n")
	obj<-matrix(c(x$independent.logL,x$dependent.logL,
		x$independent.AIC,x$dependent.AIC),2,2)
	rownames(obj)<-c("independent","dependent")
	colnames(obj)<-c("log-likelihood","AIC")
	print(obj)
	cat("\nHypothesis test result:\n")
	cat(paste("  likelihood-ratio: ",signif(x$lik.ratio,7),"\n"))
	cat(paste("  p-value: ",signif(x$P,7),"\n"))
	cat(paste("\nModel fitting method used was",x$method,"\n\n"))
}

## function borrowed from geiger to pull the Q-matrix from a fit returned by 
## fitDiscrete
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

## make the model matrix symmetric
make.sym<-function(X){
	for(i in 1:nrow(X)) for(j in i:nrow(X)) X[j,i]<-X[i,j]
	X
}