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

## S3 plot method for objects of class "fitPagel
## written by Liam J. Revell 2016
plot.fitPagel<-function(x,...){
	if(hasArg(signif)) signif<-list(...)$signif
	else signif<-3
	if(hasArg(main)) main<-list(...)$main
	else main<-NULL
	if(!is.null(main)&&length(main)==1) main<-rep(main,2)
	if(hasArg(cex.main)) cex.main<-list(...)$cex.main
	else cex.main<-1.2
	if(hasArg(cex.sub)) cex.sub<-list(...)$cex.sub
	else cex.sub<-1
	if(hasArg(cex.traits)) cex.traits<-list(...)$cex.traits
	else cex.traits<-0.9
	if(hasArg(cex.rates)) cex.rates<-list(...)$cex.rates
	else cex.rates<-0.8
	if(hasArg(lwd.by.rate)) lwd.by.rate<-list(...)$lwd.by.rate
	else lwd.by.rate<-FALSE
	if(lwd.by.rate){
		rates<-c(x$independent.Q[x$independent.Q>0],x$dependent.Q[x$dependent.Q>0])
		LWD.ind<-round(x$independent.Q/min(rates))
		LWD.dep<-round(x$dependent.Q/min(rates))
		if(hasArg(max.lwd)) max.lwd<-list(...)$max.lwd
		else max.lwd<-10
		LWD.ind[LWD.ind>max.lwd]<-max.lwd
		LWD.dep[LWD.dep>max.lwd]<-max.lwd
	} else LWD.ind<-LWD.dep<-matrix(2,nrow(x$dependent.Q),
		ncol(x$dependent.Q))
	par(mfrow=c(2,1))
	## INDEPENDENT MODEL
	plot.new()
	par(mar=c(1.1,2.1,3.1,2.1))
	plot.window(xlim=c(0,2),ylim=c(0,1),asp=1)
	mtext(if(!is.null(main)) main[1] else "a) Independent model",
		side=3,adj=0,line=1.2,cex=cex.main)
	## trait 1
	text(x=0.15,y=1,"Trait 1:",cex=cex.sub)
	arrows(x0=0.5,y0=0.15,x1=0.5,y1=0.85,
		lwd=max(LWD.ind[3,1],1),
		lty=if(LWD.ind[3,1]==0) "dashed" else "solid",
		length=0.15,lend=3,angle=20)
	arrows(x0=0.55,y0=0.85,x1=0.55,y1=0.15,
		lwd=max(LWD.ind[1,3],1),
		lty=if(LWD.ind[1,3]==0) "dashed" else "solid",
		length=0.15,lend=3,angle=20)
	text(x=0.525,y=0.95,namesplit(rownames(x$dependent.Q)[1])[1],
		cex=cex.traits)
	text(x=0.525,y=0.05,namesplit(rownames(x$dependent.Q)[3])[1],
		cex=cex.traits)
	text(x=0.60,y=0.5,round(x$independent.Q[1,3],signif),cex=cex.rates,srt=90)
	text(x=0.45,y=0.5,round(x$independent.Q[3,1],signif),cex=cex.rates,srt=90)
	## trait 2
	text(x=1.3,y=1,"Trait 2:",cex=cex.sub)	
	arrows(x0=1.65,y0=0.15,x1=1.65,y1=0.85,
		lwd=max(LWD.ind[2,1],1),
		lty=if(LWD.ind[2,1]==0) "dashed" else "solid",
		length=0.15,lend=3,angle=20)
	arrows(x0=1.70,y0=0.85,x1=1.70,y1=0.15,
		lwd=max(LWD.ind[1,2],1),
		lty=if(LWD.ind[1,2]==0) "dashed" else "solid",
		length=0.15,lend=3,angle=20)
	text(x=1.675,y=0.95,namesplit(rownames(x$dependent.Q)[1])[2],
		cex=cex.traits)
	text(x=1.675,y=0.05,namesplit(rownames(x$dependent.Q)[2])[2],
		cex=cex.traits)
	text(x=1.75,y=0.5,round(x$independent.Q[1,2],signif),cex=cex.rates,srt=90)
	text(x=1.60,y=0.5,round(x$independent.Q[2,1],signif),cex=cex.rates,srt=90)
	## DEPENDENT MODEL
	collapse<-
		if(any(sapply(strsplit(rownames(x$dependent.Q),""),length)>6)) 
		",\n" else ", "
	plot.new()
	par(mar=c(1.1,2.1,3.1,2.1))
	plot.window(xlim=c(0,2),ylim=c(0,1),asp=1)
	mtext(if(!is.null(main)) main[2] else "b) Dependent model",
		side=3,adj=0,line=1.2,cex=cex.main)
	text(x=0.15,y=0.95,"Trait 1,\nTrait 2:",cex=cex.sub)
	arrows(x0=0.5,y0=0.15,x1=0.5,y1=0.85,
		lwd=max(LWD.dep[3,1],1),
		lty=if(LWD.dep[3,1]==0) "dashed" else "solid",
		length=0.15,lend=3,angle=20)
	arrows(x0=0.55,y0=0.85,x1=0.55,y1=0.15,
		lwd=max(LWD.dep[1,3],1),
		lty=if(LWD.dep[1,3]==0) "dashed" else "solid",
		length=0.15,lend=3,angle=20)
	arrows(x0=1.45,y0=0.05,x1=0.75,y1=0.05,
		lwd=max(LWD.dep[4,3],1),
		lty=if(LWD.dep[4,3]==0) "dashed" else "solid",
		length=0.15,lend=3,angle=20)
	arrows(x0=0.75,y0=0.1,x1=1.45,y1=0.1,
		lwd=max(LWD.dep[3,4],1),
		lty=if(LWD.dep[3,4]==0) "dashed" else "solid",
		length=0.15,lend=3,angle=20)
	arrows(x0=1.65,y0=0.15,x1=1.65,y1=0.85,
		lwd=max(LWD.dep[4,2],1),
		lty=if(LWD.dep[4,2]==0) "dashed" else "solid",
		length=0.15,lend=3,angle=20)
	arrows(x0=1.7,y0=0.85,x1=1.7,y1=0.15,
		lwd=max(LWD.dep[2,4],1),
		lty=if(LWD.dep[2,4]==0) "dashed" else "solid",
		length=0.15,lend=3,angle=20)
	arrows(x0=1.45,y0=0.9,x1=0.75,y1=0.9,
		lwd=max(LWD.dep[2,1],1),
		lty=if(LWD.dep[2,1]==0) "dashed" else "solid",
		length=0.15,lend=3,angle=20)
	arrows(x0=0.75,y0=0.95,x1=1.45,y1=0.95,
		lwd=max(LWD.dep[1,2],1),
		lty=if(LWD.dep[1,2]==0) "dashed" else "solid",
		length=0.15,lend=3,angle=20)
	## add states
	text(x=0.525,y=0.95,
		paste(namesplit(rownames(x$dependent.Q)[1]),
		collapse=collapse),cex=cex.traits)
	text(x=1.675,y=0.95,
		paste(namesplit(rownames(x$dependent.Q)[2]),
		collapse=collapse),cex=cex.traits)
	text(x=1.675,y=0.05,
		paste(namesplit(rownames(x$dependent.Q)[4]),
		collapse=collapse),cex=cex.traits)
	text(x=0.525,y=0.05,
		paste(namesplit(rownames(x$dependent.Q)[3]),
		collapse=collapse),cex=cex.traits)
	## add rates
	text(x=1.1,y=1,round(x$dependent.Q[1,2],signif),
		cex=cex.rates)
	text(x=1.1,y=0.85,round(x$dependent.Q[2,1],signif),
		cex=cex.rates)
	text(x=1.6,y=0.5,round(x$dependent.Q[4,2],signif),
		cex=cex.rates,srt=90)
	text(x=1.75,y=0.5,round(x$dependent.Q[2,4],signif),
		cex=cex.rates,srt=90)
	text(x=1.1,y=0,round(x$dependent.Q[4,3],signif),
		cex=cex.rates)
	text(x=1.1,y=0.15,round(x$dependent.Q[3,4],signif),
		cex=cex.rates)
	text(x=0.45,y=0.5,round(x$dependent.Q[3,1],signif),
		cex=cex.rates,srt=90)
	text(x=0.6,y=0.5,round(x$dependent.Q[1,3],signif),
		cex=cex.rates,srt=90)
}

namesplit<-function(x){
	tmp<-strsplit(x,"")[[1]]
	ii<-which(tmp=="|")
	c(paste(tmp[1:(ii-1)],collapse=""),
		paste(tmp[(ii+1):length(tmp)],collapse=""))
}