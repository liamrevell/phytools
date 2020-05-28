## function for conditional likelihoods at nodes
## written by Liam J. Revell 2015, 2016, 2019, 2020
## with input from (& structural similarity to) function ace by E. Paradis et al. 2013

fitMk<-function(tree,x,model="SYM",fixedQ=NULL,...){
	if(hasArg(output.liks)) output.liks<-list(...)$output.liks
	else output.liks<-FALSE
	if(hasArg(q.init)) q.init<-list(...)$q.init
	else q.init<-length(unique(x))/sum(tree$edge.length)
	if(hasArg(opt.method)) opt.method<-list(...)$opt.method
	else opt.method<-"nlminb"
	if(hasArg(min.q)) min.q<-list(...)$min.q
	else min.q<-1e-12
	N<-Ntip(tree)
	M<-tree$Nnode
	if(is.matrix(x)){
		x<-x[tree$tip.label,]
		m<-ncol(x)
		states<-colnames(x)
	} else {
		x<-to.matrix(x,sort(unique(x)))
		x<-x[tree$tip.label,]
		m<-ncol(x)
		states<-colnames(x)
	}
	if(hasArg(pi)) pi<-list(...)$pi
	else pi<-"equal"
	if(pi[1]=="equal") pi<-setNames(rep(1/m,m),states)
	else if(pi[1]=="estimated"){ 
		pi<-if(!is.null(fixedQ)) statdist(fixedQ) else statdist(summary(fitMk(tree,x,model),quiet=TRUE)$Q)
		cat("Using pi estimated from the stationary distribution of Q assuming a flat prior.\npi =\n")
		print(round(pi,6))
		cat("\n")
	}	
	else pi<-pi/sum(pi)
	if(is.null(fixedQ)){
		if(is.character(model)){
			rate<-matrix(NA,m,m)
			if(model=="ER"){ 
				k<-rate[]<-1
				diag(rate)<-NA
			} else if(model=="ARD"){
				k<-m*(m-1)
				rate[col(rate)!=row(rate)]<-1:k
			} else if(model=="SYM"){
				k<-m*(m-1)/2
				ii<-col(rate)<row(rate)
				rate[ii]<-1:k
				rate<-t(rate)
				rate[ii]<-1:k
			}
		} else {
			if(ncol(model)!=nrow(model)) 
				stop("model is not a square matrix")
			if(ncol(model)!=ncol(x)) 
				stop("model does not have the right number of columns")
			rate<-model
			k<-max(rate)
		}
		Q<-matrix(0,m,m)
	} else {
		rate<-matrix(NA,m,m)
		k<-m*(m-1)
		rate[col(rate)!=row(rate)]<-1:k
		Q<-fixedQ
	}
	index.matrix<-rate
	tmp<-cbind(1:m,1:m)
	rate[tmp]<-0
	rate[rate==0]<-k+1
	liks<-rbind(x,matrix(0,M,m,dimnames=list(1:M+N,states)))
	pw<-reorder(tree,"pruningwise")
	lik<-function(Q,output.liks=FALSE,pi){
		if(any(is.nan(Q))||any(is.infinite(Q))) return(1e50)
		comp<-vector(length=N+M,mode="numeric")
		parents<-unique(pw$edge[,1])
		root<-min(parents)
		for(i in 1:length(parents)){
			anc<-parents[i]
			ii<-which(pw$edge[,1]==parents[i])
			desc<-pw$edge[ii,2]
			el<-pw$edge.length[ii]
			v<-vector(length=length(desc),mode="list")
			for(j in 1:length(v)){
				v[[j]]<-EXPM(Q*el[j])%*%liks[desc[j],]
			}
			vv<-if(anc==root) Reduce('*',v)[,1]*pi else Reduce('*',v)[,1]
			comp[anc]<-sum(vv)
			liks[anc,]<-vv/comp[anc]
		}
		if(output.liks)return(liks[1:M+N,,drop=FALSE]) 
		logL<--sum(log(comp[1:M+N]))
		return(if(is.na(logL)) Inf else logL)
	}
	if(is.null(fixedQ)){
		if(length(q.init)!=k) q.init<-rep(q.init[1],k)
		if(opt.method=="optim")
			fit<-optim(q.init,function(p) lik(makeQ(m,p,index.matrix),pi=pi),
				method="L-BFGS-B",lower=rep(min.q,k))
		else if(opt.method=="none")
			fit<-list(objective=lik(makeQ(m,q.init,index.matrix),pi=pi),
				par=q.init)
		else	
			fit<-nlminb(q.init,function(p) lik(makeQ(m,p,index.matrix),pi=pi),
				lower=rep(0,k),upper=rep(1e50,k))
		obj<-list(logLik=
			if(opt.method=="optim") -fit$value else -fit$objective,
			rates=fit$par,
			index.matrix=index.matrix,
			states=states,
			pi=pi,
			method=opt.method)
		if(output.liks) obj$lik.anc<-lik(makeQ(m,obj$rates,index.matrix),TRUE,
			pi=pi)
	} else {
		fit<-lik(Q,pi=pi)
		obj<-list(logLik=-fit,
			rates=Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],
			index.matrix=index.matrix,
			states=states,
			pi=pi)
		if(output.liks) obj$lik.anc<-lik(makeQ(m,obj$rates,index.matrix),TRUE,
			pi=pi)
	}
	lik.f<-function(q) -lik(q,output.liks=FALSE,pi=pi)
	obj$lik<-lik.f
	class(obj)<-"fitMk"
	return(obj)
}

makeQ<-function(m,q,index.matrix){
	Q<-matrix(0,m,m)
	Q[]<-c(0,q)[index.matrix+1]
	diag(Q)<-0
	diag(Q)<--rowSums(Q)
	Q
}

## print method for objects of class "fitMk"
print.fitMk<-function(x,digits=6,...){
	cat("Object of class \"fitMk\".\n\n")
	cat("Fitted (or set) value of Q:\n")
	Q<-matrix(NA,length(x$states),length(x$states))
	Q[]<-c(0,x$rates)[x$index.matrix+1]
	diag(Q)<-0
	diag(Q)<--rowSums(Q)
	colnames(Q)<-rownames(Q)<-x$states
	print(round(Q,digits))
	cat("\nFitted (or set) value of pi:\n")
	print(x$pi)
	cat(paste("\nLog-likelihood:",round(x$logLik,digits),"\n"))
	cat(paste("\nOptimization method used was \"",x$method,"\"\n\n",sep=""))
}

## summary method for objects of class "fitMk"
summary.fitMk<-function(object,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-6
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(!quiet) cat("Fitted (or set) value of Q:\n")
	Q<-matrix(NA,length(object$states),length(object$states))
	Q[]<-c(0,object$rates)[object$index.matrix+1]
	diag(Q)<-0
	diag(Q)<--rowSums(Q)
	colnames(Q)<-rownames(Q)<-object$states
	if(!quiet) print(round(Q,digits))
	if(!quiet) cat(paste("\nLog-likelihood:",round(object$logLik,digits),"\n\n"))
	invisible(list(Q=Q,logLik=object$logLik))
}

## logLik method for objects of class "fitMk"
logLik.fitMk<-function(object,...){ 
	lik<-object$logLik
	attr(lik,"df")<-length(object$rates)
	lik
}

## AIC method
AIC.fitMk<-function(object,...,k=2){
	np<-length(object$rates)
	-2*logLik(object)+np*k
}

	
## S3 plot method for objects of class "fitMk"
plot.fitMk<-function(x,...){
	if(hasArg(signif)) signif<-list(...)$signif
	else signif<-3
	if(hasArg(main)) main<-list(...)$main
	else main<-NULL
	if(hasArg(cex.main)) cex.main<-list(...)$cex.main
	else cex.main<-1.2
	if(hasArg(cex.traits)) cex.traits<-list(...)$cex.traits
	else cex.traits<-1
	if(hasArg(cex.rates)) cex.rates<-list(...)$cex.rates
	else cex.rates<-0.6
	if(hasArg(show.zeros)) show.zeros<-list(...)$show.zeros
	else show.zeros<-TRUE
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-6
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-c(1.1,1.1,3.1,1.1)
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-1
	Q<-matrix(NA,length(x$states),length(x$states))
    Q[]<-c(0,x$rates)[x$index.matrix+1]
	diag(Q)<-0
	spacer<-0.1
	plot.new()
	par(mar=mar)
	xylim<-c(-1.2,1.2)
	plot.window(xlim=xylim,ylim=xylim,asp=1)
	if(!is.null(main)) title(main=main,cex.main=cex.main)
	nstates<-length(x$states)
	step<-360/nstates
	angles<-seq(0,360-step,by=step)/180*pi
	if(nstates==2) angles<-angles+pi/2
	v.x<-cos(angles)
	v.y<-sin(angles)
	for(i in 1:nstates) for(j in 1:nstates)
		if(if(!isSymmetric(Q)) i!=j else i>j){
			dx<-v.x[j]-v.x[i]
			dy<-v.y[j]-v.y[i]
			slope<-abs(dy/dx)
			shift.x<-0.02*sin(atan(dy/dx))*sign(j-i)*if(dy/dx>0) 1 else -1
			shift.y<-0.02*cos(atan(dy/dx))*sign(j-i)*if(dy/dx>0) -1 else 1
			s<-c(v.x[i]+spacer*cos(atan(slope))*sign(dx)+
				if(isSymmetric(Q)) 0 else shift.x,
				v.y[i]+spacer*sin(atan(slope))*sign(dy)+
				if(isSymmetric(Q)) 0 else shift.y)
			e<-c(v.x[j]+spacer*cos(atan(slope))*sign(-dx)+
				if(isSymmetric(Q)) 0 else shift.x,
				v.y[j]+spacer*sin(atan(slope))*sign(-dy)+
				if(isSymmetric(Q)) 0 else shift.y)
			if(show.zeros||Q[i,j]>tol){
				if(abs(diff(c(i,j)))==1||abs(diff(c(i,j)))==(nstates-1))
					text(mean(c(s[1],e[1]))+1.5*shift.x,
						mean(c(s[2],e[2]))+1.5*shift.y,
						round(Q[i,j],signif),cex=cex.rates,
						srt=atan(dy/dx)*180/pi)
				else
					text(mean(c(s[1],e[1]))+0.3*diff(c(s[1],e[1]))+
						1.5*shift.x,
						mean(c(s[2],e[2]))+0.3*diff(c(s[2],e[2]))+
						1.5*shift.y,
						round(Q[i,j],signif),cex=cex.rates,
						srt=atan(dy/dx)*180/pi)
				arrows(s[1],s[2],e[1],e[2],length=0.05,
					code=if(isSymmetric(Q)) 3 else 2,lwd=lwd)
			}
		}
	text(v.x,v.y,x$states,cex=cex.traits,
		col=make.transparent("black",0.7))
}

## S3 plot method for objects resulting from fitDiscrete
plot.gfit<-function(x,...){
	if("mkn"%in%class(x$lik)==FALSE){
		stop("Sorry. No plot method presently available for objects of this type.")
	} else {
		obj<-list()
		QQ<-.Qmatrix.from.gfit(x)
		obj$states<-colnames(QQ)
		m<-length(obj$states)
		obj$index.matrix<-matrix(NA,m,m)
		k<-m*(m-1)
		obj$index.matrix[col(obj$index.matrix)!=row(obj$index.matrix)]<-1:k
		obj$rates<-QQ[sapply(1:k,function(x,y) which(x==y),obj$index.matrix)]
		class(obj)<-"fitMk"
		plot(obj,...)
	}
}

## wraps around expm
## written by Liam Revell 2011, 2017
EXPM<-function(x,...){
	e_x<-if(isSymmetric(x)) matexpo(x) else expm(x,...)
	dimnames(e_x)<-dimnames(x)
	e_x
}

## function to simulate multiple-rate Mk multiMk
## written by Liam J. Revell 2018
sim.multiMk<-function(tree,Q,anc=NULL,nsim=1,...){
	if(hasArg(as.list)) as.list<-list(...)$as.list
	else as.list<-FALSE
	ss<-rownames(Q[[1]])
	tt<-map.to.singleton(reorder(tree))
	P<-vector(mode="list",length=nrow(tt$edge))
	for(i in 1:nrow(tt$edge))
		P[[i]]<-expm(Q[[names(tt$edge.length)[i]]]*tt$edge.length[i])
	if(nsim>1) X<- if(as.list) vector(mode="list",length=nsim) else 
		data.frame(row.names=tt$tip.label)
	for(i in 1:nsim){
		a<-if(is.null(anc)) sample(ss,1) else anc
		STATES<-matrix(NA,nrow(tt$edge),2)
		root<-Ntip(tt)+1
		STATES[which(tt$edge[,1]==root),1]<-a
		for(j in 1:nrow(tt$edge)){
			new<-ss[which(rmultinom(1,1,P[[j]][STATES[j,1],])[,1]==1)]
			STATES[j,2]<-new
			ii<-which(tt$edge[,1]==tt$edge[j,2])
			if(length(ii)>0) STATES[ii,1]<-new
		}
		x<-as.factor(
			setNames(sapply(1:Ntip(tt),function(n,S,E) S[which(E==n)],
			S=STATES[,2],E=tt$edge[,2]),tt$tip.label))
		if(nsim>1) X[,i]<-x else X<-x
	}
	X
}

## constant-rate Mk model simulator
## written by Liam J. Revell 2018
sim.Mk<-function(tree,Q,anc=NULL,nsim=1,...){
	if(hasArg(as.list)) as.list<-list(...)$as.list
	else as.list<-FALSE
	ss<-rownames(Q)
	tt<-reorder(tree)
	P<-vector(mode="list",length=nrow(tt$edge))
	for(i in 1:nrow(tt$edge))
		P[[i]]<-expm(Q*tt$edge.length[i])
	if(nsim>1) X<- if(as.list) vector(mode="list",length=nsim) else 
		data.frame(row.names=tt$tip.label)
	for(i in 1:nsim){
		a<-if(is.null(anc)) sample(ss,1) else anc
		STATES<-matrix(NA,nrow(tt$edge),2)
		root<-Ntip(tt)+1
		STATES[which(tt$edge[,1]==root),1]<-a
		for(j in 1:nrow(tt$edge)){
			new<-ss[which(rmultinom(1,1,P[[j]][STATES[j,1],])[,1]==1)]
			STATES[j,2]<-new
			ii<-which(tt$edge[,1]==tt$edge[j,2])
			if(length(ii)>0) STATES[ii,1]<-new
		}
		x<-as.factor(
			setNames(sapply(1:Ntip(tt),function(n,S,E) S[which(E==n)],
			S=STATES[,2],E=tt$edge[,2]),tt$tip.label))
		if(nsim>1) X[[i]]<-x else X<-x
	}
	X
}

## as.Qmatrix method

as.Qmatrix<-function(x,...){
	if(identical(class(x),"Qmatrix")) return(x)
	UseMethod("as.Qmatrix")
}

as.Qmatrix.default<-function(x, ...){
	warning(paste(
		"as.Qmatrix does not know how to handle objects of class ",
		class(x),"."))
}

as.Qmatrix.fitMk<-function(x,...){
	Q<-matrix(NA,length(x$states),length(x$states))
	Q[]<-c(0,x$rates)[x$index.matrix+1]
	rownames(Q)<-colnames(Q)<-x$states
	diag(Q)<-0
	diag(Q)<--rowSums(Q)
	class(Q)<-"Qmatrix"
	Q
}

print.Qmatrix<-function(x,...){
	cat("Estimated Q matrix:\n")
	print(unclass(x))
}
