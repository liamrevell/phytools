## optimizing, graphing, and analyzing extended Mk model for discrete
## character evolution
## written by Liam J. Revell (updates in 2015, 2016, 2019, 2020, 2021, 2022, 2023)
## likelihood function (with pruning) adapted from ape::ace (Paradis et al. 2013)
## lik.func="pruning" uses phytools::pruning to compute likelihood instead

## function to simulate multiple-rate Mk multiMk
## written by Liam J. Revell 2018
sim.multiMk<-function(tree,Q,anc=NULL,nsim=1,...){
	if(hasArg(as.list)) as.list<-list(...)$as.list
	else as.list<-FALSE
	if(hasArg(internal)) internal<-list(...)$internal
	else internal<-FALSE
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
		if(internal){
			x<-as.factor(setNames(sapply(1:(Ntip(tt)+tt$Nnode),
				function(n,S,E) S[which(E==n)[1]],S=STATES,E=tt$edge),
				c(tt$tip.label,1:tt$Nnode+Ntip(tt))))
		} else{
			x<-as.factor(
				setNames(sapply(1:Ntip(tt),function(n,S,E) S[which(E==n)],
				S=STATES[,2],E=tt$edge[,2]),tt$tip.label))
		}
		if(nsim>1) X[,i]<-x else X<-x
	}
	X
}

## constant-rate Mk model simulator
## written by Liam J. Revell 2018, 2023

sim.Mk<-function(tree,Q,anc=NULL,nsim=1,...){
	if(hasArg(as.list)) as.list<-list(...)$as.list
	else as.list<-FALSE
	if(hasArg(internal)) internal<-list(...)$internal
	else internal<-FALSE
	ss<-rownames(Q)
	tt<-reorder(tree)
	P<-vector(mode="list",length=nrow(tt$edge))
	for(i in 1:nrow(tt$edge))
		P[[i]]<-expm(Q*tt$edge.length[i])
	if(nsim>1) X<- if(as.list) vector(mode="list",length=nsim) else 
		data.frame(row.names=tt$tip.label)
	for(i in 1:nsim){
		if(is.null(anc)) a<-sample(ss,1)
		else if(is.numeric(anc)) a<-sample(names(anc),1,prob=anc)
		else a<-anc
		STATES<-matrix(NA,nrow(tt$edge),2)
		root<-Ntip(tt)+1
		STATES[which(tt$edge[,1]==root),1]<-a
		for(j in 1:nrow(tt$edge)){
			new<-ss[which(rmultinom(1,1,P[[j]][STATES[j,1],])[,1]==1)]
			STATES[j,2]<-new
			ii<-which(tt$edge[,1]==tt$edge[j,2])
			if(length(ii)>0) STATES[ii,1]<-new
		}
		if(internal){
			x<-as.factor(setNames(sapply(1:(Ntip(tt)+tt$Nnode),
				function(n,S,E) S[which(E==n)[1]],S=STATES,E=tt$edge),
				c(tt$tip.label,1:tt$Nnode+Ntip(tt))))
		} else{
			x<-as.factor(
				setNames(sapply(1:Ntip(tt),function(n,S,E) S[which(E==n)],
				S=STATES[,2],E=tt$edge[,2]),tt$tip.label))
		}
		if(nsim>1) X[[i]]<-x else X<-x
	}
	X
}

anova.fitMk<-function(object,...){
	fits<-list(...)
	nm<-c(
		deparse(substitute(object)),
		if(length(fits)>0) sapply(substitute(list(...))[-1],deparse)
	)
	logL<-c(logLik(object),
		if(length(fits)>0) sapply(fits,logLik))
	df<-c(attr(logLik(object),"df"),
		if(length(fits)>0) sapply(fits,function(x) attr(logLik(x),"df")))
	AICvals<-c(AIC(object),
		if(length(fits)>0) sapply(fits,AIC))
	ww<-aic.w(AICvals)
	result<-data.frame(logL,df,AICvals,unclass(ww))
	rownames(result)<-nm
	colnames(result)<-c("log(L)","d.f.","AIC","weight")
	models<-c(list(object),fits)
	attr(result,"models")<-models
	class(result)<-c(class(result),"anova.fitMk")
	print(result)
	invisible(result)
}

fitMk<-function(tree,x,model="SYM",fixedQ=NULL,...){
	if(hasArg(opt.method)) opt.method<-list(...)$opt.method
	else opt.method<-"nlminb"
	if(hasArg(lik.func)) lik.func<-list(...)$lik.func
	else lik.func<-"lik"
	if(lik.func=="parallel"){
		if(hasArg(ncores)) ncores<-list(...)$ncores
		else ncores<-min(nrow(tree$edge),detectCores()-1)
	}	
	if(opt.method=="optimParallel"){ 
		if(hasArg(ncores)) ncores<-list(...)$ncores
		else ncores<-detectCores()
		if(is.na(ncores)) ncores<-1 
		args<-list(...)
		args$tree<-tree
		args$x<-x
		args$model<-model
		args$ncores<-ncores
		obj<-do.call(fitMk.parallel,args)
	} else {
		if(hasArg(output.liks)) output.liks<-list(...)$output.liks
		else output.liks<-FALSE
		if(hasArg(smart_start)) smart_start<-list(...)$smart_start
		else smart_start<-FALSE
		if(hasArg(q.init)) q.init<-list(...)$q.init
		else q.init<-length(unique(x))/sum(tree$edge.length)
		if(hasArg(rand_start)) rand_start<-list(...)$rand_start
		else rand_start<-FALSE
		if(hasArg(min.q)) min.q<-list(...)$min.q
		else min.q<-1e-12
		if(hasArg(max.q)) max.q<-list(...)$max.q
		else max.q<-max(nodeHeights(tree))*100
		if(hasArg(logscale)) logscale<-list(...)$logscale
		else logscale<-FALSE
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
		if(is.numeric(pi)) root.prior<-"given"
		if(pi[1]=="equal"){ 
			pi<-setNames(rep(1/m,m),states)
			root.prior<-"flat"
		} else if(pi[1]=="estimated"){ 
			pi<-if(!is.null(fixedQ)) statdist(fixedQ) else 
				statdist(summary(fitMk(tree,x,model),quiet=TRUE)$Q)
			cat(paste("Using pi estimated from the stationary",
				"distribution of Q assuming a flat prior.\npi =\n"))
			print(round(pi,6))
			cat("\n")
			root.prior<-"stationary"
		} else if(pi[1]=="fitzjohn"){ 
			root.prior<-"nuisance"
		} else if(pi[1]=="mle"){ 
			root.prior<-"it's MLE"
			if(lik.func=="lik") lik.func<-"pruning"
		}
		if(is.numeric(pi)){ 
			pi<-pi/sum(pi)
			if(is.null(names(pi))) pi<-setNames(pi,states)
			pi<-pi[states]
		} 
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
		if(lik.func%in%c("pruning","parallel")){
			MODEL<-rate
			MODEL[is.na(MODEL)]<-0
			diag(MODEL)<-0
			if(hasArg(expm.method)) expm.method<-list(...)$expm.method
			else expm.method<-"Higham08.b"
		} else expm.method<-"Higham08.b"
		tmp<-cbind(1:m,1:m)
		rate[tmp]<-0
		rate[rate==0]<-k+1
		liks<-rbind(x,matrix(0,M,m,dimnames=list(1:M+N,states)))
		pw<-reorder(tree,"postorder")
		lik<-function(Q,output.liks=FALSE,pi,...){
			if(hasArg(output.pi)) output.pi<-list(...)$output.pi
			else output.pi<-FALSE
			if(is.Qmatrix(Q)) Q<-unclass(Q)
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
				if(anc==root){
					if(is.numeric(pi)) vv<-Reduce('*',v)[,1]*pi
					else if(pi[1]=="fitzjohn"){
						D<-Reduce('*',v)[,1]
						pi<-D/sum(D)
						vv<-D*D/sum(D)
					}
				} else vv<-Reduce('*',v)[,1]
				comp[anc]<-sum(vv)
				liks[anc,]<-vv/comp[anc]
			}
			if(output.liks) return(liks[1:M+N,,drop=FALSE])
			else if(output.pi) return(pi)
			else {
				logL<--sum(log(comp[1:M+N]))
				if(is.na(logL)) logL<-Inf
				return(logL)
			}
		}
		if(is.null(fixedQ)){
			if(smart_start&&max(index.matrix,na.rm=TRUE)>1){
				MM<-index.matrix
				MM[is.na(MM)]<-0
				MM[MM>0]<-1
				q.init<-fitMk(pw,x,model=MM,pi=pi,opt.method="nlminb")$rates
			}
			if(length(q.init)!=k) q.init<-rep(q.init[1],k)
			if(rand_start) q.init<-q.init*rexp(length(q.init),1)
			q.init<-if(logscale) log(q.init) else q.init
			if(opt.method=="optim"){
				if(lik.func=="lik"){
					fit<-if(logscale) 
						optim(q.init,function(p) lik(makeQ(m,exp(p),index.matrix),pi=pi),
							method="L-BFGS-B",lower=rep(log(min.q),k),upper=rep(log(max.q),k)) else
						optim(q.init,function(p) lik(makeQ(m,p,index.matrix),pi=pi),
							method="L-BFGS-B",lower=rep(min.q,k),upper=rep(max.q,k))
				} else if(lik.func=="pruning") {
					fit<-if(logscale)
						optim(q.init,function(p) -pruning(exp(p),tree=pw,x=x,model=MODEL,pi=pi,
							expm.method=expm.method),method="L-BFGS-B",lower=rep(log(min.q),k),
							upper=rep(log(max.q),k)) else
						optim(q.init,function(p) -pruning(p,tree=pw,x=x,model=MODEL,pi=pi,
							expm.method=expm.method),method="L-BFGS-B",lower=rep(min.q,k),
							upper=rep(max.q,k))
				} else if(lik.func=="parallel") {
					mc<-makeCluster(ncores,type="PSOCK")
					registerDoParallel(cl=mc)
					fit<-if(logscale)
						optim(q.init,function(p) -parallel_pruning(exp(p),tree=pw,x=x,model=MODEL,
							pi=pi,expm.method=expm.method),method="L-BFGS-B",lower=rep(log(min.q),k),
							upper=rep(log(max.q),k)) else
						optim(q.init,function(p) -parallel_pruning(p,tree=pw,x=x,model=MODEL,pi=pi,
							expm.method=expm.method),method="L-BFGS-B",lower=rep(min.q,k),
							upper=rep(max.q,k))
				}
			} else if(opt.method=="none"){
				if(lik.func=="lik")
					fit<-list(objective=lik(makeQ(m,q.init,index.matrix),pi=pi),
						par=q.init)
				else if(lik.func=="pruning")
					fit<-list(objective=-pruning(q.init,pw,x,MODEL,pi=pi,expm.method=expm.method),
						par=q.init)
				else if(lik.func=="parallel"){
					mc<-makeCluster(ncores,type="PSOCK")
					registerDoParallel(cl=mc)
					fit<-list(objective=-parallel_pruning(q.init,pw,x,MODEL,pi=pi,
						expm.method=expm.method),par=q.init)
				}
			} else {
				if(lik.func=="lik"){
					fit<-if(logscale)
						nlminb(q.init,function(p) lik(makeQ(m,exp(p),index.matrix),pi=pi),
							lower=rep(log(min.q),k),upper=rep(log(max.q),k)) else 
						nlminb(q.init,function(p) lik(makeQ(m,p,index.matrix),
							pi=pi),lower=rep(0,k),upper=rep(max.q,k))
				} else if(lik.func=="pruning"){
					fit<-if(logscale)
						nlminb(q.init,function(p) -pruning(exp(p),tree=pw,x=x,model=MODEL,
							pi=pi,expm.method=expm.method),lower=rep(log(min.q),k),
							upper=rep(log(max.q),k)) else
						nlminb(q.init,function(p) -pruning(p,tree=pw,x=x,model=MODEL,
							pi=pi,expm.method=expm.method),lower=rep(0,k),upper=rep(max.q,k))
				} else if(lik.func=="parallel"){
					mc<-makeCluster(ncores,type="PSOCK")
					registerDoParallel(cl=mc)
					fit<-if(logscale)
						nlminb(q.init,function(p) -parallel_pruning(exp(p),tree=pw,x=x,
							model=MODEL,pi=pi,expm.method=expm.method),lower=rep(log(min.q),k),
							upper=rep(log(max.q),k)) else
						nlminb(q.init,function(p) -parallel_pruning(p,tree=pw,x=x,model=MODEL,
							pi=pi,expm.method=expm.method),lower=rep(0,k),upper=rep(max.q,k))
				}
			}
			if(logscale) fit$par<-exp(fit$par)
			if(pi[1]=="fitzjohn") pi<-if(lik.func=="parallel") setNames(
				parallel_pruning(fit$par,tree=pw,x=x,model=MODEL,pi="fitzjohn",
				expm.method=expm.method,return="pi"),states) else setNames(
				lik(makeQ(m,fit$par,index.matrix),FALSE,pi=pi,output.pi=TRUE),
				states)
			else if(pi[1]=="mle") pi<-if(lik.func=="parallel") setNames(
				parallel_pruning(fit$par,tree=pw,x=x,model=MODEL,pi="mle",
				expm.method=expm.method,return="pi"),states) else setNames(
				pruning(fit$par,tree=pw,x=x,model=MODEL,pi="mle",expm.method=expm.method,
				return="pi"),states)
			obj<-list(logLik=
				if(opt.method=="optim") -fit$value else -fit$objective,
				rates=fit$par,
				index.matrix=index.matrix,
				states=states,
				pi=pi,
				method=opt.method,
				root.prior=root.prior)
			if(opt.method=="nlminb")
				obj$opt_results<-fit[c("convergence","iterations","evaluations","message")]
			else if(opt.method=="optim")
				obj$opt_results<-fit[c("counts","convergence","message")]
			if(output.liks) obj$lik.anc<-lik(makeQ(m,obj$rates,index.matrix),TRUE,
				pi=pi)
		} else {
			if(lik.func=="lik"){
				fit<--lik(Q,pi=pi)
				if(pi[1]=="fitzjohn") pi<-setNames(lik(Q,FALSE,pi=pi,output.pi=TRUE),
					states)
			} else if(lik.func%in%c("pruning","parallel")){
				q<-Q[sapply(1:k,function(x,y) which(x==y),index.matrix)]
				if(lik.func=="pruning"){
					fit<-pruning(q,pw,x,model=MODEL,expm.method=expm.method,pi=pi)
					if(pi[1]=="fitzjohn") pi<-setNames(pruning(q,tree=pw,x=x,model=MODEL,
						pi="fitzjohn",expm.method=expm.method,return="pi"),states)
					else if(pi[1]=="mle") pi<-setNames(pruning(q,tree=pw,x=x,model=MODEL,
						pi="mle",expm.method=expm.method,return="pi"),states)
				} else if(lik.func=="parallel"){
					mc<-makeCluster(ncores,type="PSOCK")
					registerDoParallel(cl=mc)
					fit<-parallel_pruning(q,pw,x,model=MODEL,expm.method=expm.method,pi=pi)
					if(pi[1]=="fitzjohn") pi<-setNames(parallel_pruning(q,tree=pw,x=x,
						model=MODEL,pi="fitzjohn",expm.method=expm.method,return="pi"),
						states)
					else if(pi[1]=="mle") pi<-setNames(parallel_pruning(q,tree=pw,x=x,
						model=MODEL,pi="mle",expm.method=expm.method,return="pi"),states)
				}
			}
			obj<-list(logLik=fit,
				rates=Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],
				index.matrix=index.matrix,
				states=states,
				pi=pi,
				root.prior=root.prior)
			if(output.liks) obj$lik.anc<-lik(makeQ(m,obj$rates,index.matrix),TRUE,
				pi=pi)
		}
		if(lik.func=="lik")
			lik.f<-function(q) -lik(q,output.liks=FALSE,
				pi=if(root.prior=="nuisance") "fitzjohn" else pi)
		else if(lik.func%in%c("pruning","parallel")) {
			lik.f<-function(q){
				q<-sapply(1:max(MODEL), function(ind,q,MODEL) q[which(MODEL==ind)],
					q=q,MODEL=MODEL)
				if(lik.func=="pruning"){
					pruning(q,tree=pw,x=x,model=MODEL,
						pi=if(root.prior=="nuisance") "fitzjohn" else pi,
						expm.method=expm.method)
				} else if(lik.func=="parallel"){
					parallel_pruning(q,tree=pw,x=x,model=MODEL,
						pi=if(root.prior=="nuisance") "fitzjohn" else pi,
						expm.method=expm.method)
				}
			}
		}
		obj$data<-x
		obj$tree<-tree
		obj$lik<-lik.f
		class(obj)<-"fitMk"
		if(lik.func=="parallel") stopCluster(cl=mc)
	}
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
	print(round(x$pi,digits))
	cat(paste("due to treating the root prior as (a) ",x$root.prior,".\n",
		sep=""))
	cat(paste("\nLog-likelihood:",round(x$logLik,digits),"\n"))
	cat(paste("\nOptimization method used was \"",x$method,"\"\n\n",
		sep=""))
	if(!is.null(x$opt_results$convergence)){
		if(x$opt_results$convergence==0) 
			cat("R thinks it has found the ML solution.\n\n")
		else cat("R thinks optimization may not have converged.\n\n")
	}
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
	if(!is.null(object$index.matrix)) 
		attr(lik,"df")<-max(object$index.matrix,na.rm=TRUE)
	else
		attr(lik,"df")<-length(object$rates)
	lik
}

## S3 plot method for objects of class "fitMk"
plot.fitMk<-function(x,...){
	Q<-as.Qmatrix(x)
	plot(Q,...)
}

## S3 plot method for "gfit" object from geiger::fitDiscrete
plot.gfit<-function(x,...){
	if("mkn"%in%class(x$lik)==FALSE){
		stop("Sorry. No plot method presently available for objects of this type.")
		object<-NULL
	} else {
		chk<-.check.pkg("geiger")
		if(chk) object<-plot(as.Qmatrix(x),...)
		else {
			obj<-list()
			QQ<-.Qmatrix.from.gfit(x)
			obj$states<-colnames(QQ)
			m<-length(obj$states)
			obj$index.matrix<-matrix(NA,m,m)
			k<-m*(m-1)
			obj$index.matrix[col(obj$index.matrix)!=row(obj$index.matrix)]<-1:k
			obj$rates<-QQ[sapply(1:k,function(x,y) which(x==y),obj$index.matrix)]
			class(obj)<-"fitMk"
			object<-plot(obj,...)
		}
	}
	invisible(object)
}

MIN<-function(x,...) min(x[is.finite(x)],...)
MAX<-function(x,...) max(x[is.finite(x)],...)
RANGE<-function(x,...) range(x[is.finite(x)],...)
	
## S3 method for "Qmatrix" object class
plot.Qmatrix<-function(x,...){
	Q<-unclass(x)
	if(hasArg(asp)) asp<-list(...)$asp
	else asp<-1
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
	if(hasArg(umbral)) umbral<-list(...)$umbral	
	else umbral<-FALSE
	if(hasArg(ncat)) ncat<-list(...)$ncat
	else ncat<-NULL
	if(hasArg(spacer)) spacer<-list(...)$spacer
	else spacer<-0.1
	if(hasArg(color)) color<-list(...)$color
	else color<-FALSE
	if(hasArg(width)) width<-list(...)$width
	else width<-FALSE
	if(hasArg(text)) text<-list(...)$text
	else text<-TRUE
	if(hasArg(max.lwd)) max.lwd<-list(...)$max.lwd
	else max.lwd<-if(text) 5 else 8
	if(hasArg(rotate)) rotate<-list(...)$rotate
	else rotate<-NULL
	if(hasArg(add)) add<-list(...)$add
	else add<-FALSE
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-NULL
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-NULL
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-0.02
	if(hasArg(palette)) palette<-list(...)$palette
	else palette<-c("blue","purple","red")
	## set all Q<tol to zero (may remove later)
	Q[Q<tol]<-0
	## end may remove later
	if(!add) plot.new()
	par(mar=mar)
	if(is.null(xlim)) xlim<-ylim
	if(is.null(ylim)) ylim<-xlim
	if(is.null(xlim)&&is.null(ylim)){
		if(!color) xlim<-ylim<-c(-1.2,1.2)
		else { 
			xlim<-c(-1.4,1)
			ylim<-c(-1.2,1.2)
		}
	}
	plot.window(xlim=xlim,ylim=ylim,asp=asp)
	if(!is.null(main)) title(main=main,cex.main=cex.main)
	nstates<-nrow(Q)
	if(is.null(rotate)){
		if(nstates==2) rotate<--90
		else rotate<--90*(nstates-2)/(nstates)
	}
	if(color){
		col_pal<-function(qq) if(is.na(qq)) NA else 
			if(is.infinite(qq)) make.transparent("grey",0.6) else
			rgb(colorRamp(palette)(qq),maxColorValue=255)
		qq<-Q
		diag(qq)<-NA
		qq<-log(qq)
		dq<-diff(RANGE(qq,na.rm=TRUE))
		if(dq<tol){
			cols<-matrix(palette[1],nstates,nstates)
			cols[Q<tol]<-make.transparent("grey",0.6)
		} else {
			qq<-(qq-MIN(qq,na.rm=TRUE))/dq
			cols<-apply(qq,c(1,2),col_pal)
		}
	} else cols<-matrix(par("fg"),nstates,nstates)
	if(width){
		lwd_maker<-function(qq,max.qq) if(is.na(qq)) NA else 
			if(is.infinite(qq)) 0 else qq*(max.lwd-1)+1
		qq<-Q
		diag(qq)<-NA
		qq<-log(qq)
		dq<-max(qq[!is.infinite(qq)],na.rm=TRUE)-
			min(qq[!is.infinite(qq)],na.rm=TRUE)
		if(dq<tol){
			lwd<-matrix(lwd,nstates,nstates)
			lwd[Q<tol]<-0
		} else {
			qq<-(qq-(min(qq[!is.infinite(qq)],na.rm=TRUE)))/dq
			lwd<-apply(qq,c(1,2),lwd_maker,max.qq=max(qq,na.rm=TRUE))
		}
	} else lwd<-matrix(lwd,nstates,nstates)
	if(!umbral||is.null(ncat)){
		step<-360/nstates
		angles<-seq(rotate,360-step+rotate,by=step)/180*pi
		if(nstates==2) angles<-angles+pi/2
		v.x<-cos(angles)
		v.y<-sin(angles)
	} else {
		v.x<-v.y<-vector()
		for(i in 1:length(ncat)){
			Q<-Q[sort(rownames(Q)),sort(colnames(Q))]
			xp<--1+2*(i-1)/(length(ncat)-1)
			v.x<-c(v.x,rep(xp,ncat[i]))
			yp<-seq(1,-1,length.out=max(ncat))[1:ncat[i]]
			v.y<-c(v.y,yp)
		}
	}	
	for(i in 1:nstates) for(j in 1:nstates)
		if(if(!isSymmetric(Q)) i!=j else i>j){
			dx<-v.x[j]-v.x[i]
			dy<-v.y[j]-v.y[i]
			slope<-abs(dy/dx)
			shift.x<-offset*sin(atan(dy/dx))*sign(j-i)*if(dy/dx>0) 1 else -1
			shift.y<-offset*cos(atan(dy/dx))*sign(j-i)*if(dy/dx>0) -1 else 1
			s<-c(v.x[i]+spacer*cos(atan(slope))*sign(dx)+
				if(isSymmetric(Q)) 0 else shift.x,
				v.y[i]+spacer*sin(atan(slope))*sign(dy)+
				if(isSymmetric(Q)) 0 else shift.y)
			e<-c(v.x[j]+spacer*cos(atan(slope))*sign(-dx)+
				if(isSymmetric(Q)) 0 else shift.x,
				v.y[j]+spacer*sin(atan(slope))*sign(-dy)+
				if(isSymmetric(Q)) 0 else shift.y)
			if(show.zeros||Q[i,j]>tol){
				if(text){
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
				}
				arrows(s[1],s[2],e[1],e[2],length=0.05,
					code=if(isSymmetric(Q)) 3 else 2,
					lwd=if(lwd[i,j]==0) 1 else lwd[i,j],
					lty=if(lwd[i,j]==0) "dotted" else "solid",
					col=cols[i,j])
			}
		}
	text(v.x,v.y,rownames(Q),cex=cex.traits,
		col=make.transparent(par("fg"),0.9))
	if(color){
		if(dq>tol){
			h<-1.5
			LWD<-diff(par()$usr[1:2])/dev.size("px")[1]
			lines(x=rep(0.93*xlim[1]+LWD*15/2,2),y=c(-h/2,h/2))
			nticks<-6
			Y<-cbind(seq(-h/2,h/2,length.out=nticks),
				seq(-h/2,h/2,length.out=nticks))
			X<-cbind(rep(0.93*xlim[1]+LWD*15/2,nticks),
				rep(0.93*xlim[1]+LWD*15/2+0.02*h,nticks))
			for(i in 1:nrow(Y)) lines(X[i,],Y[i,])
			add.color.bar(h,sapply(seq(0,1,length.out=100),col_pal),
				title="evolutionary rate (q)",
				lims=NULL,digits=3,
				direction="upwards",
				subtitle="",lwd=15,
				x=0.93*xlim[1],y=-h/2,prompt=FALSE)
			QQ<-Q
			diag(QQ)<-0
			text(x=X[,2],y=Y[,2],signif(exp(seq(MIN(log(QQ),na.rm=TRUE),
				MAX(log(QQ),na.rm=TRUE),length.out=6)),signif),pos=4,cex=0.7)
		} else {
			BLUE<-function(...) palette[1]
			h<-1.5
			LWD<-diff(par()$usr[1:2])/dev.size("px")[1]
			lines(x=rep(0.93*xlim[1]+LWD*15/2,2),y=c(-h/2,h/2))
			nticks<-6
			Y<-cbind(seq(-h/2,h/2,length.out=nticks),
				seq(-h/2,h/2,length.out=nticks))[nticks,,drop=FALSE]
			X<-cbind(rep(0.93*xlim[1]+LWD*15/2,nticks),
				rep(0.93*xlim[1]+LWD*15/2+0.02*h,nticks))[nticks,,drop=FALSE]
			for(i in 1:nrow(Y)) lines(X[i,],Y[i,])
			add.color.bar(h,sapply(seq(0,1,length.out=100),BLUE),
				title="evolutionary rate (q)",
				lims=NULL,digits=3,
				direction="upwards",
				subtitle="",lwd=15,
				x=0.93*xlim[1],y=-h/2,prompt=FALSE)
			QQ<-Q
			diag(QQ)<-0
			text(x=X[,2],y=Y[,2],signif(exp(seq(MIN(log(QQ),na.rm=TRUE),
				MAX(log(QQ),na.rm=TRUE),length.out=1)),signif),pos=4,cex=0.7)
		}
	}
	object<-data.frame(states=rownames(Q),x=v.x,y=v.y)
	invisible(object)
}

## wraps around expm
## written by Liam Revell 2011, 2017
EXPM<-function(x,...){
	e_x<-if(isSymmetric(x)) matexpo(x) else expm(x,...)
	dimnames(e_x)<-dimnames(x)
	e_x
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

as.Qmatrix.matrix<-function(x, ...){
	if(ncol(x)!=nrow(x)){
		warning("\"matrix\" object does not appear to contain a valid Q matrix.\n")
	} else {
		diag(x)<--rowSums(x)
		class(x)<-"Qmatrix"
		return(x)
	}
}

as.Qmatrix.fitMk<-function(x,...){
	Q<-matrix(NA,length(x$states),length(x$states))
	Q[]<-c(0,x$rates)[x$index.matrix+1]
	rownames(Q)<-colnames(Q)<-x$states
	diag(Q)<--rowSums(Q,na.rm=TRUE)
	class(Q)<-"Qmatrix"
	Q
}

as.Qmatrix.ace<-function(x, ...){
	if("index.matrix"%in%names(x)){
		k<-nrow(x$index.matrix)
		Q<-matrix(NA,k,k)
		Q[]<-c(0,x$rates)[x$index.matrix+1]
		rownames(Q)<-colnames(Q)<-colnames(x$lik.anc)
		diag(Q)<--rowSums(Q,na.rm=TRUE)
		class(Q)<-"Qmatrix"
		return(Q)
	} else cat("\"ace\" object does not appear to contain a Q matrix.\n")
}

print.Qmatrix<-function(x,...){
	cat("Estimated Q matrix:\n")
	print(unclass(x),...)
}

is.Qmatrix<-function(x) "Qmatrix" %in% class(x)

