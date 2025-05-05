## fit model in which the edge rates are distributed according to a 
## discretized gamma distribution with shape parameter alpha

plot.fitgammaMk<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-3
	if(hasArg(colors)) colors<-list(...)$colors
	else colors<-c("yellow","red")
	if(hasArg(title)) title<-list(...)$title
	else title<-"relative edge rate"
	if(is.null(x$marginal)){
		stop("missing marginal likelihoods.")
	} else {
		r<-qgamma(seq(1/(2*x$nrates),1,by=1/x$nrates),x$alpha,x$alpha)
		r<-r/mean(r)
		Rates<-log(apply(x$marginal,1,function(x,y) sum(x*y),y=r))
		cols<-setNames(colorRampPalette(colors)(101),
			0:100)
		args<-list(...)
		args$title<-NULL
		if(is.null(args$type)) args$type<-"phylogram"
		if(is.null(args$direction)) args$direction<-"rightwards"
		if(is.null(args$fsize)){
			if(args$type%in%c("phylogram","cladogram")){
				if(args$direction%in%c("rightwards","leftwards"))
					args$fsize<-min(c(6*par()$pin[2]/Ntip(x$tree),1))
				else
					args$fsize<-min(c(6*par()$pin[1]/Ntip(x$tree),1))
			} else {
				args$fsize<-min(c(0.6*min(par()$pin)/sqrt(Ntip(x$tree)),1))
			}
		}
		args$plot<-FALSE
		args$tree<-x$tree
		nulo<-do.call(plotTree,args)
		pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		ss<-round((Rates-min(Rates))/diff(range(Rates))*100)
		tt<-paintBranches(x$tree,edge=x$tree$edge[1,2],state=ss[1])
		for(j in 2:length(ss)) tt<-paintBranches(tt,edge=x$tree$edge[j,2],
			state=ss[j])
		args$plot<-TRUE
		args$colors<-cols
		args$xlim<-c(-0.3*pp$x.lim[2],pp$x.lim[2])
		args$ylim<-pp$y.lim
		args$add<-TRUE
		args$split.vertical<-TRUE
		args$tree<-tt
		nulo<-do.call(plotSimmap,args)
		pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		h<-max(nodeHeights(x$tree))
		LWD<-diff(par()$usr[1:2])/dev.size("px")[1]
		Nt<-Ntip(x$tree)
		lines(x=rep(-0.25*h+LWD*15/2,2),y=c(1+1/40*Nt,Nt-1/40*Nt))
		nticks<-10
		Y<-cbind(seq(1+1/40*Nt,Nt-1/40*Nt,length.out=nticks),
			seq(1+1/40*Nt,Nt-1/40*Nt,length.out=nticks))
		X<-cbind(rep(-0.25*h+LWD*15/2,nticks),
			rep(-0.25*h+LWD*15/2+0.02*h,nticks))
		for(i in 1:nrow(Y)) lines(X[i,],Y[i,])
		add.color.bar(Nt-2/40*Nt-1,cols,
			title=title,
			lims=NULL,digits=3,
			direction="upwards",
			subtitle="",lwd=15,
			x=-0.25*h,
			y=1+1/40*Nt,prompt=FALSE)
		ticks<-exp(seq(min(Rates),max(Rates),length.out=10))
		text(x=X[,2],y=Y[,2],signif(ticks,digits),pos=4,cex=0.7)
		invisible(exp(Rates))
	}
}

as.Qmatrix.fitgammaMk<-function(x,...) as.Qmatrix.fitMk(x,...)

anova.fitgammaMk<-function(object,...) anova.fitMk(object,...)

logLik.fitgammaMk<-function(object,...){
	lik<-object$logLik
	attr(lik,"df")<-length(object$rates)+1
	lik
}

gamma_pruning<-function(par,nrates=4,tree,x,model=NULL,median=TRUE,
	expm.method="Higham08.b",...){
	if(hasArg(fn_min)) fn_min<-list(...)$fn_min
	else fn_min<--Inf
	if(hasArg(marginal)) marginal<-list(...)$marginal
	else marginal<-FALSE
	if(marginal){
		if(hasArg(edge)) edge<-list(...)$edge
		else marginal<-FALSE
		if(hasArg(rate)) rate<-list(...)$rate
		else marginal<-FALSE
	}
	q<-par[1:(length(par)-1)]
	alpha<-par[length(par)]
	if(median){
		r<-qgamma(seq(1/(2*nrates),1,by=1/nrates),alpha,alpha)
		r<-r/mean(r)
	} else {
		cat("This does not work yet.\n")
	}
	if(hasArg(return)) return<-list(...)$return
	else return<-"likelihood"
	pw<-if(!is.null(attr(tree,"order"))&&
		attr(tree,"order")=="postorder") tree else 
		reorder(tree,"postorder")
	k<-ncol(x)
	if(is.null(model)){
		model<-matrix(1,k,k)
		diag(model)<-0
	}
	if(hasArg(pi)) pi<-list(...)$pi
	else pi<-rep(1/k,k)
	Q<-matrix(0,k,k)
	Q[]<-c(0,q)[model+1]
	diag(Q)<--rowSums(Q)
	L<-rbind(x[pw$tip.label,],
		matrix(0,pw$Nnode,k,
		dimnames=list(1:pw$Nnode+Ntip(pw))))
	nn<-unique(pw$edge[,1])
	pp<-vector(mode="numeric",length=length(nn))
	root<-min(nn)
	for(i in 1:length(nn)){
		ee<-which(pw$edge[,1]==nn[i])
		PP<-matrix(NA,length(ee),k)
		for(j in 1:length(ee)){
			if(marginal){
				if(pw$edge[ee[j],2]==edge){ 
					ind<-rate
				} else ind<-1:nrates
			} else ind<-1:nrates
			P<-Reduce("+",lapply(r[ind],
				function(rr,nr,Q,edge) expm(Q*rr*edge,method=expm.method)/nr,
				nr=nrates,Q=Q,edge=pw$edge.length[ee[j]]))
			PP[j,]<-P%*%L[pw$edge[ee[j],2],]
		}
		L[nn[i],]<-apply(PP,2,prod)
		if(nn[i]==root){
			if(pi[1]=="fitzjohn") pi<-L[nn[i],]/sum(L[nn[i],])
			L[nn[i],]<-pi*L[nn[i],]
		}
		pp[i]<-sum(L[nn[i],])
		L[nn[i],]<-L[nn[i],]/pp[i]
	}
	prob<-sum(log(pp))
	if(return=="likelihood") 
		if(is.na(prob)||is.nan(prob)) 
			return(fn_min) else return(prob)
	else if(return=="conditional") L
	else if(return=="pi") pi
}

fitgammaMk<-function(tree,x,model="ER",fixedQ=NULL,nrates=8,...){
	median<-TRUE
	if(hasArg(fn_min)) fn_min<-list(...)$fn_min
	else fn_min<--Inf
	if(hasArg(marginal)) marginal<-list(...)$marginal
	else marginal<-FALSE
	if(hasArg(parallel)) parallel<-list(...)$parallel
	else parallel<-TRUE
	if(hasArg(opt.method)) opt.method<-list(...)$opt.method
	else opt.method<-"nlminb"
	if(hasArg(output.liks)) output.liks<-list(...)$output.liks
	else output.liks<-FALSE
	if(hasArg(smart_start)) smart_start<-list(...)$smart_start
	else smart_start<-FALSE
	if(hasArg(q.init)) q.init<-list(...)$q.init
	else q.init<-length(unique(x))/sum(tree$edge.length)
	if(hasArg(alpha.init)) alpha.init<-list(...)$alpha.init
	else alpha.init<-1.0
	if(hasArg(rand_start)) rand_start<-list(...)$rand_start
	else rand_start<-FALSE
	if(hasArg(min.q)) min.q<-list(...)$min.q
	else min.q<-1e-12
	if(hasArg(max.q)) max.q<-list(...)$max.q
	else max.q<-max(nodeHeights(tree))*100
	if(hasArg(min.alpha)) min.alpha<-list(...)$min.alpha
	else min.alpha<-0.1
	if(hasArg(max.alpha)) max.alpha<-list(...)$max.alpha
	else max.alpha<-1000
	if(hasArg(logscale)) logscale<-list(...)$logscale
	else logscale<-TRUE
	if(hasArg(expm.method)) expm.method<-list(...)$expm.method
	else expm.method<-"Higham08.b"
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
	} else if(pi[1]=="fitzjohn") root.prior<-"nuisance"
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
	MODEL<-rate
	MODEL[is.na(MODEL)]<-0
	diag(MODEL)<-0
	tmp<-cbind(1:m,1:m)
	rate[tmp]<-0
	rate[rate==0]<-k+1
	liks<-rbind(x,matrix(0,M,m,dimnames=list(1:M+N,states)))
	pw<-reorder(tree,"postorder")
	if(is.null(fixedQ)){
		if(smart_start&&max(index.matrix,na.rm=TRUE)>1){
			MM<-index.matrix
			MM[is.na(MM)]<-0
			MM[MM>0]<-1
			q.init<-fitMk(pw,x,model=MM,pi=pi,opt.method="nlminb")$rates
		}
		if(length(q.init)!=k) q.init<-rep(q.init[1],k)
		if(rand_start){
			q.init<-q.init*rexp(length(q.init),1)
			alpha.init<-alpha.init*rexp(1)
		}
		q.init<-if(logscale) log(q.init) else q.init
		alpha.init<-if(logscale) log(alpha.init) else alpha.init
		if(opt.method=="optim"){
			fit<-if(logscale)
				optim(c(q.init,alpha.init),function(p) 
					-gamma_pruning(exp(p),nrates=nrates,tree=pw,x=x,model=MODEL,
					median=TRUE,pi=pi,fn_min=fn_min,expm.method=expm.method),
					method="L-BFGS-B",lower=c(rep(log(min.q),k),log(min.alpha)),
					upper=c(rep(log(max.q),k),log(max.alpha))) else
				optim(c(q.init,alpha.init),function(p) -gamma_pruning(p,nrates=nrates,
					tree=pw,x=x,model=MODEL,median=TRUE,pi=pi,fn_min=fn_min,
					expm.method=expm.method),method="L-BFGS-B",
					lower=c(rep(min.q,k),min.alpha),upper=c(rep(max.q,k),max.alpha))
		}
		else if(opt.method=="none"){
			fit<-list(objective=-gamma_pruning(c(q.init,alpha.init),
				nrates=nrates,pw,x,MODEL,median=TRUE,pi=pi),par=q.init,
				expm.method=expm.method)
		} else {
			fit<-if(logscale)
				nlminb(c(q.init,alpha.init),function(p) -gamma_pruning(exp(p),
					nrates=nrates,tree=pw,x=x,model=MODEL,median=TRUE,pi=pi,
					fn_min=fn_min,expm.method=expm.method),
					lower=c(rep(log(min.q),k),log(min.alpha)),
					upper=c(rep(log(max.q),k),log(max.alpha))) else
				nlminb(c(q.init,alpha.init),function(p) -gamma_pruning(p,
					nrates=nrates,tree=pw,x=x,model=MODEL,median=TRUE,pi=pi,
					fn_min=fn_min,expm.method=expm.method),
					lower=c(rep(min.q,k),min.alpha),
					upper=c(rep(max.q,k),max.alpha))
		}
		if(logscale) fit$par<-exp(fit$par)
		if(pi[1]=="fitzjohn") pi<-setNames(
			gamma_pruning(fit$par,nrates=nrates,tree=tree,x=x,model=MODEL,
				median=TRUE,pi=pi,return="pi",expm.method=expm.method),states)
		obj<-list(logLik=
			if(opt.method=="optim") -fit$value else -fit$objective,
			rates=fit$par[1:(length(fit$par)-1)],
			index.matrix=index.matrix,
			states=states,
			pi=pi,
			method=opt.method,
			root.prior=root.prior,
			nrates=nrates,
			alpha=fit$par[length(fit$par)])
		if(opt.method=="nlminb")
			obj$opt_results<-fit[c("convergence","iterations","evaluations","message")]
		else if(opt.method=="optim")
			obj$opt_results<-fit[c("counts","convergence","message")]
		if(output.liks) obj$lik.anc<-gamma_pruning(fit$par,nrates=nrates,
			tree=tree,x=x,model=MODEL,median=TRUE,pi=pi,return="conditional",
			expm.method=expm.method)
	} else {
		fit<-gamma_pruning(c(Q[sapply(1:k,function(x,y) which(x==y),
			index.matrix)],alpha.init),nrates=nrates,tree=tree,x=x,model=MODEL,
				median=TRUE,pi=pi,expm.method=expm.method)
		if(pi[1]=="fitzjohn") pi<-setNames(gamma_pruning(
			c(Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],alpha.init),
			nrates=nrates,tree=tree,x=x,model=MODEL,median=TRUE,pi=pi,
			return="pi",expm.method=expm.method),states)
		obj<-list(logLik=fit,
			rates=Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],
			index.matrix=index.matrix,
			states=states,
			pi=pi,
			root.prior=root.prior,
			nrates=nrates,
			alpha=alpha.init)
		if(output.liks) obj$lik.anc<-gamma_pruning(
			c(Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],alpha.init),
			nrates=nrates,tree=tree,x=x,model=MODEL,median=TRUE,pi=pi,
			return="conditional",expm.method=expm.method)
	}
	if(marginal){
		## get marginal likelihoods of each rate on each edge
		cat(paste("  --\n  Computing marginal scaled likelihoods",
		if(parallel) "(in parallel)" else "(in serial)",
		"of each\n"))
		cat("  rate on each edge. Caution: this is NOT fast....\n  --\n")
		flush.console()
		if(median){
			Rates<-qgamma(seq(1/(2*nrates),1,by=1/nrates),obj$alpha,obj$alpha)
			Rates<-Rates/mean(Rates)
		} else Rates<-1:nrates
		if(parallel){
			ncores<-min(c(parallel::detectCores()-2,nrow(tree$edge)))
			mc<-makeCluster(ncores,type="PSOCK")
			registerDoParallel(cl=mc)
			tmpRATES<-foreach(i=1:nrow(tree$edge))%dopar%{
				foo<-function(X) phytools::gamma_pruning(
					par=c(obj$rates,obj$alpha),
					nrates=nrates,tree=tree,x=x,model=MODEL,median=TRUE,
					pi=pi,marginal=TRUE,edge=tree$edge[i,2],rate=X,
					expm.method=expm.method)
				sapply(1:nrates,foo)
			}
			stopCluster(cl=mc)
			RATES<-t(sapply(tmpRATES,function(x) x))
			dimnames(RATES)<-list(apply(tree$edge,1,
				function(x) paste(x[1],",",x[2],sep="")),
				round(Rates,6))
		} else {
			RATES<-matrix(NA,nrow(tree$edge),nrates,
				dimnames=list(apply(tree$edge,1,
				function(x) paste(x[1],",",x[2],sep="")),
				round(Rates,6)))
			for(i in 1:nrow(RATES)){
				for(j in 1:ncol(RATES)){
					RATES[i,j]<-gamma_pruning(c(obj$rates,obj$alpha),nrates=nrates,
						tree=tree,x=x,model=MODEL,median=TRUE,pi=pi,marginal=TRUE,
						edge=tree$edge[i,2],rate=j,expm.method=expm.method)
				}
			}
		}
		RATES<-t(apply(RATES,1,function(x) exp(x)/sum(exp(x))))
	}
	lik.f<-function(q,alpha){
		q<-sapply(1:max(MODEL), function(ind,q,MODEL) q[which(MODEL==ind)],
			q=q,MODEL=MODEL)
		gamma_pruning(c(q,alpha),nrates=nrates,tree=pw,x=x,model=MODEL,
			pi=if(root.prior=="nuisance") "fitzjohn" else pi,
			expm.method=expm.method)
	}
	obj$data<-x
	obj$tree<-tree
	if(marginal) obj$marginal<-RATES
	obj$lik<-lik.f
	class(obj)<-"fitgammaMk"
	return(obj)
}

## print method for objects of class "fitgammaMk"
print.fitgammaMk<-function(x,digits=6,...){
	cat("Object of class \"fitgammaMk\".\n\n")
	cat("Fitted (or set) value of Q:\n")
	Q<-matrix(NA,length(x$states),length(x$states))
	Q[]<-c(0,x$rates)[x$index.matrix+1]
	diag(Q)<-0
	diag(Q)<--rowSums(Q)
	colnames(Q)<-rownames(Q)<-x$states
	print(round(Q,digits))
	cat(paste("\nFitted (or set) value of alpha rate heterogeneity\n(with",
		x$nrates,"rate categories):",round(x$alpha,digits)))
	cat("\n\nFitted (or set) value of pi:\n")
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
