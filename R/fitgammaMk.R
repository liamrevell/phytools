## fit model in which the edge rates are distributed according to a 
## discretized gamma distribution with shape parameter alpha

anova.fitgammaMk<-function(object,...) anova.fitMk(object,...)

logLik.fitgammaMk<-function(object,...){
	lik<-object$logLik
	attr(lik,"df")<-length(object$rates)+1
	lik
}

gamma_pruning<-function(par,nrates=4,tree,x,model=NULL,median=TRUE,...){
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
		  P<-Reduce("+",lapply(r,function(rr,k,Q,edge) EXPM(Q*rr*edge)/k,
		    k=nrates,Q=Q,edge=pw$edge.length[ee[j]]))
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
			return(-Inf) else return(prob)
	else if(return=="conditional") L
	else if(return=="pi") pi
}

fitgammaMk<-function(tree,x,model="ER",fixedQ=NULL,nrates=4,...){
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
	else min.alpha<-1e-12
	if(hasArg(max.alpha)) max.alpha<-list(...)$max.alpha
	else max.alpha<-1000
	if(hasArg(logscale)) logscale<-list(...)$logscale
	else logscale<-TRUE
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
			    median=TRUE,pi=pi),method="L-BFGS-B",lower=c(rep(log(min.q),k),
				      log(min.alpha)),upper=c(rep(log(max.q),k),
				        log(max.alpha))) else
				optim(c(q.init,alpha.init),function(p) -pruning(p,nrates=nrates,
				  tree=pw,x=x,model=MODEL,median=TRUE,pi=pi),method="L-BFGS-B",
				  lower=c(rep(min.q,k),min.alpha),upper=c(rep(max.q,k),max.alpha))
		}
		else if(opt.method=="none"){
				fit<-list(objective=-gamma_pruning(c(q.init,alpha.init),
				  nrates=nrates,pw,x,MODEL,median=TRUE,pi=pi),par=q.init)
		} else {
			fit<-if(logscale)
				nlminb(c(q.init,alpha.init),function(p) -gamma_pruning(exp(p),
				  nrates=nrates,tree=pw,x=x,model=MODEL,median=TRUE,pi=pi),
				  lower=c(rep(log(min.q),k),log(min.alpha)),
				  upper=c(rep(log(max.q),k),log(max.alpha))) else
				nlminb(c(q.init,alpha.init),function(p) -gamma_pruning(p,
				  nrates=nrates,tree=pw,x=x,model=MODEL,median=TRUE,pi=pi),
				  lower=c(rep(min.q,k),min.alpha),upper=c(rep(max.q,k),max.alpha))
		}
		if(logscale) fit$par<-exp(fit$par)
		if(pi[1]=="fitzjohn") pi<-setNames(
			gamma_pruning(fit$par,nrates=nrates,tree=tree,x=x,model=MODEL,
			  median=TRUE,pi=pi,return="pi"),states)
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
		  tree=tree,x=x,model=MODEL,median=TRUE,pi=pi,return="conditional")
	} else {
	  fit<-gamma_pruning(c(Q[sapply(1:k,function(x,y) which(x==y),
	    index.matrix)],alpha.init),nrates=nrates,tree=tree,x=x,model=MODEL,
		    median=TRUE,pi=pi)
		if(pi[1]=="fitzjohn") pi<-setNames(gamma_pruning(
		  c(Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],alpha.init),
		  nrates=nrates,tree=tree,x=x,model=MODEL,median=TRUE,pi=pi,
		  return="pi"),states)
		obj<-list(logLik=-fit,
			rates=Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],
			index.matrix=index.matrix,
			states=states,
			pi=pi,
			root.prior=root.prior)
		if(output.liks) obj$lik.anc<-gamma_pruning(
		  c(Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],alpha.init),
		  nrates=nrates,tree=tree,x=x,model=MODEL,median=TRUE,pi=pi,
		  return="conditional")
  }
	lik.f<-function(q,alpha){
	  q<-sapply(1:max(MODEL), function(ind,q,MODEL) q[which(MODEL==ind)],
	    q=q,MODEL=MODEL)
	  gamma_pruning(c(q,alpha),nrates=nrates,tree=pw,x=x,model=MODEL,
	    pi=if(root.prior=="nuisance") "fitzjohn" else pi)
	}
	obj$data<-x
	obj$tree<-tree
	obj$lik<-lik.f
	class(obj)<-"fitgammaMk"
	return(obj)
}

## print method for objects of class "fitgammaMk"
print.fitgammaMk<-function(x,digits=6,...){
  cat("Object of class \"fitmultiMk\".\n\n")
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
