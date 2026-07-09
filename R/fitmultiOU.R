## function to fit a discrete-state-dependent multi-regime Ornstein-Uhlenbeck
## model using the discrete approximation of Boucher & Demery (2016)

fitmultiOU<-function(tree,x,y=NULL,model="ER",ncat=1,...){
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-1
	if(hasArg(maxit)) maxit<-list(...)$maxit
	else maxit<-1000
  if(hasArg(logscale)) logscale<-list(...)$logscale
	else logscale<-TRUE
	if(hasArg(rand_start)) rand_start<-list(...)$rand_start
	else rand_start<-TRUE
	levs<-if(hasArg(levs)) list(...)$levs else 100
	parallel<-if(hasArg(parallel)) list(...)$parallel else 
		FALSE
	lik.func<-if(parallel) "parallel" else "pruning"
	if(hasArg(opt.method)) opt.method<-list(...)$opt.method
	else opt.method<-"nlminb"
	null_model<-if(hasArg(null_model)) list(...)$null_model else
		FALSE
	ncores<-if(hasArg(ncores)) list(...)$ncores else 
		detectCores()-1
	## continuous character
	x<-x[tree$tip.label]
	if(hasArg(lims)) lims<-list(...)$lims
	else lims<-expand.range(x)
	dd<-diff(lims)
	tol<-1e-8*dd/levs
	bins<-cbind(seq(from=lims[1]-tol,by=(dd+2*tol)/levs,
		length.out=levs),seq(to=lims[2]+tol,by=(dd+2*tol)/levs,
			length.out=levs))
	X<-to_binned(x,bins)
#	if(hasArg(wrapped)) wrapped<-list(...)$wrapped
#	else wrapped<-FALSE
	## discrete character
	if(!is.null(y)){
		if(is.matrix(y)){
			y<-y[tree$tip.label,]
			m<-ncol(y)
			states<-colnames(y)
		} else {
			y<-to.matrix(y,sort(unique(y)))
			y<-y[tree$tip.label,]
			m<-ncol(y)
			states<-colnames(y)
		}
	} else {
		y<-matrix(1,nrow(X),1,dimnames=list(rownames(X),"0"))
	}
	if(ncat>1){
		y_tmp<-y
		for(i in 2:ncat){
			colnames(y_tmp)<-paste(colnames(y_tmp),"*",sep="")
			y<-cbind(y,y_tmp)
		}
		y<-y[,order(colnames(y))]
	}
	## combine
	nn<-x_by_y(colnames(y),colnames(X))
	XX<-matrix(0,nrow=nrow(X),ncol=ncol(y)*ncol(X),
		dimnames=list(rownames(X),nn))
	for(i in 1:ncol(y)){
		for(j in 1:nrow(X)){
			XX[j,1:levs+(i-1)*levs]<-X[j,]*y[j,i]
		}
	}
	## set pi
	if(hasArg(root)) root=list(...)$root
	else root<-"mle"
	if(root=="nuisance") pi<-"fitzjohn"
	else if(root=="mle") pi<-"mle"
	else if(root=="flat") pi<-rep(1/ncol(XX),ncol(XX))
	## build continuous model
	cmodel<-matrix(0,nrow=ncol(XX),ncol=ncol(XX),
		dimnames=list(nn,nn))
	for(i in 1:(levs-1)){
		for(j in 0:(ncol(y)-1)){
			cmodel[i+j*levs,i+1+j*levs]<-(i+j*levs-j)
			cmodel[i+1+j*levs,i+j*levs]<-(i+j*levs+ncol(y)*(levs-1)-j)				
		}
	}
	state_ind<-setNames(1:ncol(y),colnames(y))
	if(null_model){ 
		if(ncat==1){ 
		  for(i in 1:(levs-1)){
		    for(j in 0:(ncol(y)-1)){
		      cmodel[i+j*levs,i+1+j*levs]<-(i)
		      cmodel[i+1+j*levs,i+j*levs]<-(i+levs-1)				
		    }
		  }
		}
		else {
		  ## FIGURE THIS OUT LATER
			## allow hidden character to have different rates
			hs<-strsplit(nn,"")
			foo<-function(x){
				a<-paste(x[x!="*"],collapse="")
				b<-paste(x[x=="*"],collapse="")
				return(c(a,b))
			}
			hs<-sapply(hs,foo)[2,]
			un.hs<-sort(unique(hs))
			ind<-which(cmodel>1,arr.ind=TRUE)
			k<-1
			for(i in 1:ncat){
				ii<-which(hs==un.hs[i])
				ii<-ii[!((ii%%levs)==0)]
				state_ind[unique(cmodel[cbind(ii+1,ii)])]<-k
				cmodel[cbind(ii+1,ii)]<-
					cmodel[cbind(ii,ii+1)]<-k
				k<-k+1
			}
		}
	}
	## build discrete model
	dmodel<-matrix(0,nrow=ncol(XX),ncol=ncol(XX),
		dimnames=list(nn,nn))
	qmodel<-matrix(0,nrow=ncol(y),ncol=ncol(y),
		dimnames=list(colnames(y),colnames(y)))
	dn<-strsplit(colnames(qmodel),"")
	foo<-function(x){
		a<-paste(x[x!="*"],collapse="")
		b<-paste(x[x=="*"],collapse="")
		return(c(a,b))
	}
	dn<-sapply(dn,foo)
	q1<-matrix(0,length(unique(dn[1,])),
		length(unique(dn[1,])),
		dimnames=list(sort(unique(dn[1,])),
			sort(unique(dn[1,]))))
	if(is.matrix(model)){
		cust_model<-model
		model<-"custom"
	}
	if(model=="ER"){ 
		q1[]<-1
		diag(q1)<-0
	} else if(model=="SYM"){
		k<-1
		for(i in 1:(nrow(q1)-1)){
			for(j in (i+1):ncol(q1)){
				q1[i,j]<-q1[j,i]<-k
				k<-k+1
			}
		}
	} else if(model=="ARD"){
		k<-1
		for(i in 1:nrow(q1)){
			for(j in 1:ncol(q1)){
				if(i!=j){
					q1[i,j]<-k
					k<-k+1
				}
			}
		}
	} else if(model=="custom"){
		if(is.null(rownames(cust_model))) 
			colnames(cust_model)<-rownames(cust_model)<-rownames(q1)
		q1[]<-cust_model[rownames(q1),colnames(q1)]
		cat("This is the design of your custom discrete-trait model:\n")
		print(q1)
	}
	q2<-matrix(0,length(unique(dn[2,])),
		length(unique(dn[2,])),
		dimnames=list(sort(unique(dn[2,])),
			sort(unique(dn[2,]))))
	if(ncat>1){
		if(hasArg(model.hrm)) model.hrm<-list(...)$model.hrm
		else model.hrm<-if(model%in%c("ER","SYM","ARD")) model else "SYM"
	} else model.hrm<-NULL
	if(!is.null(model.hrm)){
		if(is.matrix(model.hrm)){ 
			q2[]<-model.hrm
			cat("This is the design of your custom hidden-trait model:\n")
			print(q2)
		} else if(model.hrm=="ER"){ 
			q2[]<-max(q1)+1
			diag(q2)<-0
		} else if(model.hrm=="SYM"){
			k<-max(q1)+1
			for(i in 1:(nrow(q2)-1)){
				for(j in (i+1):ncol(q2)){
					q2[i,j]<-q2[j,i]<-k
					k<-k+1
				}
			}
		} else if(model.hrm=="ARD"){
			k<-max(q1)+1
			for(i in 1:nrow(q2)){
				for(j in 1:ncol(q2)){
					if(i!=j){
						q2[i,j]<-k
						k<-k+1
					}
				}
			}
		}
	}
	for(i in 1:nrow(qmodel)){
		for(j in 1:ncol(qmodel)){
			if(dn[1,i]!=dn[1,j]&&dn[2,i]==dn[2,j]) 
				qmodel[i,j]<-q1[dn[1,i],dn[1,j]]
			else if(dn[1,i]==dn[1,j]&&dn[2,i]!=dn[2,j])
				qmodel[i,j]<-q2[which(rownames(q2)==dn[2,i]),
					which(colnames(q2)==dn[2,j])]
		}
	}
	dn<-strsplit(nn,",")
	for(i in 1:ncol(dmodel)){
		for(j in 1:ncol(dmodel)){
			if(dn[[i]][2]==dn[[j]][2]){
				dmodel[i,j]<-qmodel[dn[[i]][1],dn[[j]][1]]
			}
		}
	}
	dmodel[dmodel>0]<-dmodel[dmodel>0]+max(cmodel)
	model<-cmodel+dmodel
	plot_model<-if(hasArg(plot_model)) 
		list(...)$plot_model else FALSE
	## graph model (optional)
	if(plot_model){
		if(hasArg(asp)) asp<-list(...)$asp
		else asp<-1
		if(hasArg(mar)) mar<-list(...)$mar
		else mar<-c(1.1,1.1,4.1,1.1)
		if(hasArg(title)) title<-list(...)$title
		else title<-TRUE
		plot.new()
		dev.hold()
		par(mar=mar)
		if(hasArg(cols)){ 
			cols<-list(...)$cols
			cols<-setNames(cols,0:max(model))
		} else {
			cols<-setNames(c("#f9f9f7",sample(rainbow(n=max(model)))),
				0:max(model))
		}
		plot.window(xlim=c(-0.05*ncol(model),1.1*ncol(model)),
			ylim=c(nrow(model),-0.05*nrow(model)),
			asp=asp)
		for(i in 1:nrow(model)){
			for(j in 1:ncol(model)){
				polygon(x=c(i-1,i,i,i-1),y=c(j-1,j-1,j,j),
					border=FALSE,col=cols[as.character(model[i,j])])
			}
		}
		polygon(c(0,ncol(model),ncol(model),0),
			c(0,0,nrow(model),nrow(model)),border="black")
		for(i in 1:ncol(qmodel)) for(j in 1:nrow(qmodel)){
			polygon(c(0,levs,levs,0)+(i-1)*levs,
				c(0,0,levs,levs)+(j-1)*levs,border="grey")
		}
		for(i in 1:ncol(qmodel)){
			text(x=mean(c(0,levs)+(i-1)*levs),
				y=-0.025*ncol(model),colnames(qmodel)[i])
			text(x=-0.025*ncol(model),
				y=mean(c(0,levs)+(i-1)*levs),
				rownames(qmodel)[i],srt=90)
		}
		if(title) title("structure of discretized model",font.main=3)
		lp<-legend(x=ncol(model)/2,y=1.05*nrow(model),
			legend=names(cols),pch=15,col=cols,bty="n",xpd=TRUE,
			horiz=TRUE,plot=FALSE,cex=0.8)
		for(i in 1:ncol(qmodel)){
			tmp<-paste("[",colnames(qmodel)[i],"]",sep="")
			if(i==1) nm<-bquote(sigma^2~.(tmp))
			else nm<-c(nm,bquote(sigma^2~.(tmp)))
		}
		if(max(qmodel)>=1) nm<-c(nm,paste("q[",1:max(qmodel),"]",sep=""))
		legend(x=ncol(model),y=0,
			legend=nm,pch=15,
			col=cols[c(state_ind+1,1:max(qmodel)+max(cmodel)+1)],
			bty="n",xpd=TRUE,horiz=FALSE,cex=0.8)
		dev.flush()
	}
	xi<-rowMeans(bins)
	delta<-xi[2]-xi[1]
	## if parallel, make cluster
	if(lik.func=="parallel"){
	  if(!exists("ncores")) ncores<-min(nrow(tree$edge),
	    detectCores()-1)
	  mc<-makeCluster(ncores,type="PSOCK")
	  registerDoParallel(cl=mc)
	}
	m<-if(null_model) 1 else ncol(y)
	## likelihood function
  lik.mou<-function(par){
    qq<-vector(mode="numeric",length=max(model))
	  theta<-par[1:m]
	  alpha<-exp(par[m+1])
	  sigsq<-exp(par[m+2])
	  q<-exp(par[(m+3):length(par)])
	  for(i in 1:m){
	    mu<-alpha*(theta[i]-xi)
	    q.up<-sigsq/(2*delta^2)+pmax(mu[1:(levs-1)],0)/delta
	    q.down<-sigsq/(2*delta^2)+pmax(-mu[2:levs],0)/delta
	    qq[1:(levs-1)+(i-1)*(levs-1)]<-q.up
	    qq[1:(levs-1)+(i-1)*(levs-1)+m*(levs-1)]<-q.down
	  }
	  kk<-max(dmodel)-max(cmodel)
	  qq[1:kk+max(cmodel)]<-q
	  if(length(theta)==1){
	    pi<-rep(0,length(xi))
	    pi[which.min((theta-xi)^2)]<-1
	  }
	  if(lik.func=="pruning")
	    lnL<-pruning(qq,tree,XX,model=model,pi=pi)-Ntip(tree)*log(delta)
	  else if(lik.func=="parallel")
	    lnL<-parallel_pruning(qq,tree,XX,model=model,pi=pi)-Ntip(tree)*log(delta)
	  if(trace>0) if(ct%%100==0){
	    cat(paste(c(ct, sprintf("%.4f", c(theta, alpha, sigsq, q, lnL))),
	      collapse = "\t"),
	      "\n")
    }
	  ct<<-ct+1
	  -lnL
  }
  ## initial
  if(hasArg(init)){ 
    init<-list(...)$init
    init<-c(init[1:m],log(init[(m+1):length(init)]))
  } else {
    init<-c(
      runif(n=m,min=min(x),max=max(x)), #theta
      runif(n=1,
        min=log(log(2)/max(nodeHeights(tree))),
        max=log(10*log(2)/max(nodeHeights(tree)))), #alpha
      runif(n=1,
        min=log(0.1*var(x)/max(nodeHeights(tree))),
        max=log(10*var(x)/max(nodeHeights(tree)))), #sigsq
      runif(n=max(qmodel),
        min=log(1/sum(tree$edge.length)),
        max=log(100/sum(tree$edge.length))) #q
    )
  }
  ## fit model
  if(trace>0){
    if(null_model){
      cat(c("iter","theta","alpha","sigsq",
        paste("q[",1:max(qmodel),"]",sep=""),"log(L)\n"),sep="\t")
    } else {
      cat(c("iter",paste("the[",colnames(y),"]",sep=""),
        "alpha","sigsq",paste("q[",1:max(qmodel),"]",sep=""),
        "log(L)\n"),sep="\t")
    }
  }
  ct<-0
	fit<-optim(init,lik.mou,method="Nelder-Mead",
	  control=list(trace=0,maxit=maxit))
	if(trace>0){
	  cat(paste(c(ct,sprintf("%.4f",
	    c(fit$par[1:m],
	      exp(fit$par[m+1]),
	      exp(fit$par[m+2]),
	      exp(fit$par[(m+3):length(fit$par)]),
	     -fit$value))),collapse = "\t"),"\n")
	  cat("Done optimizing.\n\n")
	}
	## stop cluster if parallel
	if(lik.func=="parallel") stopCluster(cl=mc)
	## likelihood function (to export)
  lfunc<-function(theta,alpha,sigsq,q)
    -lik.mou(c(theta,alpha,sigsq,q))
  lnL<--fit$value
  attr(lnL,"df")<-length(fit$par)
  
  qq<-vector(mode="numeric",length=max(model))
  theta<-fit$par[1:m]
  alpha<-exp(fit$par[m+1])
  sigsq<-exp(fit$par[m+2])
  q<-exp(fit$par[(m+3):length(fit$par)])
  for(i in 1:m){
    mu<-alpha*(theta[i]-xi)
    q.up<-sigsq/(2*delta^2)+pmax(mu[1:(levs-1)],0)/delta
    q.down<-sigsq/(2*delta^2)+pmax(-mu[2:levs],0)/delta
    qq[1:(levs-1)+(i-1)*(levs-1)]<-q.up
    qq[1:(levs-1)+(i-1)*(levs-1)+m*(levs-1)]<-q.down
  }
  kk<-max(dmodel)-max(cmodel)
  qq[1:kk+max(cmodel)]<-q
  
  index.matrix<-model
  index.matrix[index.matrix==0]<-NA
  
  mk_fit<-list(
    logLik=-fit$value+Ntip(tree)*log(delta),
    rates=qq,
    index.matrix=index.matrix,
    states=colnames(XX),
    pi=pruning(qq,tree,XX,model=model,pi=pi,return="pi"),
    method="optim",
    root.prior=root,
    opt_results=list(
      convergence=fit$convergence,
      counts=fit$counts,
      message=fit$message),
    data=XX,
    tree=tree,
    lik=function(q) pruning(q,tree,XX,model=model,pi=pi)
  )
  class(mk_fit)<-"fitMk"
  
	object<-list(
	  theta=theta,
	  alpha=alpha,
		sigsq=sigsq,
		bounds=lims,
		rates=q,
		index.matrix=qmodel,
		states=colnames(y),
		ncat=levs,
		logLik=lnL,
		opt_results=list(
		  convergence=fit$convergence,
		  counts=fit$counts,
		  message=fit$message),
	  mk_fit=mk_fit,
		lik=lfunc)
	class(object)<-c("fitmultiOU","fitmultiBM","fitMk")
	object
}

print.fitmultiOU<-function(x,digits=4,...){
	cat(paste("Object of class \"fitmultiOU\" based on\n",
		"   a discretization with k =",
		x$ncat,"levels.\n\n"))
	cat("Fitted multi-theta OU model parameters:\n")
	cat(paste(" levels: [",
		paste(x$states,collapse=", "),
		"]\n"))
	cat(paste("  theta: [",
		paste(round(x$theta,digits),collapse=", "),
		"]\n"))
	cat(paste("  alpha:",round(x$alpha,digits),"\n"))
	cat(paste("  sigsq:",round(x$sigsq,digits),"\n\n"))
	print(as.Qmatrix(x))
	cat(paste("\nLog-likelihood:",round(x$logLik,digits),
		"\n\n"))
	if(x$opt_results$convergence == 0) 
		cat("R thinks it has found the ML solution.\n\n")
	else cat(
		"R thinks optimization may not have converged.\n\n")
}

