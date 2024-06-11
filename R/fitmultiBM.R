## function to fit a discrete-state-dependent multi-rate Brownian
## model using the discrete approximation of Boucher & Demery (2016)

x_by_y<-function(x,y){
	nn<-vector()
	for(i in 1:length(x)) for(j in 1:length(y))
		nn<-c(nn,paste(x[i],y[j],sep=","))
	nn
}

fitmultiBM<-function(tree,x,y=NULL,model="ER",ncat=1,...){
	levs<-if(hasArg(levs)) list(...)$levs else 100
	parallel<-if(hasArg(parallel)) list(...)$parallel else 
		FALSE
	null_model<-if(hasArg(null_model)) list(...)$null_model else
		FALSE
	ncores<-if(hasArg(ncores)) list(...)$ncores else 
		detectCores()-1
	## continuous character
	x<-x[tree$tip.label]
	lims<-expand.range(x)
	dd<-diff(lims)
	tol<-1e-8*dd/levs
	bins<-cbind(seq(from=lims[1]-tol,by=(dd+2*tol)/levs,
		length.out=levs),seq(to=lims[2]+tol,by=(dd+2*tol)/levs,
			length.out=levs))
	X<-to_binned(x,bins)
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
	## build continuous model
	cmodel<-matrix(0,nrow=ncol(XX),ncol=ncol(XX),
		dimnames=list(nn,nn))
	for(i in 1:(levs-1)){
		for(j in 0:(ncol(y)-1)){
			cmodel[i+j*levs,i+1+j*levs]<-
				cmodel[i+1+j*levs,i+j*levs]<-(j+1)
		}
	}
	state_ind<-setNames(1:ncol(y),colnames(y))
	if(null_model){ 
		if(ncat==1){ 
			cmodel[cmodel>0]<-1
			state_ind[]<-1
		}
		else {
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
	}
	q2<-matrix(0,length(unique(dn[2,])),
		length(unique(dn[2,])),
		dimnames=list(sort(unique(dn[2,])),
			sort(unique(dn[2,]))))
	if(model=="ER"){ 
		q2[]<-max(q1)+1
		diag(q2)<-0
	} else if(model=="SYM"){
		k<-max(q1)+1
		for(i in 1:(nrow(q2)-1)){
			for(j in (i+1):ncol(q2)){
				q2[i,j]<-q2[j,i]<-k
				k<-k+1
			}
		}
	} else if(model=="ARD"){
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
		plot.new()
		dev.hold()
		par(mar=c(5.1,2.1,4.1,2.1))
		cols<-setNames(c("#f9f9f7",sample(rainbow(n=max(model)))),
			0:max(model))
		plot.window(xlim=c(0,ncol(model)),ylim=c(nrow(model),0))
		for(i in 1:nrow(model)){
			for(j in 1:ncol(model)){
				polygon(x=c(i-1,i,i,i-1),y=c(j-1,j-1,j,j),
					border=FALSE,col=cols[as.character(model[i,j])])
			}
		}
		polygon(c(0,ncol(model),ncol(model),0),
			c(0,0,nrow(model),nrow(model)),border="grey")
		## box(col="grey")
		title("structure of discretized model",font.main=3)
		lp<-legend(x=ncol(model)/2,y=1.05*nrow(model),
			legend=names(cols),pch=15,col=cols,bty="n",xpd=TRUE,
			horiz=TRUE,plot=FALSE,cex=0.8)
		legend(x=ncol(model)/2-lp$rect$w/2,y=1.05*nrow(model),
			legend=names(cols),pch=15,col=cols,bty="n",xpd=TRUE,
			horiz=TRUE,cex=0.8)
		dev.flush()
	}
	## set initial values for optimization
	q.init<-runif(n=max(cmodel)+max(qmodel),0,2)*c(
		rep((1/2)*mean(pic(x,multi2di(tree))^2)*(levs/dd)^2,
			max(cmodel)),
		rep(fitMk(tree,y,model="ER")$rates,max(qmodel)))
	## optimize model
	fit<-fitMk(tree,XX,model=model,
		lik.func=if(parallel) "parallel" else "pruning",
		expm.method=if(isSymmetric(model)) "R_Eigen" else 
			"Higham08.b",
		pi="mle",logscale=TRUE,q.init=q.init)
	sig2<-2*fit$rates[1:max(cmodel)]*(dd/levs)^2
	q<-fit$rates[1:max(qmodel)+max(cmodel)]
	index.matrix<-qmodel
	lnL<-logLik(fit)-Ntip(tree)*log(dd/levs)
	attr(lnL,"df")<-max(model)+1
	## likelihood function (to export)
	lfunc<-function(sig2,q,x0="mle",...){
		if(hasArg(lik.func)) lik.func<-list(...)$lik.func
		else lik.func<-"pruning"
		if(hasArg(parallel)) parallel<-list(...)$parallel
		else parallel<-FALSE
		q<-c((sig2/2)*(levs/dd)^2,q)
		if(x0=="mle") pi<-"mle"
		else if(x0=="nuisance") pi<-"fitzjohn"
		else if(is.numeric(x0)) pi<-to_binned(x0,bins)[1,]
		if(lik.func=="pruning"){
			lnL<-pruning(q,tree,XX,model,pi=pi,
				expm.method=if(isSymmetric(model)) "R_Eigen" else 
					"Higham08.b")-
				Ntip(tree)*log(dd/levs)
		} else if(lik.func=="parallel"){
			if(!exists("ncores")) ncores<-min(nrow(tree$edge),
				detectCores()-1)
			mc<-makeCluster(ncores,type="PSOCK")
			registerDoParallel(cl=mc)
			lnL<-parallel_pruning(q,tree,XX,model,pi=pi,
				expm.method=if(isSymmetric(model)) "R_Eigen" else 
					"Higham08.b")-Ntip(tree)*log(dd/levs)
			stopCluster(cl=mc)
		}
		lnL
	}
	object<-list(
		sigsq=sig2,
		state_ind=state_ind,
		x0=sum(fit$pi*rowMeans(bins)),
		rates=q,
		index.matrix=qmodel,
		states=colnames(y),
		ncat=levs,
		logLik=lnL,
		opt_results=fit$opt_results,
		mk_fit=fit,
		lik=lfunc)
	class(object)<-c("fitmultiBM","fitMk")
	object
}

print.fitmultiBM<-function(x,digits=4,...){
	cat(paste("Object of class \"fitmultiBM\" based on\n",
		"   a discretization with k =",
		x$ncat,"levels.\n\n"))
	cat("Fitted multi-rate BM model parameters:\n")
	cat(paste(" levels: [",
		paste(names(x$state_ind),collapse=", "),
		"]\n"))
	cat(paste("  sigsq: [",
		paste(round(x$sigsq[x$state_ind],digits),collapse=", "),
		"]\n"))
	cat(paste("     x0:",round(x$x0,digits),"\n\n"))
	print(as.Qmatrix(x))
	cat(paste("\nLog-likelihood:",round(x$logLik,digits),
		"\n\n"))
	if(x$opt_results$convergence == 0) 
		cat("R thinks it has found the ML solution.\n\n")
	else cat(
		"R thinks optimization may not have converged.\n\n")
}

logLik.fitmultiBM<-function(object,...) object$logLik