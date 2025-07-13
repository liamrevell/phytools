## function to fit a discrete-state-dependent multi-rate, multi-trend Brownian
## model using the discrete approximation of Boucher & Demery (2016)

fitmultiTrend<-function(tree,x,y=NULL,model="ER",ncat=1,...){
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
	delta<-dd/levs
	tol<-1e-8*delta
	bins<-cbind(seq(from=lims[1]-tol,by=(dd+2*tol)/levs,
		length.out=levs),seq(to=lims[2]+tol,by=(dd+2*tol)/levs,
			length.out=levs))
	X<-to_binned(x,bins)
	if(hasArg(wrapped)) wrapped<-list(...)$wrapped
	else wrapped<-FALSE
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
			cmodel[i+j*levs,i+1+j*levs]<-(j+1)
			cmodel[i+1+j*levs,i+j*levs]<-(j+1)+ncol(y)				
		}
	}
	if(wrapped){
		for(i in 0:(ncol(y)-1)){
			cmodel[1+i*levs,(i+1)*levs]<-(i+1)+ncol(y)
			cmodel[(i+1)*levs,1+i*levs]<-(i+1)
		}
	}
	state_ind<-setNames(1:ncol(y),colnames(y))
	if(null_model){ 
		if(ncat==1){ 
			for(i in 2:ncol(cmodel)){
				if(cmodel[i-1,i]>0) cmodel[i-1,i]<-1
				if(cmodel[i,i-1]>0) cmodel[i,i-1]<-2
			}
			state_ind[]<-1
		} else {
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
				cmodel[cbind(ii+1,ii)]<-k
				cmodel[cbind(ii,ii+1)]<-k+ncat
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
			tmp<-paste("[",colnames(qmodel)[i],"] ",sep="")
			if(i==1) nm<-bquote(sigma^2~.(tmp)~-mu~.(tmp))
			else nm<-c(nm,bquote(sigma^2~.(tmp)~-mu~.(tmp)))
		}
		for(i in 1:nrow(qmodel)){
			tmp<-paste("[",colnames(qmodel)[i],"] ",sep="")
			nm<-c(nm,bquote(sigma^2~.(tmp)~+mu~.(tmp)))
		}
		if(max(qmodel)>=1) nm<-c(nm,paste("q[",1:max(qmodel),"]",sep=""))
		legend(x=ncol(model),y=0,
			legend=nm,pch=15,
			col=cols[2:length(cols)],
			bty="n",xpd=TRUE,horiz=FALSE,cex=0.8)
		dev.flush()
	}
	## set initial values for optimization
	qq<-if(ncol(y)>1) fitMk(tree,y,model="ER")$rates else 
		1/sum(tree$edge.length)
	q.init<-c(rep((1/2)*mean(pic(x,multi2di(tree))^2)*(levs/dd)^2,
		max(cmodel)),rep(qq,max(qmodel)))
	if(rand_start) q.init<-q.init*runif(n=length(q.init),0,2)
	max.q<-max(q.init)*1000
	## optimize model
	fit<-fitMk(tree,XX,model=model,
		lik.func=if(parallel) "parallel" else "pruning",
		expm.method=if(isSymmetric(model)) "R_Eigen" else 
			"Higham08.b",
		pi=pi,logscale=logscale,q.init=q.init,
		opt.method=opt.method,max.q=max.q,
		ncores=ncores)
	sig2<-(fit$rates[1:(max(cmodel)/2)]+
		fit$rates[1:(max(cmodel)/2)+max(cmodel)/2])*delta^2
	mu<-(fit$rates[1:(max(cmodel)/2)]-
		fit$rates[1:(max(cmodel)/2)+max(cmodel)/2])*delta
	q<-fit$rates[1:max(qmodel)+max(cmodel)]
	index.matrix<-qmodel
	lnL<-logLik(fit)-Ntip(tree)*log(delta)
	attr(lnL,"df")<-max(model)+1
	## likelihood function (to export)
	lfunc<-function(sig2,mu,q,x0="mle",...){
		if(hasArg(lik.func)) lik.func<-list(...)$lik.func
		else lik.func<-"pruning"
		if(hasArg(parallel)) parallel<-list(...)$parallel
		else parallel<-FALSE
		delta<-dd/levs
		jj<-1:(length(sig2)/2)
		kk<-1:(length(sig2)/2)+length(sig2)/2
		q<-c(sig2[jj]/(2*delta^2)+mu[jj]/(2*delta),
			sig2[kk]/(2*delta^2)-mu[kk]/(2*delta),
			q)
		if(x0=="mle") pi<-"mle"
		else if(x0=="nuisance") pi<-"fitzjohn"
		else if(is.numeric(x0)) pi<-to_binned(x0,bins)[1,]
		if(lik.func=="pruning"){
			lnL<-pruning(q,tree,XX,model,pi=pi,
				expm.method=if(isSymmetric(model)) "R_Eigen" else 
					"Higham08.b")-
				Ntip(tree)*log(delta)
		} else if(lik.func=="parallel"){
			if(!exists("ncores")) ncores<-min(nrow(tree$edge),
				detectCores()-1)
			mc<-makeCluster(ncores,type="PSOCK")
			registerDoParallel(cl=mc)
			lnL<-parallel_pruning(q,tree,XX,model,pi=pi,
				expm.method=if(isSymmetric(model)) "R_Eigen" else 
					"Higham08.b")-Ntip(tree)*log(delta)
			stopCluster(cl=mc)
		}
		lnL
	}
	object<-list(
		sigsq=sig2,
		mu=mu,
		state_ind=state_ind,
		x0=sum(fit$pi*rowMeans(bins)),
		bounds=lims,
		rates=q,
		index.matrix=qmodel,
		states=colnames(y),
		ncat=levs,
		logLik=lnL,
		opt_results=fit$opt_results,
		mk_fit=fit,
		lik=lfunc)
	class(object)<-c("fitmultiTrend","fitmultiBM","fitMk")
	object
}

print.fitmultiTrend<-function(x,digits=4,...){
	cat(paste("Object of class \"fitmultiTrend\" based on\n",
		"   a discretization with k =",
		x$ncat,"levels.\n\n"))
	cat("Fitted multi-rate trend model parameters:\n")
	cat(paste(" levels: [",
		paste(names(x$state_ind),collapse=", "),
		"]\n"))
	cat(paste("  sigsq: [",
		paste(round(x$sigsq[x$state_ind],digits),collapse=", "),
		"]\n"))
	cat(paste("  mu: [",
		paste(round(x$mu[x$state_ind],digits),collapse=", "),
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

sim.multiTrend<-function(tree,x0=0,sig2=1,mu=0,...){
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	if(inherits(tree,"simmap")){ 
		tt<-map.to.singleton(tree)
	} else { 
		tt<-tree
		tt$edge.length<-setNames(tt$edge.length,
			rep("1",length(tt$edge.length)))
	}
	ss<-sort(unique(names(tt$edge.length)))
	if(length(sig2)==1) sig2<-rep(sig2,length(ss))
	if(length(mu)==1) mu<-rep(mu,length(ss))
	if(is.null(names(sig2))) names(sig2)<-ss
	if(is.null(names(mu))) names(mu)<-ss
	sig2<-sig2[names(tt$edge.length)]
	mu<-mu[names(tt$edge.length)]
	delta<-rnorm(n=nrow(tt$edge),
		mean=tt$edge.length*mu,sd=sqrt(tt$edge.length*sig2))
	tt<-reorder(tt)
	X<-matrix(x0,nrow(tt$edge),ncol(tt$edge))
	for(i in 1:nrow(tt$edge)){
		X[i,2]<-X[i,1]+delta[i]
		ii<-which(tt$edge[,1]==tt$edge[i,2])
		if(length(ii)>0) X[ii,1]<-X[i,2]
	}
	if(plot){
		xx<-c(
			sapply(1:Ntip(tt),
			function(ii,ee,x) x[which(ee[,2]==ii),2],
			ee=tt$edge,x=X),
			X[1,1],
			sapply(2:tt$Nnode+Ntip(tt),
				function(ii,ee,x) x[which(ee[,2]==ii),2],
				ee=tt$edge,x=X))
		xx<-setNames(xx,c(tt$tip.label,
			1:tt$Nnode+Ntip(tt)))
		if(hasArg(col)) col<-list(...)$col
		else col<-palette()[1:length(ss)]
		edge_col<-if(is.null(names(col))) setNames(col,ss) else col
		tt$maps<-as.list(tt$edge.length)
		for(i in 1:length(tt$maps)) 
			tt$maps[[i]]<-setNames(tt$maps[[i]],
				names(tt$edge.length)[i])
		phenogram(tt,xx,ftype="off",colors=edge_col,lwd=1)
	}
	setNames(
		sapply(1:Ntip(tt),
			function(ii,ee,x) x[which(ee[,2]==ii),2],
			ee=tt$edge,x=X),
		tt$tip.label)
}
