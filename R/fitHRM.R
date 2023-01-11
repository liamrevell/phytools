## this function fits a hidden-rates model (Beaulieu et al. 2013)
## written by Liam J. Revell 2020, 2021, 2023

fitHRM<-function(tree,x,model="ARD",ncat=2,...){
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-0
	if(!is.factor(x)) x<-setNames(as.factor(x),names(x))
	k<-length(levels(x))
	if(hasArg(niter)) niter<-list(...)$niter
	else niter<-10
	if(length(ncat)==1) ncat<-rep(ncat,k)
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(hasArg(parallel)) parallel<-list(...)$parallel
	else parallel<-FALSE
	if(is.numeric(model)||(model%in%c("ARD","SYM","ER")==FALSE)){
		cat("Only models \"ER\", \"SYM\", and \"ARD\" are permitted.\n")
		cat("Setting model to \"ARD\".\n\n")
		model<-"ARD"
	}
	if(hasArg(umbral)) umbral<-list(...)$umbral
	else umbral<-FALSE
	if(hasArg(corHMM_model)) corHMM_model<-list(...)$corHMM_model
	else corHMM_model<-FALSE
	if(umbral&&corHMM_model){
		cat("The \"umbral\" and corHMM models are inconsistent. Setting \"corHMM_model=FALSE\".\n")
		corHMM_model<-FALSE
	}
	XX<-to.matrix(x,levels(x))
	X<-matrix(NA,nrow(XX),sum(ncat),dimnames=list(rownames(XX)))
	ii<-1
	cols<-nn<-vector()
	for(i in 1:ncol(XX)){
		for(j in 1:ncat[i]){ 
			X[,ii]<-XX[,i]
			nn[ii]<-paste(colnames(XX)[i],"_R",j,sep="")
			cols[ii]<-paste(colnames(XX)[i],paste(rep("*",j-1),
				collapse=""),sep="")
			ii<-ii+1
		}
	}
	colnames(X)<-nn
	MODEL<-model
	model<-matrix(0,ncol(X),ncol(X),dimnames=list(colnames(X),
		colnames(X)))
	if(!umbral){
		ii<-1
		for(i in 1:length(nn)){
			FROM<-strsplit(nn[i],"_")[[1]]
			TO<-strsplit(nn,"_")
			for(j in 1:length(nn)){
				if(i!=j){
					if(FROM[1]==TO[[j]][1]){
						r<-as.numeric(
							c(strsplit(FROM[2],"R")[[1]][2],
							strsplit(TO[[j]][2],"R")[[1]][2]))
						if(abs(diff(r))==1){
							model[i,j]<-ii
							ii<-ii+1
						}
					} else if(FROM[2]==TO[[j]][2]){
						model[i,j]<-ii
						ii<-ii+1
					}
				}
			}
		}
	} else {
		if(hasArg(ordered)) ordered<-list(...)$ordered
		else ordered<-TRUE
		order<-if(ordered) setNames(1:length(levels(x)),levels(x)) else 
			setNames(rep(1,length(levels(x))),levels(x))
		ii<-1
		for(i in 1:length(nn)){
			FROM<-strsplit(nn[i],"_")[[1]]
			TO<-strsplit(nn,"_")
			for(j in 1:length(nn)){
				if(i!=j){
					if(FROM[1]==TO[[j]][1]){
						r<-as.numeric(
							c(strsplit(FROM[2],"R")[[1]][2],
							strsplit(TO[[j]][2],"R")[[1]][2]))
						if(abs(diff(r))==1){
							model[i,j]<-ii
							ii<-ii+1
						}
					} else {
						r<-order[c(FROM[1],TO[[j]][1])]
						if(abs(diff(r))<=1){
							if(as.numeric(strsplit(FROM[2],"R")[[1]][2])==1&&
								as.numeric(strsplit(TO[[j]][2],"R")[[1]][2])==1){
								model[i,j]<-ii
								ii<-ii+1
							}
						
						}
					}
				}
			}
		}
	}
	if(MODEL=="SYM"){
		for(i in 1:nrow(model)) for(j in i:ncol(model)) model[i,j]<-model[j,i]
		oo<-sort(unique(as.vector(model)))
		ii<-0:length(oo)
		for(i in 1:length(oo)) model[model==oo[i]]<-ii[i]
	} else if(MODEL=="ER"){
		for(i in 1:nrow(model)) for(j in 1:ncol(model)){
			if(model[i,j]!=0){
				nr<-rownames(model)[i]
				nc<-rownames(model)[j]
				if(strsplit(nr,"_")[[1]][2]==strsplit(nc,"_")[[1]][2]){
					model[i,j]<-as.numeric(strsplit(nr,"_R")[[1]][2])
				} else model[i,j]<-if(umbral) 2 else max(ncat)+1
			}
		}
	}
	if(corHMM_model){
		mm<-model
		ind<-which(model>0,arr.ind=TRUE)
		for(i in 1:max(ncat)){
			ii<-grep(paste("R",i,sep=""),rownames(model))
			for(j in 1:max(ncat)){
				if(i!=j){
					jj<-grep(paste("R",j,sep=""),rownames(model))
					if(any(mm[ii,jj]>0)){
						kk<-ind[which(ind[,1]%in%ii),,drop=FALSE]
						kk<-kk[which(kk[,2]%in%jj),,drop=FALSE]
						mm[kk]<-min(model[kk])
					}
				}
			}
		}
		old<-sort(unique(as.vector(mm)))
		new<-0:length(old)
		for(i in 1:length(old)) mm[which(mm==old[i],arr.ind=TRUE)]<-new[i]
		model<-mm
	}
	colnames(model)<-rownames(model)<-colnames(X)<-cols
	if(!quiet){
		cat("\nThis is the design matrix of the fitted model.\nDoes it make sense?\n\n")
		print(model)
		cat("\n")
		flush.console()
	}
	fits<-list()
	args<-list(...)
	args$model<-model
	args$tree<-tree
	args$x<-X
	if(hasArg(logscale)) logscale<-rep(list(...)$logscale,niter)
	else logscale<-sample(c(TRUE,FALSE),niter,replace=TRUE)
	if(hasArg(opt.method)) opt.method<-rep(list(...)$opt.method,niter)
	else opt.method<-sample(c("nlminb","optim"),niter,replace=TRUE)
	## parallelize
	if(!parallel){
		for(i in 1:niter){
			args$logscale<-logscale[i]
			args$opt.method<-opt.method[i]
			args$q.init<-rexp(n=max(model),rate=sum(tree$edge.length)/(1e3*k))
			fits[[i]]<-do.call(fitMk,args)
			if(trace>0) print(fits[[i]])
			logL<-sapply(fits,logLik)
			if(!quiet){
				cat(paste("log-likelihood from current iteration:",
					round(logLik(fits[[i]]),4),"\n"))
				cat(paste(" --- Best log-likelihood so far:",round(max(logL),4),
					"---\n"))
				flush.console()
			}
		}
	} else if(parallel){
		if(hasArg(ncores)) ncores<-list(...)$ncores
		else ncores<-detectCores()-1
		mc<-makeCluster(ncores,type="PSOCK")
		registerDoParallel(cl=mc)
		if(!quiet){
			cat(paste("Opened cluster with",ncores,"cores.\n"))
			cat("Running optimization iterations in parallel.\n")
			cat("Please wait....\n")
			flush.console()
		}
		fits<-foreach(i=1:niter)%dopar%{
			args$logscale<-logscale[i]
			args$opt.method<-opt.method[i]
			args$q.init<-rexp(n=max(model),rate=sum(tree$edge.length)/(1e3*k))
			do.call(phytools::fitMk,args)
		}
		logL<-sapply(fits,logLik)
		stopCluster(cl=mc)
	}
	obj<-fits[[which(logL==max(logL))[1]]]
	obj$ncat<-ncat
	obj$model<-MODEL
	obj$umbral<-umbral
	obj$all.fits<-fits
	obj$data<-X
	class(obj)<-c("fitHRM","fitMk")
	obj	
}

## print method for objects of class "fitHRM"
print.fitHRM<-function(x,digits=6,...){
	cat("Object of class \"fitHRM\".\n\n")
	ss<-unique(sapply(x$states,function(x) strsplit(x,"*",fixed=TRUE)[[1]][1]))
	cat(paste("Observed states: [ ",paste(ss,collapse=", ")," ]\n",
		sep=""))
	cat(paste("Number of rate categories per state: [ ",
		paste(x$ncat,collapse=", ")," ]\n\n",sep=""))
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
}

isOdd<-function(x) (x%%2)==1

## S3 plot method for objects of class "fitHRM"
plot.fitHRM<-function(x,...){
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-6
	umbral<-x$umbral
	ncat<-x$ncat
	Q<-as.Qmatrix(x)
	II<-x$index.matrix
	for(i in 1:nrow(II)) for(j in 1:ncol(II)) 
		if(Q[i,j]<tol&&II[i,j]!=0) Q[i,j]<-10*tol
	nn<-sort(rownames(Q))
	mm<-nn
	for(i in 1:length(x$ncat)){
		cs<-if(i>1) sum(x$ncat[1:(i-1)]) else 0
		mm[1:x$ncat[i]+cs]<-if(isOdd(i)) nn[1:x$ncat[i]+cs] else nn[x$ncat[i]:1+cs]
	}
	if(isOdd(length(x$ncat))) mm<-nn
	Q<-Q[mm,mm]
	class(Q)<-"Qmatrix"
	args.list<-list(...)
	args.list$show.zeros<-FALSE
	plot(Q,show.zeros=FALSE,umbral=umbral,ncat=ncat,...)
}

as.Qmatrix.corhmm<-function(x,...){
	Q<-x$solution
	Q[is.na(Q)]<-0
	diag(Q)<--rowSums(Q)
	class(Q)<-"Qmatrix"
	Q
}
