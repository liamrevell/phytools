## functions to do the pruning algorithm and reconstruct ancestral states

joint_asr<-function(q,tree,X,pi,model=NULL,tips=FALSE,tol=1e-12){
	tol<-tol*max(nodeHeights(tree))
	tipn<-1:Ntip(tree)
	tt<-tree
	for(i in 1:length(tipn)) 
		tt<-paintSubTree(tt,tipn[i],state="2",stem=tol)
	tt<-map.to.singleton(tt)
	pw<-reorder(tt,"postorder")
	X<-X[pw$tip.label,]
	k<-ncol(X)
	if(is.null(model)){
		model<-matrix(1,k,k)
		diag(model)<-0
	}
	Q<-matrix(0,k,k)
	Q[]<-c(0,q)[model+1]
	diag(Q)<--rowSums(Q)
	C<-matrix(NA,Nnode(pw)+Ntip(pw),ncol(X),
		dimnames=list(1:(Ntip(pw)+Nnode(pw))))
	colnames(C)<-colnames(X)
	L<-matrix(0,nrow(C),k,dimnames=dimnames(C))
	for(i in 1:Ntip(pw)){
		ind<-which(pw$edge[,2]==i)
		P<-expm(Q*pw$edge.length[ind])
		L[i,]<-log(X[i,]%*%P)
	}
	nn<-unique(pw$edge[,1])
	chk<-FALSE
	for(i in 1:(length(nn)-1)){
		P<-expm(Q*pw$edge.length[which(pw$edge[,2]==nn[i])])
		daughters<-pw$edge[which(pw$edge[,1]==nn[i]),2]
		for(j in 1:k){
			pp<-log(P[j,])+colSums(L[daughters,,drop=FALSE])
			L[nn[i],j]<-max(pp)
			if(length(which(pp==max(pp)))>1) chk<-TRUE
			C[nn[i],j]<-sample(x=names(pp)[which(pp==max(pp))],size=1)
		}
	}
	if(chk){
		cat("Note:\n")
		cat("	 there may be more than one equally likely\n")
		cat("	 joint reconstruction.\n")
	}
	ROOT<-Ntip(pw)+1
	daughters<-pw$edge[which(pw$edge[,1]==ROOT),2]
	L[ROOT,]<-log(pi)+colSums(L[daughters,])
	C[ROOT,]<-names(which(L[ROOT,]==max(L[ROOT,])))
	cw<-reorder(pw)
	nn<-unique(as.vector(cw$edge))
	states<-setNames(rep(NA,length(nn)),nn)
	states[1]<-C[ROOT,1]
	logL<-max(L[ROOT,])
	for(i in 2:length(nn)){
		mother<-cw$edge[which(cw$edge[,2]==nn[i]),1]
		states[i]<-C[nn[i],states[as.character(mother)]]
	}
	M<-matchNodes(tree,tt,method="distances")
	anc<-setNames(states[as.character(M[,2])],M[,1])
	if(tips){
		pp<-sapply(1:Ntip(tt),getParent,tree=tt)
		anc<-c(setNames(states[as.character(pp)],tt$tip.label[1:Ntip(tt)]),
			anc)
	}
	attr(logL,"df")<-length(anc)
	object<-list(ace=as.factor(anc),logLik=logL)
	object
}

## plot method for "ancr" object class
plot.ancr<-function(x,args.plotTree=list(...),args.nodelabels=list(...),...){
	TREE<-attr(x,"tree")
	TYPE<-attr(x,"type")
	args.plotTree$tree<-TREE
	if(is.null(args.plotTree$type))
		args.plotTree$type<-"phylogram"
	if(is.null(args.plotTree$direction)) args.plotTree$direction<-"rightwards"
	if(is.null(args.plotTree$fsize)){
		if(args.plotTree$type%in%c("phylogram","cladogram")){
			if(args.plotTree$direction%in%c("rightwards","leftwards"))
				args.plotTree$fsize<-min(c(6*par()$pin[2]/Ntip(TREE),1))
			else
				args.plotTree$fsize<-min(c(6*par()$pin[1]/Ntip(TREE),1))
		} else {
			args.plotTree$fsize<-min(c(0.6*min(par()$pin)/sqrt(Ntip(TREE)),1))
		}
	}
	if(is.null(args.plotTree$ftype))
		args.plotTree$ftype<-"i"
	if(is.null(args.plotTree$lwd))
		args.plotTree$lwd<-1
	if(is.null(args.plotTree$offset)){
		if(args.plotTree$type%in%c("phylogram","cladogram")){
			if(args.plotTree$direction%in%c("rightwards","leftwards"))
				args.plotTree$offset<-0.4
			else args.plotTree$offset<-1
		} else args.plotTree$offset<-2
	}
	do.call(plotTree,args.plotTree)
	if(TYPE=="marginal"){
		if(nrow(x$ace)==Nnode(TREE)){
			args.nodelabels$pie<-x$ace
			data<-attr(x,"data")
		} else if(nrow(x$ace)==(Ntip(TREE)+Nnode(TREE))){
			args.nodelabels$pie<-x$ace[1:Nnode(TREE)+Ntip(TREE),,drop=FALSE]
			data<-x$ace[1:Ntip(TREE),,drop=FALSE]
		} else cat("Warning: wrong number of rows in input object.\n\n")
	} else if(TYPE=="joint"){
		if(length(x$ace)==Nnode(TREE)){
			data<-attr(x,"data")
			args.nodelabels$pie<-to.matrix(x$ace,colnames(data))
		} else if(nrow(x$ace)==(Ntip(TREE)+Nnode(TREE))){
			args.nodelabels$pie<-to.matrix(x$ace[1:Nnode(TREE)+Ntip(TREE)],
				colnames(attr(x,"data")))
			data<-to.matrix(x$ace[1:Ntip(TREE)],colnames(attr(x,"data")))
		} else cat("Warning: wrong number of elements in input object.\n\n")
	}
	if(is.null(args.nodelabels$cex)){
		if(args.plotTree$type%in%c("phylogram","cladogram")){
			if(args.plotTree$direction%in%c("rightwards","leftwards"))
				args.nodelabels$cex<-min(c(6*par()$pin[2]/Ntip(TREE),0.5))
			else
				args.nodelabels$cex<-min(c(6*par()$pin[1]/Ntip(TREE),0.5))
		} else {
			args.nodelabels$cex<-min(c(0.5*min(par()$pin)/sqrt(Ntip(TREE)),0.5))
		}
	}
	k<-if(TYPE=="marginal") ncol(x$ace) else ncol(attr(x,"data"))
	if(is.null(args.nodelabels$piecol)){
		args.nodelabels$piecol<-palette.colors(k,"Polychrome 36",recycle=TRUE)
		if(k>36) 
			cat("Warning: maximum number of colors reached. Some colors are recycled.\n\n")
	}
	old_fg<-par()$fg
	par(fg="transparent")
	do.call(nodelabels,args.nodelabels)
	if(hasArg(args.tiplabels)) args.tiplabels<-list(...)$args.tiplabels
	else args.tiplabels<-list(...)
	args.tiplabels$piecol<-args.nodelabels$piecol
	args.tiplabels$pie<-data[TREE$tip.label,,drop=FALSE]
	if(is.null(args.tiplabels$cex)){
		if(args.plotTree$type%in%c("phylogram","cladogram")){
			if(args.plotTree$direction%in%c("rightwards","leftwards"))
				args.tiplabels$cex<-min(c(2*par()$pin[2]/Ntip(TREE),0.2))
			else
				args.tiplabels$cex<-min(c(2*par()$pin[1]/Ntip(TREE),0.2))
		} else
			args.tiplabels$cex<-min(c(0.25*min(par()$pin)/sqrt(Ntip(TREE)),0.25))
	}
	do.call(tiplabels,args.tiplabels)
	par(fg=old_fg)
	if(hasArg(legend)) legend<-list(...)$legend
	else legend<-"bottomleft"
	if(legend!=FALSE){
		legend(x=legend,legend=colnames(attr(x,"data")),
			pch=16,col=args.nodelabels$piecol,
			pt.cex=1.2,cex=0.8)
	}
	object<-list(
		fsize=args.plotTree$fsize,
		piecol=args.nodelabels$piecol,
		node_cex=args.nodelabels$cex,
		tip_cex=args.nodelabels$cex,
		legend=legend)
	invisible(object)
}

ancr<-function(object,...) UseMethod("ancr")

ancr.default<-function(object,...){
	warning(paste(
		"ancr does not know how to handle objects of class ",
		class(object),".\n"))
}

logLik.ancr<-function(object,...){
	lik<-object$logLik
	attr(lik,"df")<-attr(object$logLik,"df")
	lik
}

print.ancr<-function(x,digits=6,printlen=6,...){
	Nnode<-Nnode(attr(x,"tree"))
	if(!is.null(printlen)) if(printlen>Nnode) printlen<-Nnode
	if(attr(x,"type")=="marginal"){
		cat("Marginal ancestral state estimates:\n")
		if (is.null(printlen)) 
			print(round(x$ace,digits))
		else {
			ii<-if(nrow(x$ace)>Nnode) 
				1:printlen+Ntip(attr(x,"tree")) else 
				1:printlen
			print(round(x$ace[ii,],digits))
			cat("...\n")
		}
	} else if(attr(x,"type")=="joint"){
		cat("Joint ancestral state estimates:\n")
		tmp<-data.frame(state=x$ace)
		if(is.null(printlen)) print(tmp)
		else {
			ii<-if(nrow(tmp)>Nnode) 
				1:printlen+Ntip(attr(x,"tree")) else 
				1:printlen
			print(tmp[ii,,drop=FALSE])
			cat("...\n")
		}
	}
	cat(paste("\nLog-likelihood =",round(logLik(x),
		digits),"\n\n"))
}

## marginal ancestral states for "anova.fitMk" object
ancr.anova.fitMk<-function(object,...){
	if(hasArg(weighted)) weighted<-list(...)$weighted
	else weighted<-TRUE
	if(hasArg(type)) type<-list(...)$type
	else type<-"marginal"
	if(type!="marginal"){
		cat("\nOnly type=\"marginal\" supported for objects of class \"anova.fitMk\".\n")
		cat("Updating type.\n\n")
		type<-"marginal"
	}
	if(weighted){
		w<-object$weight
		fits<-attr(object,"models")
		anc<-lapply(fits,function(x,...) ancr(x,...)$ace,...)
		ss<-sort(unique(unlist(lapply(anc,colnames))))
		if(any(sapply(attr(object,"models"),class)=="fitHRM")){
				for(i in 1:length(anc)){
						tmp<-anc[[i]]
						anc[[i]]<-matrix(0,nrow(tmp),length(ss),
							dimnames=list(rownames(tmp),ss))
						anc[[i]][rownames(tmp),colnames(tmp)]<-tmp
				}
		}
		anc<-mapply("*",anc,w,SIMPLIFY=FALSE)
		anc<-Reduce("+",anc)
		foo<-function(obj) unclass(as.Qmatrix(obj))
		Q<-lapply(fits,foo)
		if(any(sapply(attr(object,"models"),class)=="fitHRM")){
			for(i in 1:length(Q)){
				tmp<-Q[[i]]
				Q[[i]]<-matrix(0,length(ss),length(ss),
					dimnames=list(ss,ss))
				Q[[i]][rownames(tmp),colnames(tmp)]<-tmp
			}
		}
		Q<-mapply("*",Q,w,SIMPLIFY=FALSE)
		Q<-Reduce("+",Q)
		model<-matrix(0,nrow(Q),ncol(Q))
		k<-nrow(Q)*(nrow(Q)-1)
		model[col(model)!=row(model)]<-1:k
		q<-sapply(1:k,function(i,Q,model) Q[which(model==i)],
			Q=Q,model=model)
		TREE<-fits[[1]]$tree
		DATA<-fits[[1]]$data
		if(any(sapply(attr(object,"models"),class)=="fitHRM")){
			levs<-unique(gsub("*","",ss,fixed=TRUE))
			tmp<-DATA[,levs]
			DATA<-matrix(0,nrow(tmp),length(ss),
				dimnames=list(rownames(tmp),ss))
			for(i in 1:nrow(tmp)){
				for(j in 1:length(levs)){
					DATA[i,grep(levs[j],ss)]<-tmp[i,levs[j]]
				}
			}
		}
		log_lik<-pruning(q,TREE,DATA,model=model)
		attr(log_lik,"df")<-max(model)
		obj<-list(ace=anc,logLik=log_lik)
		attr(obj,"tree")<-TREE
		attr(obj,"data")<-DATA
		attr(obj,"type")<-"marginal"
		class(obj)<-"ancr"
		return(obj)
	} else {
		best<-which(object$AIC==min(object$AIC))
		fits<-attr(object,"models")
		return(ancr(fits[[best]],...))
	}
}

## marginal ancestral states for "fitHRM" object
ancr.fitHRM<-function(object,...) ancr.fitMk(object,...)

hide.hidden<-function(object,...) UseMethod("hide.hidden")

hide.hidden.default<-function(object,...){
	warning(paste(
		"hide.hidden does not know how to handle objects of class ",
		class(object),".\n"))
}

hide.hidden.ancr<-function(object,...){
	ss<-colnames(object$ace)
	ii<-grep("*",ss,fixed=TRUE)
	if(length(ii)>0) ss<-ss[-ii]
	anc<-matrix(0,nrow(object$ace),length(ss),
		dimnames=list(rownames(object$ace),ss))
	for(i in 1:length(ss)){
		anc[,ss[i]]<-rowSums(object$ace[,grep(ss[i],
			colnames(object$ace)),drop=FALSE])
	}
	anc
}

## marginal ancestral states for "fitpolyMk" object
ancr.fitpolyMk<-function(object,...) ancr.fitMk(object,...)

## marginal ancestral states for "fitMk" object
ancr.fitMk<-function(object,...){
	if(hasArg(type)) type<-list(...)$type
	else type<-"marginal"
	if(hasArg(tips)) tips<-list(...)$tips
	else tips<-FALSE
	if(hasArg(lik.func)) lik.func<-list(...)$lik.func
	else lik.func<-"pruning"
	if(lik.func=="parallel"){
		if(hasArg(ncores)) ncores<-list(...)$ncores
		else ncores<-min(nrow(object$tree$edge),detectCores()-1)
		mc<-makeCluster(ncores,type="PSOCK")
		registerDoParallel(cl=mc)
	}
	if(hasArg(expm.method)) expm.method<-list(...)$expm.method
	else expm.method<-"Higham08.b"
	if(hasArg(parallel)) parallel<-list(...)$parallel
	else parallel<-FALSE ## only for lik.func="eigen"
	if(lik.func=="eigen"){
		QQ<-object$index.matrix
		diag(QQ)<--rowSums(QQ)
		eQQ<-eigen(QQ)
	}
	x<-object$data
	tree<-object$tree
	q<-object$rates
	model<-object$index.matrix
	model[is.na(model)]<-0
	pi<-object$pi
	if(type=="marginal"){
		if(lik.func=="pruning"){
			plik<-pruning(q,tree,x,model=model,pi=pi,
				return="conditional",expm.method=expm.method)
		} else if(lik.func=="parallel"){
			plik<-parallel_pruning(q,tree,x,model=model,pi=pi,
				return="conditional",expm.method=expm.method)
		} else if(lik.func=="eigen"){
			plik<-eigen_pruning(q,tree,x,eQQ,pi=pi,return="conditional")
		}
		ace<-marginal_asr(q,tree,plik,model,tips,
			parallel=if((lik.func=="parallel")||parallel) TRUE else FALSE,
			expm.method=expm.method,
			eigen=if(lik.func=="eigen") TRUE else FALSE)
		result<-if(lik.func=="parallel") list(ace=ace,
			logLik=parallel_pruning(q,tree,x,model=model,pi=pi,
			expm.method=expm.method)) else 
			if(lik.func=="pruning") list(ace=ace,logLik=pruning(q,
			tree,x,model=model,pi=pi,expm.method=expm.method)) else 
			if(lik.func=="eigen") list(ace=ace,
			logLik=eigen_pruning(q,tree,x,eigenQ=eQQ))
		if(lik.func=="parallel") stopCluster(mc)
		attr(result$logLik,"df")<-max(model)
		attr(result,"type")<-"marginal"
		attr(result,"tree")<-tree
		attr(result,"data")<-x
		class(result)<-"ancr"
	} else if(type=="joint"){
		if(hasArg(tol)) tol<-list(...)$tol
		else tol<-1e-12
		result<-joint_asr(q,tree,x,pi,model,tips,tol)
		attr(result,"type")<-"joint"
		attr(result,"tree")<-tree
		attr(result,"data")<-x
		class(result)<-"ancr"
	}
	result
}

pruning<-function(q,tree,x,model=NULL,...){
	if(hasArg(return)) return<-list(...)$return
	else return<-"likelihood"
	if(hasArg(expm.method)) expm.method<-list(...)$expm.method
	else expm.method<-"Higham08.b"
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
			P<-expm(Q*pw$edge.length[ee[j]],method=expm.method)
			PP[j,]<-P%*%L[pw$edge[ee[j],2],]
		}
		L[nn[i],]<-apply(PP,2,prod)
		if(nn[i]==root){
			if(pi[1]=="fitzjohn") pi<-L[nn[i],]/sum(L[nn[i],])
			else if(pi[1]=="mle") pi<-as.numeric(L[nn[i],]==max(L[nn[i],]))
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

marginal_asr<-function(q,tree,L,model=NULL,tips=FALSE,
	parallel=FALSE,expm.method="Higham08.b",eigen=FALSE){
	pw<-reorder(tree,"postorder")
	k<-ncol(L)
	if(is.null(model)){
		model<-matrix(1,k,k)
		diag(model)<-0
	}
	Q<-matrix(0,k,k)
	Q[]<-c(0,q)[model+1]
	diag(Q)<--rowSums(Q)
	nn<-unique(pw$edge[,1])
	if(parallel&&(!eigen)){
		P.all<-foreach(i=1:nrow(pw$edge))%dopar%{ 
			expm(Q*pw$edge.length[i],method=expm.method)
		}
	}
	if(eigen){
		QQ<-model
		diag(QQ)<--rowSums(QQ)
		eigQ<-eigen(Q)
		V<-eigQ$vectors
		Vi<-t(V)
		Vals<-eigQ$values
		Expm<-function(t,q) Re(V%*%(exp(q*t*Vals)*Vi))
		if(parallel){
			P.all<-foreach(i=1:nrow(pw$edge))%dopar%{
				Expm(pw$edge.length[i],q)
			}
		} else P.all<-lapply(pw$edge.length,Expm,q=q)
	}
	for(i in length(nn):1){
		ee<-which(pw$edge[,1]==nn[i])
		for(j in 1:length(ee)){
			if(parallel||eigen) P<-P.all[[ee[j]]]
			else P<-expm(Q*pw$edge.length[ee[j]])
			pp<-t(L[nn[i],]/(P%*%L[pw$edge[ee[j],2],]))
			pp[is.nan(pp)]<-0
			L[pw$edge[ee[j],2],]<-(pp%*%P)*
				L[pw$edge[ee[j],2],]
		}
	}
	anc<-L/rep(rowSums(L),k)
	if(tips) anc else anc[1:Nnode(tree)+Ntip(tree),]
}

## marginal ancestral states for "fitgammaMk" object
ancr.fitgammaMk<-function(object,...){
	if(hasArg(type)) type<-list(...)$type
	else type<-"marginal"
	if(hasArg(tips)) tips<-list(...)$tips
	else tips<-FALSE
	x<-object$data
	tree<-object$tree
	q<-object$rates
	alpha<-object$alpha
	nrates<-object$nrates
	model<-object$index.matrix
	model[is.na(model)]<-0
	pi<-object$pi
	if(type=="marginal"){
		plik<-gamma_pruning(c(q,alpha),nrates,tree,x,
			model=model,median=TRUE,pi=pi,return="conditional")
		ace<-marginal_asr_gamma(q,alpha,nrates,tree,plik,
			model,tips)
		result<-list(ace=ace,
			logLik=gamma_pruning(c(q,alpha),nrates,tree,x,
			model=model,pi=pi))
		attr(result$logLik,"df")<-max(model)
		attr(result,"type")<-"marginal"
		attr(result,"tree")<-tree
		attr(result,"data")<-x
		class(result)<-"ancr"
	} else if(type=="joint"){
		## turned off for now
		# if(hasArg(tol)) tol<-list(...)$tol
		# else tol<-1e-12
		# result<-joint_asr(q,tree,x,pi,model,tips,tol)
		# attr(result,"type")<-"joint"
		# attr(result,"tree")<-tree
		# attr(result,"data")<-x
		# class(result)<-"ancr"
		msg<-paste("type = \"joint\" does not yet work",
			"for this object class.")
		stop(msg)
	}
	result
}

marginal_asr_gamma<-function(q,alpha,nrates,tree,L,
	model=NULL,tips=FALSE){
	median<-TRUE
	if(median){
		r<-qgamma(seq(1/(2*nrates),1,by=1/nrates),alpha,alpha)
		r<-r/mean(r)
	} else {
		stop("This does not work yet.\n")
	}
	pw<-reorder(tree,"postorder")
	k<-ncol(L)
	if(is.null(model)){
		model<-matrix(1,k,k)
		diag(model)<-0
	}
	Q<-matrix(0,k,k)
	Q[]<-c(0,q)[model+1]
	diag(Q)<--rowSums(Q)
	nn<-unique(pw$edge[,1])
	for(i in length(nn):1){
		ee<-which(pw$edge[,1]==nn[i])
		for(j in 1:length(ee)){
			P<-Reduce("+",lapply(r,
				function(rr,k,Q,edge) EXPM(Q*rr*edge)/k,
				k=nrates,Q=Q,edge=pw$edge.length[ee[j]]))
			pp<-t(L[nn[i],]/(P%*%L[pw$edge[ee[j],2],]))
			pp[is.nan(pp)]<-0
			L[pw$edge[ee[j],2],]<-(pp%*%P)*
				L[pw$edge[ee[j],2],]
		}
	}
	anc<-L/rep(rowSums(L),k)
	if(tips) anc else anc[1:Nnode(tree)+Ntip(tree),]
}

parallel_pruning<-function(q,tree,x,model=NULL,...){
	if(hasArg(return)) return<-list(...)$return
	else return<-"likelihood"
	if(hasArg(expm.method)) expm.method<-list(...)$expm.method
	else expm.method<-"Higham08.b"
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
	P.all<-foreach(i=1:nrow(pw$edge))%dopar%{ 
		expm(Q*pw$edge.length[i],method=expm.method)
	}
	for(i in 1:length(nn)){
		ee<-which(pw$edge[,1]==nn[i])
		PP<-matrix(NA,length(ee),k)
		for(j in 1:length(ee)){
			P<-P.all[[ee[j]]]
			PP[j,]<-P%*%L[pw$edge[ee[j],2],]
		}
		L[nn[i],]<-apply(PP,2,prod)
		if(nn[i]==root){
			if(pi[1]=="fitzjohn") pi<-L[nn[i],]/sum(L[nn[i],])
			else if(pi[1]=="mle") pi<-as.numeric(L[nn[i],]==max(L[nn[i],]))
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