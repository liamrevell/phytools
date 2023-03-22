## function creates a stochastic character mapped tree as a modified "phylo" object
## written by Liam Revell 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023

## S3 method for Mk models & objects of various classes

simmap<-function(object,...) UseMethod("simmap")

simmap.default<-function(object,...){
	warning(paste(
		"simmap does not know how to handle objects of class ",
		class(object),".\n"))
}

simmap.simmap<-function(object,...){
	args<-list(...)
	args$tree<-as.phylo(object)
	args$x<-getStates(object,"tips")
	if(hasArg(Q)) Q<-list(...)$Q
	else {
		args$Q<-if(!is.null(object$Q)) object$Q
	}
	args$pi<-if(hasArg(pi)) list(...)$pi else "fitzjohn"
	if(hasArg(nsim)){
		nsim<-list(...)$nsim
		args$nsim<-nsim
	} else args$nsim<-100
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-0
	if(trace>0) args$message<-TRUE
	else args$message<-FALSE
	do.call(make.simmap,args)
}

simmap.fitpolyMk<-function(object,...) simmap.fitMk(object,...)

simmap.fitMk<-function(object,...){
	args<-list(...)
	args$tree<-object$tree
	args$x<-object$data
	args$Q<-as.Qmatrix(object)
	args$pi<-if(object$root.prior=="fitzjohn") 
		"fitzjohn" else object$pi
	if(is.null(args$nsim)) args$nsim<-100
	if(hasArg(trace)) trace<-list(...)$trace
	else trace<-0
	if(trace>0) args$message<-TRUE
	else args$message<-FALSE
	do.call(make.simmap,args)
}

simmap.anova.fitMk<-function(object,...){
	if(hasArg(weighted)) weighted<-list(...)$weighted
	else weighted<-TRUE
	if(hasArg(nsim)) nsim<-list(...)$nsim
	else nsim<-100
	if(weighted){
		w<-object$weight
		mods<-sample(1:length(w),size=nsim,replace=TRUE,prob=w)
	} else {
		best<-which(object$AIC==min(object$AIC))
		mods<-sample(best,size=nsim,replace=TRUE)
	}
	fits<-attr(object,"models")[mods]
	foo<-function(obj,args){
		args$object<-obj
		do.call(simmap,args)
	}
	args<-list(...)
	args$nsim<-1
	tt<-lapply(fits,foo,args=args)
	class(tt)<-c("multiSimmap","multiPhylo")
	return(tt)
}

make.simmap<-function(tree,x,model="SYM",nsim=1,...){
	if(inherits(tree,"multiPhylo")){
		ff<-function(yy,x,model,nsim,...){
			zz<-make.simmap(yy,x,model,nsim,...)
			if(nsim>1) class(zz)<-NULL
			return(zz)
		}	
		if(nsim>1){
			mtrees<-unlist(sapply(tree,ff,x,model,nsim,...,simplify=FALSE),
				recursive=FALSE)
		} else mtrees<-sapply(tree,ff,x,model,nsim,...,simplify=FALSE)
		class(mtrees)<-c("multiSimmap","multiPhylo")
	} else {
		## get optional arguments
		if(hasArg(pi)) pi<-list(...)$pi
		else pi<-"equal"
		if(hasArg(message)) pm<-list(...)$message
		else pm<-TRUE
		if(hasArg(tol)) tol<-list(...)$tol
		else tol<-0
		if(hasArg(Q)) Q<-list(...)$Q
		else Q<-"empirical"
		if(hasArg(burnin)) burnin<-list(...)$burnin
		else burnin<-1000
		if(hasArg(samplefreq)) samplefreq<-list(...)$samplefreq
		else samplefreq<-100
		if(hasArg(vQ)) vQ<-list(...)$vQ
		else vQ<-0.1
		prior<-list(alpha=1,beta=1,use.empirical=FALSE)
		if(hasArg(prior)){ 
			pr<-list(...)$prior
			prior[names(pr)]<-pr
		}
		## done optional arguments
		# check
		if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
		# if vector convert to binary matrix
		if(!is.matrix(x)) xx<-to.matrix(x,sort(unique(x)))
		else xx<-x
		xx<-xx[tree$tip.label,]
		xx<-xx/rowSums(xx)
		# reorder to cladewise
		tree<-bt<-reorder.phylo(tree,"cladewise")
		if(!is.binary(bt)) bt<-multi2di(bt,random=FALSE)
		# some preliminaries
		N<-Ntip(tree)
		m<-ncol(xx)
		root<-N+1
		# get conditional likelihoods & model
		if(is.character(Q)&&Q=="empirical"){				
			XX<-getPars(bt,xx,model,Q=NULL,tree,tol,m,pi=pi,args=list(...))
			L<-XX$L
			Q<-XX$Q
			logL<-XX$loglik
			pi<-XX$pi
			if(pi[1]=="equal") pi<-setNames(rep(1/m,m),colnames(L)) # set equal
			else if(pi[1]=="estimated") pi<-statdist(Q) # set from stationary distribution
			else if(pi[1]=="fitzjohn") pi<-"fitzjohn"
			else pi<-pi/sum(pi) # obtain from input
			if(pm) printmessage(Q,pi,method="empirical")
			mtrees<-replicate(nsim,smap(tree,x,N,m,root,L,Q,pi,logL),simplify=FALSE)
		} else if(is.character(Q)&&Q=="mcmc"){
			if(prior$use.empirical){
				qq<-fitMk(bt,xx,model)$rates
				prior$alpha<-qq*prior$beta
			}
			get.stationary<-if(pi[1]=="estimated") TRUE else FALSE
			if(pi[1]%in%c("equal","estimated"))
				pi<-setNames(rep(1/m,m),colnames(xx)) # set equal
			else if(pi[1]=="fitzjohn") pi<-"fitzjohn"
			else pi<-pi/sum(pi) # obtain from input
			XX<-mcmcQ(bt,xx,model,tree,tol,m,burnin,samplefreq,nsim,vQ,prior,pi=pi)
			L<-lapply(XX,function(x) x$L)
			Q<-lapply(XX,function(x) x$Q)
			logL<-lapply(XX,function(x) x$loglik)
			pi<-if(get.stationary) lapply(Q,statdist) else 
				if(pi[1]=="fitzjohn") lapply(XX,function(x) x$pi) else 
				lapply(1:nsim,function(x,y) y,y=pi)
			if(pm) printmessage(Reduce('+',Q)/length(Q),Reduce('+',pi)/length(pi),
				method="mcmc")
			mtrees<-if(nsim>1) mapply(smap,L=L,Q=Q,pi=pi,logL=logL,MoreArgs=
				list(tree=tree,x=x,N=N,m=m,root=root),SIMPLIFY=FALSE) else
				list(smap(tree=tree,x=x,N=N,m=m,root=root,L=L[[1]],Q=Q[[1]],pi=pi[[1]],
				logL=logL[[1]]))
		} else if(is(Q,"Qmatrix")||is.matrix(Q)){
			if(is(Q,"Qmatrix")) Q<-unclass(Q)
			XX<-getPars(bt,xx,model,Q=Q,tree,tol,m,pi=pi,args=list(...))
			L<-XX$L
			logL<-XX$loglik
			pi<-XX$pi
			if(pi[1]=="equal") pi<-setNames(rep(1/m,m),colnames(L)) # set equal
			else if(pi[1]=="estimated") pi<-statdist(Q) # set from stationary distribution
			else if(pi[1]=="fitzjohn") pi<-"fitzjohn"
			else pi<-pi/sum(pi) # obtain from input
			if(pm) printmessage(Q,pi,method="fixed")
			mtrees<-replicate(nsim,smap(tree,x,N,m,root,L,Q,pi,logL),simplify=FALSE)
		}
		if(length(mtrees)==1) mtrees<-mtrees[[1]]
		else class(mtrees)<-c("multiSimmap","multiPhylo")
	}	
	(if(hasArg(message)) list(...)$message else TRUE)
	if((if(hasArg(message)) list(...)$message else TRUE)&&inherits(tree,"phylo")) message("Done.")
	return(mtrees)
}

# print message
# written by Liam J. Revell 2013
printmessage<-function(Q,pi,method){
	if(method=="empirical"||method=="fixed")
		cat("make.simmap is sampling character histories conditioned on\nthe transition matrix\n\nQ =\n")
	else if(method=="mcmc"){
		cat("make.simmap is simulating with a sample of Q from\nthe posterior distribution\n")
		cat("\nMean Q from the posterior is\nQ =\n")
	}
	print(Q)
	if(method=="empirical") cat("(estimated using likelihood);\n")
	else if(method=="fixed") cat("(specified by the user);\n")
	cat("and (mean) root node prior probabilities\npi =\n")
	if(is.list(pi)) pi<-Reduce("+",pi)/length(pi)
	print(pi)
	flush.console()
}

# mcmc for Q used in Q="mcmc"
# written by Liam J. Revell 2013
mcmcQ<-function(bt,xx,model,tree,tol,m,burnin,samplefreq,nsim,vQ,prior,pi,args=list()){
	update<-function(x){
		x<-abs(x+rnorm(n=np,mean=0,sd=sqrt(vQ)))
		return(x)
	}
	# get model matrix
	if(is.character(model)){
		rate<-matrix(NA,m,m)
		if(model=="ER"){
			np<-rate[]<-1
			diag(rate)<-NA
		}
		if(model=="ARD"){
			np<-m*(m-1)
			rate[col(rate)!=row(rate)]<-1:np
		}
		if(model=="SYM") {
			np<-m*(m-1)/2
			sel<-col(rate)<row(rate)
			rate[sel]<-1:np
			rate<-t(rate)
			rate[sel]<-1:np
		}
	} else {
		if(ncol(model)!=nrow(model)) stop("the matrix given as 'model' is not square")
		if(ncol(model)!=m) 
			stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
		rate<-model
		np<-max(rate)
	}
	# burn-in
	p<-rgamma(np,prior$alpha,prior$beta)
	Q<-matrix(c(0,p)[rate+1],m,m)
	diag(Q)<--rowSums(Q,na.rm=TRUE)
	yy<-getPars(bt,xx,model,Q,tree,tol,m,pi=pi,args=args)
	cat("Running MCMC burn-in. Please wait....\n")
	flush.console()
	for(i in 1:burnin){
		pp<-update(p)
		Qp<-matrix(c(0,pp)[rate+1],m,m)
		diag(Qp)<--rowSums(Qp,na.rm=TRUE)
		zz<-getPars(bt,xx,model,Qp,tree,tol,m,FALSE,pi=pi,args)
		p.odds<-exp(zz$loglik+sum(dgamma(pp,prior$alpha,prior$beta,log=TRUE))-
			yy$loglik-sum(dgamma(p,prior$alpha,prior$beta,log=TRUE)))
		if(p.odds>=runif(n=1)){
			yy<-zz
			p<-pp
		}
	}
	# now run MCMC generation, sampling at samplefreq
	cat(paste("Running",samplefreq*nsim,"generations of MCMC, sampling every",samplefreq,"generations.\nPlease wait....\n\n"))
	flush.console()
	XX<-vector("list",nsim)
	for(i in 1:(samplefreq*nsim)){
		pp<-update(p)
		Qp<-matrix(c(0,pp)[rate+1],m,m)
		diag(Qp)<--rowSums(Qp,na.rm=TRUE)
		zz<-getPars(bt,xx,model,Qp,tree,tol,m,FALSE,pi=pi,args)
		p.odds<-exp(zz$loglik+sum(dgamma(pp,prior$alpha,prior$beta,log=TRUE))-
			yy$loglik-sum(dgamma(p,prior$alpha,prior$beta,log=TRUE)))
		if(p.odds>=runif(n=1)){
			yy<-zz
			p<-pp
		}
		if(i%%samplefreq==0){
			Qi<-matrix(c(0,p)[rate+1],m,m)
			diag(Qi)<--rowSums(Qi,na.rm=TRUE)
			XX[[i/samplefreq]]<-getPars(bt,xx,model,Qi,tree,tol,m,TRUE,pi=pi,args)
		}
	}
	return(XX)
}

# get pars
# written by Liam J. Revell 2013, 2017, 2021
getPars<-function(bt,xx,model,Q,tree,tol,m,liks=TRUE,pi,args=list()){
	if(!is.null(args$pi)) args$pi<-NULL
	args<-c(list(tree=bt,x=xx,model=model,fixedQ=Q,output.liks=liks,pi=pi),args)
	obj<-do.call(fitMk,args)
	N<-length(bt$tip.label)
	pi<-obj$pi
	II<-obj$index.matrix+1
	lvls<-obj$states
	if(liks){
		L<-obj$lik.anc
		rownames(L)<-N+1:nrow(L)
		if(!is.binary(tree)){
			ancNames<-matchNodes(tree,bt)
			L<-L[as.character(ancNames[,2]),]
			rownames(L)<-ancNames[,1]
		}
		L<-rbind(xx,L)
		rownames(L)[1:N]<-1:N
	} else L<-NULL	
	Q<-matrix(c(0,obj$rates)[II],m,m,dimnames=list(lvls,lvls))
	if(any(rowSums(Q,na.rm=TRUE)<tol)){
		message(paste("\nWarning: some rows of Q not numerically distinct from 0; setting to",tol,"\n"))
		ii<-which(rowSums(Q,na.rm=TRUE)<tol)
		for(i in 1:length(ii)) Q[ii[i],setdiff(1:ncol(Q),ii[i])]<-tol/(ncol(Q)-1)
	}
	diag(Q)<--rowSums(Q,na.rm=TRUE)
	return(list(Q=Q,L=L,loglik=logLik(obj),pi=pi))
}

# convert vector of x to binary matrix
# written by Liam J. Revell 2012
to.matrix<-function(x,seq){
	X<-matrix(0,length(x),length(seq),dimnames=list(names(x),seq))
	for(i in 1:length(seq)) X[x==seq[i],i]<-1
	return(X)
}

## function does the stochastic mapping, conditioned on our model & given the conditional likelihoods
## written by Liam J. Revell 2013, 2017
smap<-function(tree,x,N,m,root,L,Q,pi,logL){
	# create the map tree object
	mtree<-tree
	mtree$maps<-list()
	mtree$mapped.edge<-matrix(0,nrow(tree$edge),m,dimnames=list(paste(tree$edge[,1],",",
		tree$edge[,2],sep=""),colnames(L)))
	# now we want to simulate the node states & histories by pre-order traversal
	NN<-matrix(NA,nrow(tree$edge),2) # our node values
	NN[which(tree$edge[,1]==root),1]<-rstate(L[as.character(root),]/
		sum(L[as.character(root),])) # assign root
	for(j in 1:nrow(tree$edge)){
		# conditioned on the start value, assign end value of node (if internal)
		p<-EXPM(Q*tree$edge.length[j])[NN[j,1],]*L[as.character(tree$edge[j,2]),]
		NN[j,2]<-rstate(p/sum(p))
		NN[which(tree$edge[,1]==tree$edge[j,2]),1]<-NN[j,2]
		# now simulate on the branches
		accept<-FALSE	
		while(!accept){
			map<-sch(NN[j,1],tree$edge.length[j],Q)
			if(names(map)[length(map)]==NN[j,2]) accept=TRUE
		}	
		mtree$maps[[j]]<-map
		for(k in 1:length(mtree$maps[[j]])) 
			mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]<-mtree$mapped.edge[j,
				names(mtree$maps[[j]])[k]]+mtree$maps[[j]][k]
	}
	mtree$Q<-Q
	mtree$logL<-logL
	if(!inherits(mtree,"simmap")) class(mtree)<-c("simmap",setdiff(class(mtree),
		"simmap"))
	attr(mtree,"map.order")<-"right-to-left"
	return(mtree)
}

# function generates a history along a branch
# written by Liam J. Revell 2013
sch<-function(start,t,Q){
	tol<-t*1e-12
	dt<-setNames(0,start)
	while(sum(dt)<(t-tol)){
		s<-names(dt)[length(dt)]
		dt[length(dt)]<-if(-Q[s,s]>0) rexp(n=1,rate=-Q[s,s]) else t-sum(dt)
		if(sum(dt)<(t-tol)){
			dt<-c(dt,0)
			if(sum(Q[s,][-match(s,colnames(Q))])>0)
				names(dt)[length(dt)]<-rstate(Q[s,][-match(s,colnames(Q))]/sum(Q[s,][-match(s,colnames(Q))]))
			else names(dt)[length(dt)]<-s
		} else dt[length(dt)]<-dt[length(dt)]-sum(dt)+t
	}
	return(dt)
}

# function uses numerical optimization to solve for the stationary distribution
# written by Liam J. Revell 2013
statdist<-function(Q){
	foo<-function(theta,Q){
		Pi<-c(theta[1:(nrow(Q)-1)],1-sum(theta[1:(nrow(Q)-1)]))
		sum((Pi%*%Q)^2)
	}
	k<-nrow(Q)
	if(nrow(Q)>2){ 
		fit<-optim(rep(1/k,k-1),foo,Q=Q,control=list(reltol=1e-16))
		return(setNames(c(fit$par[1:(k-1)],1-sum(fit$par[1:(k-1)])),rownames(Q)))
	} else {
		fit<-optimize(foo,interval=c(0,1),Q=Q)
		return(setNames(c(fit$minimum,1-fit$minimum),rownames(Q)))
	}
}

## S3 print method for objects of class "simmap" & multiSimmap
## based on print.phylo in ape
print.simmap<-function(x,printlen=6,...){
	N<-Ntip(x)
    	M<-x$Nnode
	cat(paste("\nPhylogenetic tree with",N,"tips and",M,"internal nodes.\n\n"))
	cat("Tip labels:\n")
	if(N>printlen) cat(paste("\t",paste(x$tip.label[1:printlen],collapse=", "),", ...\n",sep=""))
	else print(x$tip.label)
	ss<-sort(unique(c(getStates(x,"tips"),getStates(x,"nodes"))))
	cat(paste("\nThe tree includes a mapped, ",length(ss),"-state discrete character\nwith states:\n",
		sep=""))
	if(length(ss)>printlen) cat(paste("\t",paste(ss[1:printlen],collapse=", "),", ...\n",sep=""))
	else cat(paste("\t",paste(ss,collapse=", "),"\n",sep=""))
    	rlab<-if(is.rooted(x)) "Rooted" else "Unrooted"
	cat("\n",rlab,"; includes branch lengths.\n",sep="")
}
print.multiSimmap<-function(x,details=FALSE,...){
	N<-length(x)
	cat(N,"phylogenetic trees with mapped discrete characters\n")
    	if(details){
		n<-sapply(x,Ntip)
		s<-sapply(x,function(x) length(unique(c(getStates(x,"tips"),getStates(x,"nodes")))))
		for(i in 1:N) cat("tree",i,":",n[i],"tips,",s[i],"mapped states\n")
	}
}

## S3 summary method for objects of class "simmap" & "multiSimmap"
summary.simmap<-function(object,...) describe.simmap(object,...)
summary.multiSimmap<-function(object,...) describe.simmap(object,...)

## for backward compatibility with any function using apeAce internally
apeAce<-function(tree,x,model,fixedQ=NULL,...){
	if(hasArg(output.liks)){ 
		output.liks<-list(...)$output.liks
		return(fitMk(tree,x,model,fixedQ,...))
	} else { 
		output.liks<-TRUE
		return(fitMk(tree,x,model,fixedQ,output.liks=TRUE,...))
	}
}

## S3 logLik methods for "simmap" & "multiSimmap"
logLik.simmap<-function(object,...) object$logL
logLik.multiSimmap<-function(object,...)
	sapply(object,function(x) x$logL)
	
## S3 density method for "multiSimmap"
density.multiSimmap<-function(x,...){
	if(hasArg(method)) method<-list(...)$method
	else method<-"changes"
	if(!method%in%c("densityMap","changes")){
		cat("method not recognized. Setting to default method.\n\n")
		method<-"changes"
	}
	if(method=="densityMap") obj<-densityMap(x,plot=FALSE,...)
	else if(method=="changes"){
		if(hasArg(bw)) bw<-list(...)$bw
		else bw<-1
		tmp<-summary(x)
		ab<-lapply(2:ncol(tmp$count),function(i,x) x[,i],x=tmp$count)
		names(ab)<-sapply(strsplit(colnames(tmp$count)[2:ncol(tmp$count)],
			","),function(x) paste(x,collapse="->"))
		ab<-lapply(ab,function(x){
			class(x)<-"mcmc"
			x })
		## if(.check.pkg("coda")){
			## hpd.ab<-lapply(ab,HPDinterval)
		## } else {
			## cat("  HPDinterval requires package coda.\n")
			## cat("  Computing 95% interval from samples only.\n\n")
			## hpd95<-function(x){
				## obj<-setNames(c(sort(x)[round(0.025*length(x))],
					## sort(x)[round(0.975*length(x))]),
					## c("lower","upper"))
				## attr(obj,"Probability")<-0.95
				## obj
			## }
			## hpd.ab<-lapply(ab,hpd95)
		## }
		hpd.ab<-lapply(ab,HPDinterval)
		minmax<-range(unlist(ab))
		pcalc<-function(x,mm)
			hist(x,breaks=seq(mm[1]-1.5,mm[2]+1.5,bw),plot=FALSE)
		p.ab<-lapply(ab,pcalc,mm=minmax)
		states<-colnames(tmp$ace)
		trans<-names(ab)
		obj<-list(hpd=hpd.ab,
			p=p.ab,
			states=states,trans=trans,
			bw=bw,
			mins=sapply(ab,min),
			meds=sapply(ab,median),
			means=sapply(ab,mean),
			maxs=sapply(ab,max))
		class(obj)<-"changesMap"
	}
	else if(method=="timings"){
		cat("This method doesn't work yet.\n")
		obj<-NULL
	}
	obj
}

## S3 plot method for "changesMap" object from density.multiSimmap
plot.changesMap<-function(x,...){
	if(hasArg(bty)) bty<-list(...)$bty
	else bty<-"l"
	if(hasArg(alpha)) alpha<-list(...)$alpha
	else alpha<-0.3
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-NULL
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-NULL
	if(hasArg(main)) main<-list(...)$main
	else main<-NULL
	if(hasArg(colors)){ 
		colors<-list(...)$colors
		nn<-names(colors)
		colors<-setNames(make.transparent(colors,alpha),nn)
	} else { 
		colors<-if(length(x$trans)==2) 
			setNames(make.transparent(c("blue","red"),alpha),x$trans)
		else
			setNames(rep(make.transparent("blue",alpha),length(x$trans)),
				x$trans)
	}
	if(length(colors)<length(x$trans))
		colors<-rep(colors,ceiling(length(x$trans)/length(colors)))[1:length(x$trans)]
	if(is.null(names(colors))) colors<-setNames(colors,x$trans)
	if(hasArg(transition)){ 
		transition<-list(...)$transition
		if(length(transition)>1){
			cat("transition should be of length 1; truncating to first element.\n")
			transition<-transition[1]
		}
	} else transition<-NULL
	p<-x$p
	hpd<-x$hpd
	bw<-x$bw
	if(length(x$trans)==2&&is.null(transition)){
		plot(p[[1]]$mids,p[[1]]$density,xlim=if(is.null(xlim)) 
			c(min(x$mins)-1,max(x$maxs)+1) else xlim,
			ylim=if(is.null(ylim)) c(0,1.2*max(c(p[[1]]$density,
			p[[2]]$density))) else ylim,
			type="n",xlab="number of changes",
			ylab="relative frequency across stochastic maps",
			bty=bty)
		y2<-rep(p[[1]]$density,each=2)
		y2<-y2[-length(y2)]
		x2<-rep(p[[1]]$mids-bw/2,each=2)[-1]
		x3<-c(min(x2),x2,max(x2))
		y3<-c(0,y2,0)
		polygon(x3,y3,col=colors[x$trans[1]],border=FALSE)
		lines(p[[1]]$mids-bw/2,p[[1]]$density,type="s")
		y2<-rep(p[[2]]$density,each=2)
		y2<-y2[-length(y2)]
		x2<-rep(p[[2]]$mids-bw/2,each=2)[-1]
		x3<-c(min(x2),x2,max(x2))
		y3<-c(0,y2,0)
		polygon(x3,y3,col=colors[x$trans[2]],border=FALSE)
		lines(p[[2]]$mids-bw/2,p[[2]]$density,type="s")
		dd<-0.01*diff(par()$usr[3:4])
		lines(hpd[[1]],rep(max(p[[1]]$density)+dd,2))
		lines(rep(hpd[[1]][1],2),c(max(p[[1]]$density)+dd,
			max(p[[1]]$density)+dd-0.005))
		lines(rep(hpd[[1]][2],2),c(max(p[[1]]$density)+dd,
			max(p[[1]]$density)+dd-0.005))
		CHARS<-strsplit(x$trans[1],"->")[[1]]
		CHARS[1]<-paste("HPD(",CHARS[1],collapse="")
		CHARS[2]<-paste(CHARS[2],")",collapse="")
		T1<-bquote(.(CHARS[1])%->%.(CHARS[2]))
		text(mean(hpd[[1]]),max(p[[1]]$density)+dd,
			T1,pos=3)
		lines(hpd[[2]],rep(max(p[[2]]$density)+dd,2))
		lines(rep(hpd[[2]][1],2),c(max(p[[2]]$density)+dd,
			max(p[[2]]$density)+dd-0.005))
		lines(rep(hpd[[2]][2],2),c(max(p[[2]]$density)+dd,
			max(p[[2]]$density)+dd-0.005))
		CHARS<-strsplit(x$trans[2],"->")[[1]]
		CHARS[1]<-paste("HPD(",CHARS[1],collapse="")
		CHARS[2]<-paste(CHARS[2],")",collapse="")
		T2<-bquote(.(CHARS[1])%->%.(CHARS[2]))
		text(mean(hpd[[2]]),max(p[[2]]$density)+dd,
			T2,pos=3)
		CHARS<-strsplit(x$trans[1],"->")[[1]]
		T1<-bquote(.(CHARS[1])%->%.(CHARS[2]))
		CHARS<-strsplit(x$trans[2],"->")[[1]]
		T2<-bquote(.(CHARS[1])%->%.(CHARS[2]))
		legend("topleft",legend=c(T1,T2),pch=22,pt.cex=2.2,bty="n",
			pt.bg=colors[x$trans])
	} else {
		k<-if(is.null(transition)) length(x$states) else 1
		if(k>1) par(mfrow=c(k,k))
		ii<-if(is.null(transition)) 1 else which(x$trans==transition)
		max.d<-max(unlist(lapply(p,function(x) x$density)))
		for(i in 1:k){
			for(j in 1:k){
				if(i==j&&is.null(transition)) plot.new()
				else {
					CHARS<-strsplit(x$trans[ii],"->")[[1]]
					MAIN<-if(is.null(main)) bquote(.(CHARS[1])%->%.(CHARS[2])) else
						main
					plot(p[[ii]]$mids,p[[ii]]$density,xlim=if(is.null(xlim)) 
						c(min(x$mins)-1,max(x$maxs)+1) else xlim,
						ylim=if(is.null(ylim)) c(0,1.2*max.d) else ylim,
						type="n",xlab="number of changes",
						ylab="relative frequency",main=MAIN,font.main=1,
						bty=bty)
					y2<-rep(p[[ii]]$density,each=2)
					y2<-y2[-length(y2)]
					x2<-rep(p[[ii]]$mids-bw/2,each=2)[-1]
					x3<-c(min(x2),x2,max(x2))
					y3<-c(0,y2,0)
					polygon(x3,y3,col=colors[x$trans[ii]],border=FALSE)
					lines(p[[ii]]$mids-bw/2,p[[ii]]$density,type="s")
					dd<-0.03*diff(par()$usr[3:4])
					lines(hpd[[ii]],rep(max(p[[ii]]$density)+dd,2))
					text(mean(hpd[[ii]]),max(p[[ii]]$density)+dd,"HPD",pos=3)
					ii<-ii+1
				}
			}
		}
	}
}

print.changesMap<-function(x, ...){
	if(hasArg(signif)) signif<-list(...)$signif
	else signif<-2
	cat("\nDistribution of changes from stochastic mapping:\n")
	NROW<-ceiling(length(x$trans)/2)
	if(NROW>1) cat("\n")
	for(i in 1:NROW){
		ii<-2*i-1
		jj<-ii+1
		cat(paste("\t",x$trans[ii],"\t\t",x$trans[jj],"\n",sep=""))
		cat(paste("\tMin.   :",round(x$mins[ii],signif),
			"\tMin.   :",round(x$mins[jj],signif),"\n",sep=""))
		cat(paste("\tMedian :",round(x$meds[ii],signif),
			"\tMedian :",round(x$meds[jj],signif),"\n",sep=""))
		cat(paste("\tMean   :",round(x$means[ii],signif),
			"\tMean   :",round(x$means[jj],signif),"\n",sep=""))
		cat(paste("\tMax.   :",round(x$maxs[ii],signif),
			"\tMax.   :",round(x$maxs[jj],signif),"\n\n",sep=""))
	}
	for(i in 1:length(x$trans))
		cat("95% HPD interval(",x$trans[i],"): [",x$hpd[[i]][1],", ",
			x$hpd[[i]][2],"]\n",sep="")
	cat("\n")
}

