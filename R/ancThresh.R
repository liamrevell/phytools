## function performs ancestral character estimation under the threshold model
## written by Liam J. Revell 2012, 2013, 2014, 2017

ancThresh<-function(tree,x,ngen=10000,sequence=NULL,method="mcmc",model=c("BM","OU","lambda"),control=list(),...){

	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	
	# check method
	if(method!="mcmc") stop(paste(c("do not recognize method =",method,",quitting")))

	# get model for the evolution of liability
	model<-model[1]

	# check x
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)){
		X<-x[tree$tip.label,]
		if(is.null(sequence)){
			message("**** NOTE: no sequence provided, column order in x")
 			seq<-colnames(X)
		} else seq<-sequence
	} else if(is.vector(x)){
		x<-x[tree$tip.label]
		if(is.null(sequence)){
			message("**** NOTE: no sequence provided, using alphabetical or numerical order")
 			seq<-sort(levels(as.factor(x)))
		} else seq<-sequence
		X<-to.matrix(x,seq)
	}
	# row scale X
	X<-X/apply(X,1,sum)
	X<-X[,seq] # order columns by seq

	# ok, now set starting thresholds
	th<-c(1:length(seq))-1; names(th)<-seq
	x<-to.vector(X)
	l<-sapply(x,function(x) runif(n=1,min=th[x]-1,max=th[x])) # set plausible starting liability
	if(model=="OU") alpha<-0.1*max(nodeHeights(tree))
	if(model=="lambda") lambda<-1.0

	# for MCMC
	n<-length(tree$tip)
	m<-length(th)
	npar<-if(model=="BM") tree$Nnode+n+m-2 else tree$Nnode+n+m-1
	

	# populate control list
	PrA<-matrix(1/m,tree$Nnode,m,dimnames=list(1:tree$Nnode+n,seq))
	if(!is.null(control$pr.anc)){
		if(!is.matrix(control$pr.anc)){
			message("**** NOTE: prior on ancestral states must be in matrix form; using default prior")
			control$pr.anc<-NULL
		} else {
			control$pr.anc<-control$pr.anc[,seq,drop=FALSE]
			PrA[rownames(control$pr.anc),]<-control$pr.anc
			control$pr.anc<-PrA	
		}
	}
	con=list(sample=1000,
		propliab=0.5*max(nodeHeights(tree)),
		propthresh=0.05*max(nodeHeights(tree)),
		propalpha=0.1*max(nodeHeights(tree)),
		proplambda=0.01,
		pr.anc=PrA,
		pr.th=0.01,
		burnin=round(0.2*ngen),
		plot=FALSE,
		print=TRUE,
		piecol=setNames(palette()[1:length(seq)],seq),
		tipcol="input",
		quiet=FALSE)
	con[(namc<-names(control))]<-control
	con<-con[!sapply(con,is.null)]

	# now set ancestral liabilities, by first picking ancestral states from their prior
	temp<-apply(con$pr.anc,1,rstate)
	# assign random liabilities consistent with the starting thresholds
	a<-sapply(temp,function(x) runif(n=1,min=th[x]-1,max=th[x]))

	# now change the upper limit of th to Inf
	th[length(th)]<-Inf

	# compute some matrices & values
	V<-if(model=="BM") vcvPhylo(tree) 
	   else if(model=="OU") vcvPhylo(tree,model="OU",alpha=alpha) 
	   else if(model=="lambda") vcvPhylo(tree,model="lambda",lambda=lambda)
	# check to make sure that V will be non-singular
	if(any(tree$edge.length<=(10*.Machine$double.eps)))
		stop("some branch lengths are 0 or nearly zero")
	invV<-solve(V)
	detV<-determinant(V,logarithm=TRUE)$modulus[1]
	lik1<-likLiab(l,a,V,invV,detV)+log(probMatch(X,l,th,seq))
	
	# store
	A<-matrix(NA,ngen/con$sample+1,tree$Nnode,dimnames=list(NULL,n+1:tree$Nnode))
	B<-if(model=="BM") matrix(NA,ngen/con$sample+1,m+2,dimnames=list(NULL,c("gen",names(th),"logLik")))
	   else if(model=="OU") matrix(NA,ngen/con$sample+1,m+3,dimnames=list(NULL,c("gen",names(th),"alpha","logLik")))
	   else if(model=="lambda") matrix(NA,ngen/con$sample+1,m+3,dimnames=list(NULL,c("gen",names(th),"lambda","logLik")))

	C<-matrix(NA,ngen/con$sample+1,tree$Nnode+n,dimnames=list(NULL,c(tree$tip.label,1:tree$Nnode+n)))
	A[1,]<-sapply(a,threshState,thresholds=th)
	B[1,]<-if(model=="BM") c(0,th,lik1) else if(model=="OU") c(0,th,alpha,lik1) else if(model=="lambda") c(0,th,lambda,lik1)
	C[1,]<-c(l[tree$tip.label],a[as.character(1:tree$Nnode+n)])

	# run MCMC
	message("MCMC starting....")
	logL<-lik1<-likLiab(l,a,V,invV,detV)+log(probMatch(X,l,th,seq))
	for(i in 1:ngen){
		lik1<-logL
		d<-i%%npar
		if(ngen>=1000) if(i%%1000==0) if(con$print) message(paste("gen",i))
		ap<-a; lp<-l; thp<-th
		if(model=="OU") alphap<-alpha
		if(model=="lambda") lambdap<-lambda
		Vp<-V; invVp<-invV; detVp<-detV
		if(d<=tree$Nnode&&d!=0){
			# update node liabilities
			ind<-d%%tree$Nnode; if(ind==0) ind<-tree$Nnode
			ap[ind]<-a[ind]+rnorm(n=1,sd=sqrt(con$propliab))
		} else {
			if((d>tree$Nnode&&d<=(tree$Nnode+n))||(npar==(tree$Nnode+n)&&d==0)){
				# update tip liabilities
				if(d==0) ind<-n
				else { ind<-(d-tree$Nnode)%%n; if(ind==0) ind<-n }
				lp[ind]<-l[ind]+rnorm(n=1,sd=sqrt(con$propliab))
			} else if(d>(tree$Nnode+n)&&d<=(tree$Nnode+n+m-2)||(npar==(tree$Nnode+n+m-2)&&d==0)) {
				# update thresholds
				if(d) ind<-(d-tree$Nnode-n)%%m+1
				else ind<-m-1
				thp[ind]<-bounce(th[ind],rnorm(n=1,sd=sqrt(con$propthresh)),c(th[ind-1],th[ind+1]))
			} else {
				if(model=="OU"){
					alphap<-bounce(alpha,rnorm(n=1,sd=sqrt(con$propalpha)),c(0,Inf))
					Vp<-vcvPhylo(tree,model="OU",alpha=alphap)
				} else if(model=="lambda"){
					lambdap<-bounce(lambda,rnorm(n=1,sd=sqrt(con$proplambda)),c(0,1))
					Vp<-vcvPhylo(tree,model="lambda",lambda=lambdap)
				}
				invVp<-solve(Vp)
				detVp<-determinant(Vp,logarithm=TRUE)$modulus[1]
			}
		}
		lik2<-likLiab(lp,ap,Vp,invVp,detVp)+log(probMatch(X,lp,thp,seq))
		p.odds<-min(c(1,exp(lik2+logPrior(sapply(ap,threshState,thresholds=thp),thp,con)-lik1-logPrior(sapply(a,threshState,thresholds=th),th,con))))

		if(p.odds>runif(n=1)){
			a<-ap; l<-lp; th<-thp
			V<-Vp; detV<-detVp; invV<-invVp
			if(model=="OU") alpha<-alphap
			if(model=="lambda") lambda<-lambdap
			logL<-lik2
		} else logL<-lik1
		if(i%%con$sample==0){ 
			A[i/con$sample+1,]<-sapply(a,threshState,thresholds=th)
			B[i/con$sample+1,]<-if(model=="BM") c(i,th[colnames(B)[1+1:m]],logL) else if(model=="OU") c(i,th[colnames(B)[1+1:m]],alpha,logL) else if(model=="lambda") c(i,th[colnames(B)[1+1:m]],lambda,logL)
			C[i/con$sample+1,]<-c(l[tree$tip.label],a[as.character(1:tree$Nnode+n)])
		}
	}
	mcmc<-as.data.frame(A)
	param<-as.data.frame(B)
	liab<-as.data.frame(C)
	ace<-matrix(0,tree$Nnode,m,dimnames=list(colnames(A),seq))
	burnin<-which(param[,"gen"]==con$burnin)
	for(i in 1:tree$Nnode){
		temp<-summary(mcmc[burnin:nrow(mcmc),i])/(nrow(mcmc)-burnin+1)
		ace[i,names(temp)]<-temp
	}
	obj<-list(ace=ace,mcmc=mcmc,par=param,liab=liab,
		tree=tree,x=x,model=model,
		seq=seq,
		ngen=ngen,sample=con$sample,
		burnin=con$burnin)
	class(obj)<-"ancThresh"
	if(con$plot) plot(obj)
	obj
}

## some S3 methods (added in 2017)

print.ancThresh<-function(x,...){
	cat("\nObject containing the results from an MCMC analysis\nof the threshold model using ancThresh.\n\n")
	cat("List with the following components:\n")
	cat(paste("ace:\tmatrix with posterior probabilities assuming",x$burnin,
		"\n\tburn-in generations.\n"))
	cat("mcmc:\tposterior sample of liabilities at tips & internal\n")
	cat(paste("\tnodes (a matrix with",nrow(x$mcmc),"rows &",ncol(x$mcmc),"columns).\n"))
	cat("par:\tposterior sample of the relative positions of the\n")
	cat(paste("\tthresholds, the log-likelihoods, and any other\n",
		"\tmodel variables (a matrix with",nrow(x$par),"rows).\n\n"))
	cat("The MCMC was run under the following conditions:\n")
	cat(paste("\tseq =",paste(x$seq,collapse=" <-> "),
		"\n\tmodel =",x$model,"\n\tnumber of generations =",x$ngen,
		"\n\tsample interval=",x$sample,
		"\n\tburn-in =",x$burnin,"\n\n"))
}

plot.ancThresh<-function(x,...){
	if(hasArg(burnin)) burnin<-list(...)$burnin 
	else burnin<-x$burnin
	args<-list(...)
	if(is.null(args$lwd)) args$lwd<-1
	if(is.null(args$ylim)) args$ylim<-c(-0.1*Ntip(x$tree),Ntip(x$tree))
	if(is.null(args$offset)) args$offset<-0.5
	if(is.null(args$ftype)) args$ftype="i"
	args$tree<-x$tree	
	do.call(plotTree,args)
	ii<-which(x$par[,1]==burnin)+1
	LIAB<-as.matrix(x$liab)[ii:nrow(x$liab),]
	THRESH<-as.matrix(x$par)[ii:nrow(x$par),1:length(x$seq)+1]
	STATES<-matrix(NA,nrow(LIAB),ncol(LIAB),dimnames=dimnames(LIAB))
	for(i in 1:nrow(LIAB)) STATES[i,]<-threshState(LIAB[i,],THRESH[i,])
	PP<-t(apply(STATES,2,function(x,levs) summary(factor(x,levels=levs))/length(x),
		levs=x$seq))
	if(hasArg(piecol)) piecol<-list(...)$piecol
	else piecol<-setNames(colorRampPalette(c("blue",
		"yellow"))(length(x$seq)),x$seq)
	if(hasArg(node.cex)) node.cex<-list(...)$node.cex
	else node.cex<-0.6
	nodelabels(pie=PP[1:x$tree$Nnode+Ntip(x$tree),],
		piecol=piecol,cex=node.cex)
	if(hasArg(tip.cex)) tip.cex<-list(...)$tip.cex
	else tip.cex<-0.4
	tiplabels(pie=PP[x$tree$tip.label,],piecol=piecol,
		cex=tip.cex)
	legend(x=par()$usr[1],y=par()$usr[1],legend=x$seq,pch=21,pt.bg=piecol,
		pt.cex=2.2,bty="n")
}

# plots ancestral states from the threshold model
# written by Liam J. Revell 2012, 2014

plotThresh<-function(tree,x,mcmc,burnin=NULL,piecol,tipcol="input",legend=TRUE,...){

	if(is.logical(legend)||is.vector(legend)){
		if(is.logical(legend)&&legend==TRUE) leg<-setNames(names(piecol),names(piecol))
		else if(is.vector(legend)){ 
			leg<-legend[names(piecol)]
			legend<-TRUE
		}
	}

	# plot tree
	par(lend=2)
	plotTree(tree,ftype="i",lwd=1,ylim=if(legend) c(-0.1*length(tree$tip.label),length(tree$tip.label)) else NULL,...)
	if(legend){
		zz<-par()$cex; par(cex=0.6)
		for(i in 1:length(piecol))
			add.simmap.legend(leg=leg[i],colors=piecol[i],prompt=FALSE,x=0.02*max(nodeHeights(tree)),y=-0.1*length(tree$tip.label),vertical=FALSE,shape="square",fsize=1)
		par(cex=zz)
	}
	# pull matrices from mcmc
	ace<-mcmc$ace
	liab<-mcmc$liab
	param<-mcmc$par

	# get burnin
	if(is.null(burnin)) burnin<-round(0.2*max(param[,"gen"]))
	burnin<-which(param[,"gen"]==burnin)

	# check x
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)) X<-x[tree$tip.label,]
	else if(is.vector(x)){
		x<-x[tree$tip.label]
		X<-to.matrix(x,names(piecol))
	}
	# row scale X
	X/apply(X,1,sum)->X

	# plot node labels
	nodelabels(pie=ace,piecol=piecol[colnames(ace)],cex=0.6)

	# plot tip labels
	if(tipcol=="input") tiplabels(pie=X,piecol=piecol[colnames(X)],cex=0.6)
	else if(tipcol=="estimated") {
		XX<-matrix(NA,nrow(liab),length(tree$tip),dimnames=list(rownames(liab),colnames(liab)[1:length(tree$tip)]))
		for(i in 1:nrow(liab)) XX[i,]<-sapply(liab[i,1:length(tree$tip)],threshState,thresholds=param[i,1:ncol(X)+1])
		X<-t(apply(XX,2,function(x) summary(factor(x,levels=colnames(X)))))
		tiplabels(pie=X/rowSums(X),piecol=piecol[colnames(X)],cex=0.6)
	}
}

# computes DIC for threshold model
# written by Liam J. Revell 2012, 2014

threshDIC<-function(tree,x,mcmc,burnin=NULL,sequence=NULL,method="pD"){
	## identify model
	if(any(colnames(mcmc$par)=="alpha")) model<-"OU"
	else if(any(colnames(mcmc$par)=="lambda")) model<-"lambda"
	else model<-"BM"
	# check x
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)){
		X<-x[tree$tip.label,]
		if(is.null(sequence)){
			message("**** NOTE: no sequence provided, column order in x")
 			seq<-colnames(X)
		} else seq<-sequence
	} else if(is.vector(x)){
		x<-x[tree$tip.label]
		if(is.null(sequence)){
			message("**** NOTE: no sequence provided, using alphabetical or numerical order")
 			seq<-sort(levels(as.factor(x)))
		} else seq<-sequence
		X<-to.matrix(x,seq)
	}
	# row scale X
	X<-X/apply(X,1,sum)
	X<-X[,seq] # order columns by seq
	# convert burnin to starting row
	if(is.null(burnin)) burnin<-0.2*max(mcmc$par[,"gen"])
	start<-which(mcmc$par[,"gen"]==burnin)+1
	# compute
	k<-if(model=="BM") 1 else 2
	thBar<-colMeans(mcmc$par[start:nrow(mcmc$par),2:(ncol(mcmc$par)-k)])
	liabBar<-colMeans(mcmc$liab[start:nrow(mcmc$liab),])
	if(model=="BM")	V<-vcvPhylo(tree)
	else if(model=="OU") V<-vcvPhylo(tree,model="OU",alpha=mean(mcmc$par[start:nrow(mcmc$par),"alpha"]))
	else if(model=="lambda") V<-vcvPhylo(tree,model="lambda",lambda=mean(mcmc$par[start:nrow(mcmc$par),"lambda"]))
	Dtheta<--2*(likLiab(liabBar[tree$tip.label],liabBar[as.character(1:tree$Nnode+length(tree$tip))],V,solve(V),determinant(V,logarithm=TRUE)$modulus[1])+log(probMatch(X[tree$tip.label,],liabBar[tree$tip.label],thBar,seq)))
	D<--2*mcmc$par[start:nrow(mcmc$par),"logLik"]
	Dbar<-mean(D)
	if(method=="pD"){		
		pD<-Dbar-Dtheta
		DIC<-pD+Dbar
		result<-setNames(c(Dbar,Dtheta,pD,DIC),c("Dbar","Dhat","pD","DIC"))
	} else if(method=="pV"){
		pV<-var(D)/2
		DIC<-pV+Dbar
		result<-setNames(c(Dbar,Dtheta,pV,DIC),c("Dbar","Dhat","pV","DIC"))
	}
	return(result)
}

# internal functions for ancThresh, plotThresh, and threshDIC

## returns a state based on position relative to thresholds
## threshStateC is a function from phangorn>=2.3.1
threshState<-if(packageVersion("phangorn")>='2.3.1'){
	function(x,thresholds) names(thresholds)[threshStateC(x,thresholds)]
} else function(x,thresholds){
	t<-c(-Inf,thresholds,Inf)
	names(t)[length(t)]<-names(t)[length(t)-1] 
	i<-1
	while(x>t[i]) i<-i+1
	names(t)[i]
}

# likelihood function for the liabilities
likLiab<-function(l,a,V,invV,detV){
	y<-c(l,a[2:length(a)])-a[1]
	logL<--y%*%invV%*%y/2-nrow(V)*log(2*pi)/2-detV/2
	return(logL)
}

# function for the log-prior
logPrior<-function(a,t,control){
	pp<-sum(log(diag(control$pr.anc[names(a),a])))+
		if(length(t)>2) sum(dexp(t[2:(length(t)-1)],rate=control$pr.th,log=TRUE)) else 0				
	return(pp)		
}

# check if the liability predictions match observed data
allMatch<-function(x,l,thresholds){
	result<-all(sapply(l,threshState,thresholds=thresholds)==x)
	if(!is.na(result)) return(result)
	else return(FALSE)
}

# check if the liability predictions match observed data & return a probability
# (this allows states to be uncertain)
probMatch<-function(X,l,thresholds,sequence){
	Y<-to.matrix(sapply(l,threshState,thresholds=thresholds),sequence)
	return(prod(rowSums(X*Y)))
}

# bounds parameter by bouncing
bounce<-function(start,step,bounds){
	x<-start+step
	while(x>bounds[2]||x<bounds[1]){
		if(x>bounds[2]) x<-2*bounds[2]-x
		if(x<bounds[1]) x<-2*bounds[1]-x
	}
	return(x)
}

# convert vector of x to binary matrix
to.matrix<-function(x,seq){
	X<-matrix(0,length(x),length(seq),dimnames=list(names(x),seq))
	for(i in 1:length(seq)) X[x==seq[i],i]<-1
	return(X)
}

# convert binary matrix to vector
to.vector<-function(X) apply(X,1,rstate)
