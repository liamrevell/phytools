# Simulates BM evolution more quickly.
# A trend can be simulated by mu!=0.
# mu=0 is standard BM; mu<0 downward trend; mu>0 upward trend.
# Bounds can be simulated by bounds=c(>-Inf,<Inf).
# OU can be simulated by alpha>0.
# Written by Liam J. Revell 2011, 2013, 2015

fastBM<-function(tree,a=0,mu=0,sig2=1,bounds=c(-Inf,Inf),internal=FALSE,nsim=1,...){
	# some minor error checking
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	## check to see if alpha & theta
	if(hasArg(alpha)) alpha<-list(...)$alpha
	else alpha<-NULL
	if(hasArg(theta)) theta<-list(...)$theta
	else theta<-NULL
	if(!is.null(alpha)&&is.null(theta)){
		cat("Warning: alpha but not theta specified in OU model, setting theta to a.\n")
		theta<-a
	}
	## check for OU w. trend (not permitted)
	if(!is.null(alpha)&&mu!=0)
		cat("Warning: OU with a trend not permitted. Trend parameter will be ignored.\n")
	## check for OU w. bounds (not permitted)
	if(!is.null(alpha)&&(bounds[1]!=-Inf||bounds[2]!=Inf))
		cat("Warning: OU with bounds not permitted. Bounds will be ignored.\n")
	## if BM
	if(is.null(alpha)) x<-simBM(tree,a,mu,sig2,bounds,internal,nsim)
	else x<-if(nsim==1) simOU(tree,alpha,sig2,theta,a,internal) else replicate(nsim,simOU(tree,alpha,sig2,theta,a,internal))
	x
}

## internal function does BM simulation
## written by Liam J. Revell 2011, 2013
simBM<-function(tree,a,mu,sig2,bounds,internal,nsim){
	if(bounds[2]<bounds[1]){
		warning("bounds[2] must be > bounds[1]. Simulating without bounds.")
		bounds<-c(-Inf,Inf)
	}
	if(bounds[1]==-Inf&&bounds[2]==Inf) no.bounds=TRUE
	else no.bounds=FALSE
	if(a<bounds[1]||a>bounds[2]){
		warning("a must be bounds[1]<a<bounds[2]. Setting a to midpoint of bounds.")
		a<-bounds[1]+(bounds[2]-bounds[1])/2
	}
	if(sig2<0){
		warning("sig2 must be > 0.  Setting sig2 to 1.0.")
		sig2=1.0
	}
	# function for reflection off bounds
	reflect<-function(yy,bounds){
		while(yy<bounds[1]||yy>bounds[2]){
			if(yy<bounds[1]) yy<-2*bounds[1]-yy
			if(yy>bounds[2]) yy<-2*bounds[2]-yy
		}
		return(yy)
	}
	# how many species?
	n<-length(tree$tip)
	# first simulate changes along each branch
	x<-matrix(data=rnorm(n=length(tree$edge.length)*nsim,mean=rep(mu*tree$edge.length,nsim),sd=rep(sqrt(sig2*tree$edge.length),nsim)),length(tree$edge.length),nsim)
	# now add them up
	y<-array(0,dim=c(nrow(tree$edge),ncol(tree$edge),nsim))
	for(i in 1:nrow(x)){
		if(tree$edge[i,1]==(n+1))
			y[i,1,]<-a
		else
			y[i,1,]<-y[match(tree$edge[i,1],tree$edge[,2]),2,]

		y[i,2,]<-y[i,1,]+x[i,]
		if(!no.bounds) y[i,2,]<-apply(as.matrix(y[i,2,]),1,function(yy) reflect(yy,bounds))
	}
	rm(x); x<-matrix(data=rbind(y[1,1,],as.matrix(y[,2,])),length(tree$edge.length)+1,nsim)
	rownames(x)<-c(n+1,tree$edge[,2])
	x<-as.matrix(x[as.character(1:(n+tree$Nnode)),])
	rownames(x)[1:n]<-tree$tip.label
	# return simulated data
	if(internal==TRUE)
		return(x[1:nrow(x),]) # include internal nodes
	else
		return(x[1:length(tree$tip.label),]) # tip nodes only
}

## internal function does BM simulation
## written by Liam J. Revell 2013
simOU<-function(tree,alpha,sig2,theta,a0,internal){
	tree<-reorder(tree,"cladewise")
	X<-matrix(0,nrow(tree$edge),ncol(tree$edge))
	root<-length(tree$tip.label)+1
	X[which(tree$edge[,1]==root),1]<-a0
	for(i in 1:nrow(X)){
		t<-tree$edge.length[i]
		s2<-sig2*(1-exp(-2*alpha*t))/(2*alpha)
		X[i,2]<-exp(-alpha*t)*X[i,1]+(1-exp(-alpha*t))*theta+rnorm(n=1,sd=sqrt(s2))
		ii<-which(tree$edge[,1]==tree$edge[i,2])
		if(length(ii)>0) X[ii,1]<-X[i,2]
	}
	x<-sapply(1:max(tree$edge),function(x,y,tree) y[which(tree$edge==x)[1]],y=X,tree=tree)
	x<-setNames(x,c(tree$tip.label,1:tree$Nnode+length(tree$tip.label)))
	if(internal==TRUE)
		return(x) # include internal nodes
	else
		return(x[1:length(tree$tip.label)]) # tip nodes only
}
