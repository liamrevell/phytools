# implements method of Ives et al. 2007 for PGLS regression with sampling error
# written by Liam J. Revell 2012, 2013, 2015

pgls.Ives<-function(tree,X,y,Vx=NULL,Vy=NULL,Cxy=NULL,lower=c(1e-8,1e-8)){

	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	
	# likelihood function
	lik<-function(theta,C,x,y,Mx,My,Mxy){
		sig2x<-theta[1]; sig2y<-theta[2]; b1<-theta[3]
		a<-theta[4:5]
		n<-nrow(C)
		Psi<-matrix(0,2*n,2*n)
		Psi[1:n,1:n]<-sig2x*C+diag(Mx)
		Psi[n+1:n,1:n]<-Psi[1:n,n+1:n]<-b1*sig2x*C+diag(Mxy)
		Psi[n+1:n,n+1:n]<-b1^2*sig2x*C+sig2y*C+diag(My)
		z<-c(X,y)
		D<-kronecker(diag(rep(1,2)),matrix(rep(1,n)))
		## L<--2*n/2*log(2*pi)-(1/2)*determinant(Psi,logarithm=TRUE)$modulus[1]-(1/2)*t(z-D%*%a)%*%solve(Psi)%*%(z-D%*%a)
		L<-dmnorm(z,(D%*%a)[,1],Psi,log=TRUE)
		return(-L)
	}

	# check data input format
	Xbar<-ybar<-NULL
	if((length(X)>length(unique(names(X))))&&is.null(Vx)){
		a<-aggregate(X,by=list(names(X)),FUN=mean)
		Xbar<-setNames(a[,2],a[,1])
		a<-aggregate(X,by=list(names(X)),FUN=var)
		Vx<-setNames(a[,2],a[,1])
		nx<-summary(as.factor(names(X)))[names(Vx)]
		Vx<-Vx/nx
		if(any(is.na(Vx))){
			warning("Some species contain only one sample. Substituting mean variance.",call.=FALSE)
			Vx[which(is.na(Vx))]<-mean(Vx*nx,na.rm=TRUE)
		}
		rm(a,nx)
	}
	if((length(y)>length(unique(names(y))))&&is.null(Vy)){
		a<-aggregate(y,by=list(names(y)),FUN=mean)
		ybar<-setNames(a[,2],a[,1])
		a<-aggregate(y,by=list(names(y)),FUN=var)
		Vy<-setNames(a[,2],a[,1])
		ny<-summary(as.factor(names(y)))[names(Vy)]
		Vy<-Vy/ny
		if(any(is.na(Vy))){
			warning("Some species contain only one sample. Substituting mean variance.",call.=FALSE)
			Vy[which(is.na(Vy))]<-mean(Vy*ny,na.rm=TRUE)
		}
		rm(a,ny)
	}
	if(is.null(Cxy)){
		a<-aggregate(X,by=list(names(X)),FUN=mean); a<-setNames(a[,2],a[,1])
		b<-aggregate(y,by=list(names(y)),FUN=mean); b<-setNames(b[,2],b[,1])
		c<-(X-a[names(X)])*(y-b[names(y)])
		d<-aggregate(c,by=list(names(c)),FUN=sum); d<-setNames(d[,2],d[,1])
		nxy<-summary(as.factor(names(y)))[names(d)]
		Cxy<-d/(nxy-1)/nxy
		if(any(is.na(Cxy))){
			warning("Some species contain only one sample. Substituting mean variance.",call.=FALSE)
			Cxy[which(is.na(Cxy))]<-mean(Cxy*nxy,na.rm=TRUE)
		}
		rm(a,b,c,nxy)
	}
	if(!is.null(Xbar)) X<-Xbar
	if(!is.null(ybar)) y<-ybar

	# perform calculation & organization
	C<-vcv.phylo(tree)
	X<-X[tree$tip.label]
	y<-y[tree$tip.label]
	Vx<-Vx[tree$tip.label]
	Vy<-Vy[tree$tip.label]
	Cxy<-Cxy[tree$tip.label]
	
	# get some reasonable starting values for optimization
	b<-runif(n=1,min=0,max=2)*lm(pic(y,tree)~pic(X,tree))$coefficients[2]; names(b)<-NULL
	sig2x<-runif(n=1,min=0,max=2)*mean(pic(X,tree)^2)
	sig2y<-runif(n=1,min=0,max=2)*mean(pic(y,tree)^2)
	a<-runif(n=2,min=-1,max=1)*c(mean(X),mean(y))

	# optimize regression model
	r<-optim(c(sig2x,sig2y,b,a),lik,C=C,x=X,y=y,Mx=Vx,My=Vy,Mxy=Cxy,method="L-BFGS-B",lower=c(lower,-Inf,-Inf,-Inf),control=list(factr=1e10))

	# return r
	return(list(beta=c(r$par[5]-r$par[3]*r$par[4],r$par[3]),sig2x=r$par[1],sig2y=r$par[2],a=r$par[4:5],logL=-r$value,convergence=r$convergence,message=r$message))
}

## simpler function to take sampling error into account for y only
## written by Liam J. Revell 2017

pgls.SEy<-function(model,data,corClass=corBrownian,tree=tree,
	se=NULL,method=c("REML","ML"),interval=c(0,1000),...){
	corfunc<-corClass
	## preliminaries
	data<-data[tree$tip.label,]
	if(is.null(se)) se<-setNames(rep(0,Ntip(tree),
		tree$tip.label))
	## likelihood function
	lk<-function(sig2e,data,tree,model,ve){
		tree$edge.length<-tree$edge.length*sig2e
		ii<-sapply(1:Ntip(tree),function(x,e) which(e==x),
			e=tree$edge[,2])
		tree$edge.length[ii]<-tree$edge.length[ii]+
			ve[tree$tip.label]
		v<-diag(vcv(tree))
		vf<-varFixed(~v)
		COR<-corfunc(1,tree,...)
		fit<-gls(model,data=data,correlation=COR,weights=vf)
		-logLik(fit)
	}
	## estimate sig2[e]
	fit<-optimize(lk,interval=interval,
		data=data,tree=tree,model=model,ve=se^2)
	tree$edge.length<-tree$edge.length*fit$minimum
	ii<-sapply(1:Ntip(tree),function(x,e) which(e==x),
		e=tree$edge[,2])
	tree$edge.length[ii]<-tree$edge.length[ii]+
		se[tree$tip.label]^2
	v<-diag(vcv(tree))
	vf<-varFixed(~v)
	## fit & return model
	gls(model,data,correlation=corfunc(1,tree),weights=vf,
		method=method)
}