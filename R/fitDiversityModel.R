# this function fits a "diversity-dependent-evolutionary-diversification" model (similar to Mahler et al. 2010)
# written by Liam Revell, 2010/2011/2012

fitDiversityModel<-function(tree,x,d=NULL,showTree=TRUE,tol=1e-6){
	# some minor error checking
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)) x<-x[,1]
	if(is.null(names(x))){
		if(length(x)==length(tree$tip)){
			message("x has no names; assuming x is in the same order as tree$tip.label")
			names(x)<-tree$tip.label
		} else
			stop("x has no names and is a different length than tree$tip.label")
	}
	if(any(is.na(match(tree$tip.label,names(x))))){
		message("some species in tree are missing from data, dropping missing taxa from the tree")
		tree<-drop.tip(tree,tree$tip.label[-match(names(x),tree$tip.label)])
	}
	if(any(is.na(match(names(x),tree$tip.label)))){
		message("some species in data are missing from tree, dropping missing taxa from the data")
		x<-x[tree$tip.label]
	}
	if(any(is.na(x))){
		message("some data given as 'NA', dropping corresponding species from tree")
		tree<-drop.tip(tree,names(which(is.na(x))))
	}
	if(!is.null(d)){
		if(is.data.frame(d)) d<-as.matrix(d)
		if(is.matrix(d)) d<-d[,1]
		if(is.null(names(d))){
			if(length(d)==tree$Nnode){
				message("d has no names; assuming d is in node number order of the resolved tree")
				names(d)<-c(length(tree$tip)+tol:tree$Nnode)
			} else
				stop("d has no names and is a different length than tree$Nnode for the resolved tree")
		}
	} else {
		message("no values for lineage density provided; computing assuming single biogeographic region")
		# compute lineage diversity at each node
		ages<-branching.times(tree)
		d<-vector()
		for(i in 1:length(ages)) d[i]<-sum(ages>ages[i])
		names(d)<-names(ages)
	}
	maxd<-max(d)
	d<-d/(maxd+tol)
	# likelihood function
	lik<-function(theta,y,phy,diversity){
		scaled.psi<-theta
		for(i in 1:nrow(phy$edge)){
			vi<-phy$edge.length[i]
			phy$edge.length[i]<-vi+vi*scaled.psi*diversity[as.character(phy$edge[i,1])]
		}
		D<-vcv(phy)
		D<-D[names(y),names(y)]
		Dinv<-solve(D)
		a<-as.numeric(colSums(Dinv)%*%y/sum(Dinv))
		sig0<-as.numeric(t(y-a)%*%Dinv%*%(y-a)/nrow(D))
		Dinv<-Dinv/sig0; D<-D*sig0
		logL<-as.numeric(-t(y-a)%*%Dinv%*%(y-a)/2-determinant(D)$modulus[1]/2-length(y)*log(2*pi)/2)
		if(showTree) plot(phy)
		return(logL)
	}
	# optimize
	res<-optimize(lik,c(-1,1),y=x,phy=tree,diversity=d,maximum=TRUE)
	# compute condition sig0
	compute_sig0<-function(scaled.psi,y,phy,diversity){
		for(i in 1:nrow(phy$edge)){
			vi<-phy$edge.length[i]
			phy$edge.length[i]<-vi+vi*scaled.psi*diversity[as.character(phy$edge[i,1])]
		}
		D<-vcv(phy)
		D<-D[names(y),names(y)]
		Dinv<-solve(D)
		a<-as.numeric(colSums(Dinv)%*%y/sum(Dinv))
		sig0<-as.numeric(t(y-a)%*%Dinv%*%(y-a)/nrow(D))
		return(sig0)
	}
	sig0=compute_sig0(res$maximum,x,tree,d)
	# compute the Hessian
	compute_Hessian<-function(scaled.psi,sig0,y,phy,d,maxd){
		psi<-scaled.psi*sig0/(maxd+tol)
		likHessian<-function(theta,y,phy,d,maxd){
			sig0<-theta[1]
			psi<-theta[2]
			for(i in 1:nrow(phy$edge)){
				vi<-phy$edge.length[i]
				phy$edge.length[i]<-vi*(sig0+psi*d[as.character(phy$edge[i,1])]*(maxd+tol))
			}
			D<-vcv(phy)
			D<-D[names(y),names(y)]
			Dinv<-solve(D)
			a<-as.numeric(colSums(Dinv)%*%y/sum(Dinv))
			logL<-as.numeric(-t(y-a)%*%Dinv%*%(y-a)/2-determinant(D)$modulus[1]/2-length(y)*log(2*pi)/2)
			return(logL)
		}
		H<-hessian(likHessian,c(sig0,psi),y=y,phy=phy,d=d,maxd=maxd)
		return(H)
	}
	H<-compute_Hessian(res$maximum,sig0,x,tree,d,maxd)
	# return results to user
	if(var(d)>0)
		return(list(logL=res$objective,sig0=sig0,psi=sig0*res$maximum/(maxd+tol),vcv=matrix(solve(-H),2,2,dimnames=list(c("sig0","psi"),c("sig0","psi")))))
	else {
		message("psi not estimable because diversity is constant through time.")
		return(list(logL=res$objective,sig0=sig0,vcv=matrix(-1/H[1,1],1,1,dimnames=list(c("sig0"),c("sig0")))))
	}
}
