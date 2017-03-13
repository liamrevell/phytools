# simulates with multiple evolutionary rates in different parts of the tree
# written by Liam J. Revell 2011, 2013, 2015

sim.rates<-function(tree,sig2,anc=0,nsim=1,internal=F,plot=F){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.null(tree$mapped.edge)){
		message("tree does not contain a mapped discrete character history, using fastBM")
		X<-fastBM(tree,a=anc,sig2=sig2[1],nsim=nsim,internal=internal)
	} else {
		# first name (if necessary) and reorder sig2
		if(is.null(names(sig2))){
			message("names absent from sig2: assuming same order as $mapped.edge")
			if(length(sig2)==ncol(tree$mapped.edge)) names(sig2)<-colnames(tree$mapped.edge)
			else stop("the number of elements in sig2 should match the number of rows in mapped.edge")
		}
		sig2<-sig2[colnames(tree$mapped.edge)]
		# now create a tree for simulation
		edge.length<-rep(0,nrow(tree$edge))
		# scale the edge lengths by the rate
		for(i in 1:ncol(tree$mapped.edge))
			edge.length<-edge.length+sig2[i]*tree$mapped.edge[,i]
		names(edge.length)<-NULL
		tree<-list(Nnode=tree$Nnode,edge=tree$edge,tip.label=tree$tip.label,edge.length=edge.length)
		class(tree)<-"phylo"
		if(plot) plot(tree)
		# simulate
		X<-fastBM(tree,a=anc,nsim=nsim,internal=internal)
	}
	X
}

## function simulates multiple OU regimes using a difference equation approximation
## written by Liam J. Revell 2017

multiOU<-function(tree,alpha,sig2,theta=NULL,a0=NULL,nsim=1,internal=FALSE,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(nsim>1) x<-replicate(nsim,multiOU(tree,alpha,sig2,theta,a0,nsim=1,internal,...))
	else {
		rt<-Ntip(tree)+1
		if(!inherits(tree,"simmap")){
			tree<-paintSubTree(tree,rt,"1")
			names(alpha)<-names(sig2)<-"1"
		}
		if(hasArg(dt)) dt<-list(...)$dt
		else dt<-1/1000*max(nodeHeights(tree))
		ss<-sort(unique(c(getStates(tree,"tips"),getStates(tree,"nodes"))))
		if(is.null(theta)) theta<-setNames(rep(0,length(ss)),ss)
		if(is.null(a0)) a0<-0
		tree<-reorder(tree,"cladewise")
		S<-matrix(NA,nrow(tree$edge),2)
		S[which(tree$edge[,1]==rt),1]<-a0
		for(i in 1:nrow(tree$edge)){
			x1<-S[i,1]
			for(j in 1:length(tree$maps[[i]])){
				t<-tree$maps[[i]][j]
				ALPHA<-alpha[names(t)]
				SIG2<-sig2[names(t)]
				THETA<-theta[names(t)]
				t<-c(rep(dt,floor(t/dt)),t%%dt)
				for(k in 1:length(t))
					x1<-x1+ALPHA*(THETA-x1)*t[k]+
						rnorm(n=1,sd=sqrt(SIG2*t[k]))
			}
			S[i,2]<-x1
			if(any(tree$edge[,1]==tree$edge[i,2]))
				S[which(tree$edge[,1]==tree$edge[i,2]),1]<-S[i,2]
		}
		x<-setNames(c(S[1,1],S[,2]),c(tree$edge[1,1],
			tree$edge[,2]))[as.character(1:(Ntip(tree)+tree$Nnode))]
		names(x)[1:Ntip(tree)]<-tree$tip.label
		if(!internal) x<-x[tree$tip.label]
	}
	x
}
