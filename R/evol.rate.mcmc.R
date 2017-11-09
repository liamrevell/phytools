## these functions uses a Bayesian MCMC approach to estimate heterogeneity in the evolutionary rate for a
## continuous character (Revell, Mahler, Peres-Neto, & Redelings. 2012.)
## code written by Liam J. Revell 2010, 2011, 2013, 2015, 2017

## function for Bayesian MCMC
## written by Liam J. Revell 2010, 2011, 2017
evol.rate.mcmc<-function(tree,x,ngen=10000,control=list(),...){	
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	# some minor error checking
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
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
	# first, try and obtain reasonable estimates for control parameters
	# and starting values for the MCMC
	C<-vcv.phylo(tree)
	C<-C[names(x),names(x)]
	n<-nrow(C)
	one<-matrix(1,n,1)
	a<-colSums(solve(C))%*%x/sum(solve(C)) # MLE ancestral value, used to start MCMC
	sig1<-as.numeric(t(x-one%*%a)%*%solve(C)%*%(x-one%*%a)/n) # MLE sigma-squared, used to start MCMC
	sig2<-sig1 # used to start MCMC
	flipped=FALSE # used to start MCMC
	# populate control list
	con=list(sig1=sig1,sig2=sig2,a=as.numeric(a),sd1=0.2*sig1,sd2=0.2*sig2,
		sda=0.2*abs(as.numeric(a)),kloc=0.2*mean(diag(C)),sdlnr=1,
		rand.shift=0.05,print=100,sample=100)
	# also might use: sig1mu=1000,sig2mu=1000
	con[(namc <- names(control))] <- control
	con<-con[!sapply(con,is.null)]
	# print control parameters to screen
	if(!quiet){
		message("Control parameters (set by user or default):")
		str(con)
	}
	# now detach the starting parameter values (to be compatible with downstream code)
	sig1<-con$sig1
	sig2<-con$sig1
	a<-con$a
	# all internal functions start here
	# function to return the index of a random edge
	random.node<-function(phy){
		# sum edges cumulatively
		cum.edge<-vector(mode="numeric")
		index<-vector(mode="numeric")
		for(i in 1:length(phy$edge.length)){
			if(i==1) cum.edge[i]<-phy$edge.length[1]
			else cum.edge[i]<-cum.edge[i-1]+phy$edge.length[i]
			index[i]<-phy$edge[i,2]
		}
		# pick random position
		pos<-runif(1)*cum.edge[length(phy$edge.length)]
		edge<-1
		while(pos>cum.edge[edge]) edge<-edge+1
		return (index[edge])
	}
	# return the indices of a vector that match a scalar
	match.all<-function(s,v){
		result<-vector(mode="numeric")
		j=1
		for(i in 1:length(v)){
			if(s==v[i]){
				result[j]=i
				j=j+1
			}
		}
		if(j==1) result<-NA
		return(result)
	}
	# function to return a matrix of the descendant tips from each internal & terminal node
	compute.descendant.species<-function(phy){
		D<-dist.nodes(phy)
		ntips<-length(phy$tip.label)
		Cii<-D[ntips+1,]
		C<-D
		C[,]<-0
		counts<-vector()
		for(i in 1:nrow(D)) for(j in 1:ncol(D)) C[i,j]<-(Cii[i]+Cii[j]-D[i,j])/2
		tol<-1e-10
		descendants<-matrix(0,nrow(D),ntips,dimnames=list(rownames(D)))
		for(i in 1:nrow(C)){
			k<-0
			for(j in 1:ntips){
				if(C[i,j]>=(C[i,i]-tol)){
					k<-k+1
					descendants[i,k]<-phy$tip.label[j]
				}
			}
			counts[i]<-k
		}
		names(counts)<-rownames(descendants)
		return(descendants=list(species=descendants,counts=counts))
	}
	# take a step on the tree (used in MCMC)
	tree.step<-function(phy,node,bp,step,up=NA,flip=FALSE){
		if(step<0)
			step=-step # if user has given -step, positivize
		if(is.na(up))
			up=(runif(1)>0.5) # if up/down is not assigned, we must be in the middle of a branch
		# decide to go up (1) or down (0) with equal probability
		if(up){
			# pick a new position along the branch (or go to the end)
			new.bp<-min(bp+step,phy$edge.length[match(node,phy$edge[,2])])
			# adjust step length to remain step
			step=step-(new.bp-bp)
			# check to see if we're done
			if(step<1e-6){
				return(list(node=node,bp=new.bp,flip=flip))
			} else {
				# we're going up, so get the daughters
				daughters<-phy$edge[match.all(node,phy$edge[,1]),2]
				# pick a random daughter
				new.node<-daughters[ceiling(runif(1)*length(daughters))]
				if(is.na(new.node)){
					location<-tree.step(phy,node,phy$edge.length[match(node,phy$edge[,2])],step,up=FALSE,flip) # we're at a tip
				} else {
					location<-tree.step(phy,new.node,0,step,up=TRUE,flip)
				}
			}
		} else {
			# pick a new position along the branch (or go to the start)
			new.bp<-max(bp-step,0)
			# adjust step length
			step=step-(bp-new.bp)
			# check to see if we're done
			if(step<1e-6){
				return(list(node=node,bp=new.bp,flip=flip))
			} else {
				# we're going down so find out who the parent is
				parent<-phy$edge[match(node,phy$edge[,2]),][1]
				# find the other daughter(s)				
				daughters<-phy$edge[match.all(parent,phy$edge[,1]),2]
				# don't use parent if root
				if(parent==(length(phy$tip.label)+1)){
					parent=NULL # if at the base of the tree
				}
				# create a vector of the possible nodes: the parent, and sister(s)				
				possible.nodes<-c(parent,daughters[-match(node,daughters)])
				# now pick randomly
				new.node<-possible.nodes[ceiling(runif(1)*length(possible.nodes))]
				# if parent
				if(is.null(parent)==FALSE&&new.node==parent){
					location<-tree.step(phy,new.node,phy$edge.length[match(new.node,phy$edge[,2])],step,up=FALSE,flip)
				} else {
					location<-tree.step(phy,new.node,0,step,up=TRUE,flip=TRUE)
				}
			}
		}
	}
	# log-likelihood function
	likelihood<-function(y,phy,C,descendants,sig1,sig2,a,loc){
		C1<-matrix(0,nrow(C),ncol(C),dimnames=dimnames(C))
		C2<-matrix(0,nrow(C),ncol(C),dimnames=dimnames(C))
		n<-length(y)
		D<-matrix(1,n,1)
		if(loc$node>length(phy$tip.label)){
			tr1<-extract.clade(phy,loc$node)
			tr1$root.edge<-phy$edge.length[match(loc$node,phy$edge[,2])]-loc$bp
			temp<-vcv.phylo(tr1)+tr1$root.edge
		} else {
			temp<-matrix(phy$edge.length[match(loc$node,phy$edge[,2])]-loc$bp,1,1,dimnames=list(c(phy$tip.label[loc$node]),c(phy$tip.label[loc$node])))
		}
		C2[rownames(temp),colnames(temp)]<-temp
		C1<-C-C2
		# tips<-phy$tip.label[-match(rownames(temp),phy$tip.label)]
		tips<-rownames(temp)		
		V<-sig1*C1+sig2*C2
		logL<-as.numeric(-t(y-D%*%a)%*%solve(V)%*%(y-D%*%a)/2-n*log(2*pi)/2-determinant(V)$modulus[1]/2)
		return(list(logL=logL,tips=tips))	
	}
	# prior probability function
	log.prior<-function(s1,s2,a,location){
		# logpr<-dexp(s1,rate=1/con$sig1mu,log=TRUE)+dexp(s2,rate=1/con$sig2mu,log=TRUE) # exponential prior
		logpr<-dnorm(log(s1)-log(s2),mean=0,sd=con$sdlnr,log=TRUE)-log(s1*s2) # log-normal
		return(logpr)
	}
	# proposal on sig1 or sig2
	propose.sig<-function(sig,scale){
		# if normal
		sig.prime<-abs(sig+rnorm(n=1,sd=scale)) # normal proposal distribution
		# if cauchy
		# sig.prime<-abs(sig+rcauchy(n=1,scale=scale))
		return(sig.prime)
	}
	# proposal on a
	propose.a<-function(a,scale){
		# if normal
		a.prime<-a+rnorm(n=1,sd=scale) # normal proposal distribution
		return(a.prime)
	}
	# proposal on loc
	propose.loc<-function(phy,loc,k,r){
		loc.prime<-list()
		loc.prime$flip=FALSE
		if(runif(1)>r){
			loc.prime<-tree.step(phy,loc$node,loc$bp,step=rexp(n=1,rate=1/k)) # update node & bp by random walk: rexp()
			# loc.prime<-tree.step(phy,loc$node,loc$bp,step=abs(rnorm(n=1,sd=sqrt(2)/k))) # update node & bp by random walk: rnorm()
		} else {
			loc.prime$node<-random.node(phy) # pick random branch
			loc.prime$bp<-runif(1)*phy$edge.length[match(loc.prime$node,phy$edge[,2])]
			if((runif(1)>0.5)) loc.prime$flip=TRUE
		}
		return(loc.prime)
	}
	# obtain remaining starting values for the MCMC
	location<-list()
	location$node<-random.node(tree)
	location$bp<-runif(1)*tree$edge.length[match(location$node,tree$edge[,2])]
	location$flip<-FALSE
	descendants<-compute.descendant.species(tree) # compute descendants
	logL<-likelihood(x,tree,C,descendants,sig1,sig2,a,location)$logL
	logpr<-log.prior(sig1,sig2,a,location)
	# create matrix for results
	results<-matrix(NA,floor(ngen/con$sample)+1,7,dimnames=list(c(0,1:(ngen/con$sample)),c("state","sig1","sig2","a","node","bp","likelihood")))
	curr.gen<-matrix(NA,1,7,dimnames=list("curr",c("state","sig1","sig2","a","node","bp","likelihood")))
	results[1,]<-c(0,sig1,sig2,a,location$node,location$bp,logL) # populate the first row
	curr.gen[1,]<-results[1,]
	group.tips<-list()
	tips<-list()
	group.tips[[1]]<-likelihood(x,tree,C,descendants,sig1,sig2,a,location)$tips
	tips[[1]]<-group.tips[[1]]
	message("Starting MCMC run....")	
	if(!quiet) print(results[1,])
	j<-2
	# now run Markov-chain
	for(i in 1:ngen){
		if(i%%4==1) sig1.prime<-propose.sig(sig1,scale=con$sd1) # update sig1
		else sig1.prime<-sig1
		if(i%%4==2) sig2.prime<-propose.sig(sig2,scale=con$sd2) # update sig2
		else sig2.prime<-sig2
		if(i%%4==3) a.prime<-propose.a(a,scale=con$sda) # update a
		else a.prime<-a
		if(i%%4==0){ 
			location.prime<-propose.loc(phy=tree,loc=location,k=con$kloc,r=con$rand.shift)
			if(location.prime$flip==TRUE){
				flipped.prime<-!flipped # flip the sigmas
			} else flipped.prime<-flipped
		} else { 
			location.prime<-location
			flipped.prime<-flipped
		}
		if(!flipped.prime){ 
			temp<-likelihood(x,tree,C,descendants,sig1.prime,sig2.prime,a.prime,location.prime)
			logpr.prime<-log.prior(sig1.prime,sig2.prime,a.prime,location.prime)
			if(exp(temp$logL+logpr.prime-curr.gen[1,"likelihood"]-logpr)>runif(1)){
				sig1<-sig1.prime
				sig2<-sig2.prime
				a<-a.prime
				location<-location.prime
 				logL<-temp$logL
				logpr<-logpr.prime
				flipped<-flipped.prime
				group.tips[[i+1]]<-temp$tips
			} else group.tips[[i+1]]<-group.tips[[i]]
		} else { 
			temp<-likelihood(x,tree,C,descendants,sig2.prime,sig1.prime,a.prime,location.prime)
			logpr.prime<-log.prior(sig2.prime,sig1.prime,a.prime,location.prime)
			if(exp(temp$logL+logpr.prime-curr.gen[1,"likelihood"]-logpr)>runif(1)){
				sig1<-sig1.prime
				sig2<-sig2.prime
				a<-a.prime
				location<-location.prime
 				logL<-temp$logL
				logpr<-logpr.prime
				flipped<-flipped.prime
				group.tips[[i+1]]<-setdiff(tree$tip.label,temp$tips)
			} else group.tips[[i+1]]<-group.tips[[i]]
		}	
		rm(temp)
		curr.gen[1,]<-c(i,sig1,sig2,a,location$node,location$bp,logL)
		if(i%%con$print==0)
			if(!quiet) print(curr.gen[1,])
		if(i%%con$sample==0){
			results[j,]<-curr.gen
			tips[[j]]<-group.tips[[i+1]]
			j<-j+1
		}
	} 
	message("Done MCMC run.")
	# return results
	obj<-list(mcmc=results,tips=tips,ngen=ngen,sample=con$sample)
	class(obj)<-"evol.rate.mcmc"
	obj
}

## S3 print method
print.evol.rate.mcmc<-function(x, ...){
	cat("\nObject of class \"evol.rate.mcmc\" containing the results from a\n")
	cat("the Bayesian MCMC analysis of a Brownian-motion rate-shift model.\n\n")
	cat(paste("MCMC was conducted for",x$ngen,"generations sampling every",x$sample,"\n"))
	cat("generations.\n\n")
	cat("The most commonly sampled rate shift(s) occurred on the edge(s)\n")
	pp<-table(x$mcmc[,"node"])/nrow(x$mcmc)
	node<-names(pp)[which(pp==max(pp))]
	if(length(node)==1) cat(paste("to node(s) ",node,".\n\n",sep=""))
	else cat(paste("leading to node(s) ",paste(node,collapse=", "),".\n\n",sep=""))
	cat("Use the functions posterior.evolrate and minSplit for more detailed\n")
	cat("analysis of the posterior sample from this analysis.\n\n")
}

# this function finds the split with the minimum the distance to all the other splits in the sample
# written by Liam J. Revell 2010, 2015
minSplit<-function(tree,split.list,method="sum",printD=FALSE){
	# some minor error checking
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	# determine is split.list is a list from our mcmc run, or a matrix
	if(class(split.list)=="list") split.list<-split.list$mcmc[,c("node","bp")]
	else if(class(split.list)=="matrix") split.list<-split.list[,c("node","bp")]
	# start by creating a matrix of the nodes of the tree with their descendant nodes
	D<-dist.nodes(tree)
	ntips<-length(tree$tip.label)
	Cii<-D[ntips+1,]
	C<-D
	C[,]<-0
	for(i in 1:nrow(D)) for(j in 1:ncol(D)) C[i,j]<-(Cii[i]+Cii[j]-D[i,j])/2
	tol<-1e-10
	descendants<-matrix(0,nrow(D),ncol(D),dimnames=list(rownames(D)))
	for(i in 1:nrow(C)){
		k<-1
		for(j in 1:ncol(C)){
			if(C[i,j]>=(C[i,i]-tol)){
				descendants[i,k]<-j
				k<-k+1
			}
		}
	}
	distances<-matrix(0,nrow(split.list),nrow(split.list))
	for(i in 1:nrow(split.list)){
		for(j in i:nrow(split.list)){
			if(i!=j){
				# first, if on the same branch as current average - then just compute the difference
				if(split.list[i,1]==split.list[j,1]){
					distances[i,j]<-abs(split.list[i,2]-split.list[j,2])
				} else {
					distances[i,j]<-D[split.list[i,1],split.list[j,1]]
					# is the split j downstream
					if(split.list[j,1]%in%descendants[split.list[i,1],]){
						# downstream
						distances[i,j]<-distances[i,j]-(tree$edge.length[match(split.list[j,1],tree$edge[,2])]-split.list[j,2])
						distances[i,j]<-distances[i,j]+(tree$edge.length[match(split.list[i,1],tree$edge[,2])]-split.list[i,2])
					} else if(split.list[i,1]%in%descendants[split.list[j,1],]){
						# ancestral
						distances[i,j]<-distances[i,j]-(tree$edge.length[match(split.list[i,1],tree$edge[,2])]-split.list[i,2])
						distances[i,j]<-distances[i,j]+(tree$edge.length[match(split.list[j,1],tree$edge[,2])]-split.list[j,2])
					} else {
						# neither
						distances[i,j]<-distances[i,j]-(tree$edge.length[match(split.list[i,1],tree$edge[,2])]-split.list[i,2])
						distances[i,j]<-distances[i,j]-(tree$edge.length[match(split.list[j,1],tree$edge[,2])]-split.list[j,2])
					}
				}
				distances[j,i]<-distances[i,j]					
			}
		}
	}
	if(method=="sumsq") distances<-distances^2
	if(method!="sumsq"&&method!="sum") message("allowable methods are 'sum' and 'sumsq' - using default method ('sum')")
	if(printD) print(distances)
	sum.dist<-colSums(distances)
	ind<-which.min(sum.dist) # this is the index of the minimum split
	# return minimum split
	return(list(node=split.list[ind,1],bp=split.list[ind,2]))
}

# function to analyze the posterior from evol.rate.mcmc()
# written by Liam Revell 2011
posterior.evolrate<-function(tree,ave.shift,mcmc,tips,showTree=FALSE){
	result<-matrix(NA,nrow(mcmc),7,dimnames=list(NULL,c("state","sig1","sig2","a","node","bp","likelihood")))
	tree$node.label<-NULL
	for(i in 1:nrow(mcmc)){
		shift=list(node=mcmc[i,"node"],bp=mcmc[i,"bp"])
		temp<-ave.rates(tree,shift,tips[[i]],mcmc[i,"sig1"],mcmc[i,"sig2"],ave.shift,showTree=showTree)
		result[i,]<-c(mcmc[i,"state"],temp[[1]],temp[[2]],mcmc[i,"a"],mcmc[i,"node"],mcmc[i,"bp"],mcmc[i,"likelihood"])
	}
	return(result)
}

# average the posterior rates
# written by Liam Revell 2011
ave.rates<-function(tree,shift,tips,sig1,sig2,ave.shift,showTree=TRUE){
	# first split and scale at shift
	unscaled<-splitTree(tree,shift)
	# now scale
	scaled<-unscaled
	if(length(setdiff(scaled[[1]]$tip.label,tips))!=1){
		scaled[[1]]$edge.length<-scaled[[1]]$edge.length*sig1
		scaled[[2]]$edge.length<-scaled[[2]]$edge.length*sig2
		if(!is.null(scaled[[2]]$root.edge)) scaled[[2]]$root.edge<-scaled[[2]]$root.edge*sig2
	} else {
		scaled[[1]]$edge.length<-scaled[[1]]$edge.length*sig2
		scaled[[2]]$edge.length<-scaled[[2]]$edge.length*sig1
		if(!is.null(scaled[[2]]$root.edge)) scaled[[2]]$root.edge<-scaled[[2]]$root.edge*sig1
	}
	# now bind
	tr.scaled<-paste.tree(scaled[[1]],scaled[[2]])
	if(showTree==TRUE) plot(tr.scaled)
	# now split tr.scaled and tree at ave.shift
	unscaled<-splitTree(tree,ave.shift)
	scaled<-splitTree(tr.scaled,ave.shift)
	# now compute the sig1 and sig2 to return
	sig1<-sum(scaled[[1]]$edge.length)/sum(unscaled[[1]]$edge.length)
	sig2<-sum(scaled[[2]]$edge.length)/sum(unscaled[[2]]$edge.length)
	return(list(sig1=sig1,sig2=sig2))
}




