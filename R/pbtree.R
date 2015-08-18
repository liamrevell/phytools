# function to simulate a pure-birth phylogenetic tree or trees
# written by Liam J. Revell 2011-2015

pbtree<-function(b=1,d=0,n=NULL,t=NULL,scale=NULL,nsim=1,type=c("continuous","discrete"),...){
	# get arguments
	if(hasArg(ape)) ape<-list(...)$ape
	else ape<-TRUE
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(hasArg(extant.only)) extant.only<-list(...)$extant.only
	else extant.only<-FALSE
	if(hasArg(max.count)) max.count<-list(...)$max.count
	else max.count<-1e5
	if(hasArg(method)) method<-list(...)$method
	else method<-"rejection"
	if(hasArg(tip.label)){
		tip.label<-list(...)$tip.label
		if(!is.null(tip.label)){
			if(is.null(n)){ 
				tip.label<-NULL
				cat("Warning: tip.label not allowed for n=NULL.\n")
				cat("         using default labels\n")
			} else if(length(tip.label)==n){ 
				if(d>0){ 
					cat("Warning: only using labels in tip.label for extant tips.\n")
					cat("         extinct tips will be labeled X1, X2, etc.\n")
				}
			} else if(length(tip.label)!=n) {
				cat("Warning: length(tip.label) and n do not match.\n")
				cat("         using default labels\n")
				tip.label<-NULL
			}
		}
	} else tip.label<-NULL
	type<-matchType(type[1],c("continuous","discrete"))
	if(type=="discrete"){
		if((b+d)>1){ 
			cat("Warning:\n  b + d cannot exceed 1.0 in discrete-time simulations\n")
			cat(paste("  setting b & d to",b/(b+d),"and",d/(b+d),"respectively\n"))
			b<-b/(b+d)
			d<-d/(b+d)
		}
	}
	tol<-1e-12
	# done get arguments
	# if nsim > 1 replicate nsim times
	if(nsim>1){
		trees<-replicate(nsim,pbtree(b,d,n,t,scale,type=type,ape=ape,quiet=quiet,extant.only=extant.only,method=method,tip.label=tip.label),simplify=FALSE)
		class(trees)<-"multiPhylo"
		return(trees)
	} else {
		if(!is.null(n)) NN<-n else NN<-NULL
		if(!is.null(n)&&!is.null(t)){
			if(method=="rejection"){
				# simulate taxa & time stop using rejection sampling to max.count
				if(!quiet){
					cat("simulating with both taxa-stop (n) and time-stop (t) is\n")
					cat("performed via rejection sampling & may be slow\n\n")
				}
				N<-0; T<-0; count<--1
				while((N!=n||T<t)&&count<max.count){ 
					tree<-pbtree(b,d,t=t,type=type,ape=FALSE,extant.only=extant.only,quiet=TRUE)
					if(is.null(tree)) N<-T<-0
					else {
						N<-if(d>0&&extant.only==FALSE) length(getExtant(tree)) else length(tree$tip.label)
						T<-max(nodeHeights(tree))
					}
					count<-count+1
				} 
				if(N==n&&T>=(t-tol)){ 
					if(!quiet) cat(paste("  ",count," trees rejected before finding a tree\n\n",sep=""))
				} else { 
					if(!quiet) cat(paste("  max count of ",count," reached without finding a tree\n\n",sep=""))
					tree<-NULL
				}
				# done simulate taxa & time stop
			} else if(method=="direct"){
				# simulate using direct sampling (experimental)
				if(!quiet){ 
					cat("simulating with both taxa-stop (n) & time-stop (t) using\n")
					cat("'direct' sampling. this is experimental\n")
				}
				m<-2
				while(m[length(m)]!=n){
					ll<-bd<-vector(); m<-2; i<-1
					while(sum(ll)<t){
						ll[i]<-rexp(n=1,rate=m[i]*(b+d))
						if(sum(ll)<t) bd[i]<-2*rbinom(n=1,size=1,prob=b/(b+d))-1 else bd[i]<-0
						m[i+1]<-m[i]+bd[i]; i<-i+1
						if(m[i]==0) break
					}
				}
				ll[length(ll)]<-ll[length(ll)]+t-sum(ll)
				bd<-bd
				node<-1; dead<-1
				edge<-matrix(c(node,NA,node,NA),2,2,byrow=T)
				edge.length<-c(0,0)
				node<-node+1
				for(i in 1:length(bd)){
					o<-is.na(edge[,2])
					p<-which(o)
					l<-ll[i]
					birth<-if(bd[i]==1) TRUE else FALSE
					q<-if(length(p)>1) sample(p,size=1) else p
					if(birth){
						# new edge
						edge[q,2]<-node
						edge<-rbind(edge,matrix(c(node,NA,node,NA),2,2,byrow=T))
						node<-node+1
					} else {
						edge[q,2]<--dead
						dead<-dead+1
					}	
					edge.length[p]<-edge.length[p]+l
					if(birth) edge.length<-c(edge.length,rep(0,2))
				}
				edge[edge[,2]<0,2]<-NA
				o<-is.na(edge[,2])
				n<-sum(o)
				edge<-edge+n
				p<-which(o)
				edge[o,2]<-1:sum(is.na(edge[,2]))	
				# build 'phylo' object
				tree<-list(edge=edge,edge.length=edge.length,tip.label=paste("t",1:n,sep=""),Nnode=n-1)
				class(tree)<-"phylo"
				# done simulate using direct sampling (experimental)
			}	
		} else {
			if(!is.null(t)){
				# simulation time stop
				node<-1; dead<-1
				edge<-matrix(c(node,NA,node,NA),2,2,byrow=T)
				edge.length<-c(0,0)
				node<-node+1; tt<-0
				while(tt<t){
					o<-is.na(edge[,2])
					if(!any(o)) break
					p<-which(o)
					if(type=="discrete"){
						l<-rgeom(n=sum(o),prob=b+d)+1
						l<-l[which(l==min(l))]
					} else l<-rexp(n=1,sum(o)*(b+d))
					tt<-tt+l[1]
					if(tt>=t) l<-l-tt+t
					else {
						birth<-sapply(runif(n=length(l)),function(x) if(x<b/(b+d)) TRUE else FALSE)
						q<-if(length(p)>1) sample(p)[1:min(length(p),length(l))] else p
						for(i in 1:length(l)){
							if(birth[i]){
								# new edge
								edge[q[i],2]<-node
								edge<-rbind(edge,matrix(c(node,NA,node,NA),2,2,byrow=T))
								node<-node+1
							} else {
								edge[q[i],2]<--dead
								dead<-dead+1
							}
						}
						
					}
					edge.length[p]<-edge.length[p]+l[1]
					edge.length<-c(edge.length,rep(0,2*length(l)))
				}
				edge[edge[,2]<0,2]<-NA
				o<-is.na(edge[,2])
				n<-sum(o)
				edge<-edge+n
				p<-which(o)
				edge[o,2]<-1:sum(is.na(edge[,2]))
				# done unique part of time stop
			} else if(!is.null(n)) {
				# simulate taxa stop
				node<-1
				edge<-matrix(c(node,NA,node,NA),2,2,byrow=T)
				edge.length<-c(0,0)
				node<-node+1; dead<-1; nn<-2
				while(nn<n){
					o<-is.na(edge[,2])
					if(!any(o)) break
					p<-which(o)
					if(type=="discrete"){
						l<-rgeom(n=sum(o),prob=b+d)+1
						l<-l[which(l==min(l))]
					} else l<-rexp(n=1,sum(o)*(b+d))
					birth<-birth<-sapply(runif(n=length(l)),function(x) if(x<b/(b+d)) TRUE else FALSE)
					q<-if(length(p)>1) sample(p)[1:min(length(p),length(l))] else p
					for(i in 1:length(l)){
						if(birth[i]){
							# new edge
							edge[q[i],2]<-node
							edge<-rbind(edge,matrix(c(node,NA,node,NA),2,2,byrow=T))
							node<-node+1
						} else {
							edge[q[i],2]<--dead
							dead<-dead+1
						}
					}
					edge.length[p]<-edge.length[p]+l[1]
					edge.length<-c(edge.length,rep(0,2*length(l)))
					nn<-sum(is.na(edge[,2]))
				}
				edge[edge[,2]<0,2]<-NA
				o<-is.na(edge[,2])
				nn<-sum(o)
				edge<-edge+nn
				p<-which(o)
				l<-if(type=="discrete") min(rgeom(n=sum(o),prob=(b+d))+1) else rexp(n=1,sum(o)*(b+d))			
				edge.length[p]<-edge.length[p]+l		
				edge[is.na(edge[,2]),2]<-1:sum(is.na(edge[,2]))
				if((nn-dead+1)>n&&!quiet){ 
					# this might happen in discrete time only
					cat("Warning:\n  due to multiple speciation events in the final time interval\n")
					cat("  realized n may not equal input n\n\n")
					if(!is.null(tip.label)){
						cat("Warning: length(tip.label) and n do not match.\n")
						cat("         using default labels\n")
						tip.label<-NULL
					}
				}
				n<-nn
				# done unique part of taxa stop
			}
			# build 'phylo' object with temporary labels
			tree<-list(edge=edge,edge.length=edge.length,tip.label=1:n,Nnode=n-1)
			class(tree)<-"phylo"
			if(!is.null(scale)){
				# rescale if scale!=NULL
				h<-max(nodeHeights(tree))
				tree$edge.length<-scale*tree$edge.length/h
			}
			if(d>0&&extant.only){
				# prune extinct tips if extant.only==TRUE
				if(length(getExtinct(tree))==(length(tree$tip.label)-1)){ 
					if(!quiet) cat("Warning:\n  no extant tips, tree returned as NULL\n")
					tree<-NULL
				} else tree<-drop.tip(tree,getExtinct(tree))
			} 
			# if tree!=NULL assign final tip labels
			if(!is.null(tree)) tree$tip.label<-paste("t",1:length(tree$tip.label),sep="")
			
		}
		# if ape==TRUE make sure 'phylo' is consistent with ape
		if(ape&&is.null(tree)==FALSE) tree<-read.tree(text=write.tree(tree))
		if(!is.null(tip.label)){
			th<-max(nodeHeights(tree))
			if(length(getExtant(tree,tol=1e-08*th))!=NN){
				# simulation must have gone extint before reaching NN
				tree$tip.label<-paste("X",1:length(tree$tip.label),sep="")
			} else {
				ll<-getExtant(tree,tol=1e-08*th)
				ii<-sapply(ll,function(x,y) which(x==y),y=tree$tip.label)
				tree$tip.label[ii]<-tip.label
				tree$tip.label[-ii]<-paste("X",1:(length(tree$tip.label)-length(ii)),sep="")
			}
		}
		# done
		return(tree)
	}
}
