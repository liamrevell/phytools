## function simulates stochastic character history under some model
## written by Liam J. Revell 2011, 2013, 2014, 2016, 2020

sim.history<-function(tree,Q,anc=NULL,nsim=1,direction=c("column_to_row","row_to_column"),
	...){
	if(!inherits(tree,"phylo")) 
		stop("tree should be an object of class \"phylo\".")
	if(hasArg(message)) message<-list(...)$message
	else message<-TRUE
	direction<-direction[1]
	direction<-strsplit(direction,"_")[[1]][1]
	# reorder to cladewise
	tree<-reorder.phylo(tree,"cladewise")
	# check Q
	if(!isSymmetric(Q)) if(message){
		if(direction=="column") 
			cat("Note - the rate of substitution from i->j should be given by Q[j,i].\n")
		else if(direction=="row")
			cat("Note - the rate of substitution from i->j should be given by Q[i,j].\n")
	}
	if(direction=="column"){
		if(!all(round(colSums(Q),10)==0)){
			if(all(round(rowSums(Q),10)==0)&&!isSymmetric(Q)){
				if(message){ 
					cat("Detecting that rows, not columns, of Q sum to zero :\n")
					cat("Transposing Q for internal calculations.\n")
				}
				Q<-t(Q)
			} else {
				if(message) 
					cat("Some columns (or rows) of Q don't sum to 0.0. Fixing.\n")
				diag(Q)<-0
				diag(Q)<--colSums(Q,na.rm=TRUE)
			}
		}
	} else if(direction=="row"){
		Q<-t(Q)
		if(!all(round(colSums(Q),10)==0)){
			if(all(round(rowSums(Q),10)==0)&&!isSymmetric(Q)){
				if(message){ 
					cat("Detecting that columns, not rows, of Q sum to zero :\n")
					cat("Transposing Q for internal calculations.\n")
				}
				Q<-t(Q)
			} else {
				if(message) 
					cat("Some columns (or rows) of Q don't sum to 0.0. Fixing.\n")
				diag(Q)<-0
				diag(Q)<--colSums(Q,na.rm=TRUE)
			}
		}
	}
	# does Q have names?
	if(is.null(dimnames(Q))) dimnames(Q)<-list(1:nrow(Q),1:ncol(Q))
	# create "multiPhylo" object
	mtrees<-vector(mode="list",length=nsim)
	class(mtrees)<-c("multiSimmap","multiPhylo")
	## deal with ancestral state
	if(is.null(anc)) 
		anc<-setNames(rep(1/ncol(Q),ncol(Q)),colnames(Q))
	if(is.character(anc)){ 
		anc<-colSums(to.matrix(anc,colnames(Q)))
		anc<-anc/sum(anc)
	}	
	# now loop
	for(i in 1:nsim){
		# set root state
		a<-rstate(anc)
		# create the map tree object
		mtree<-tree
		mtree$maps<-vector(mode="list",length=nrow(tree$edge))
		# now we want to simulate the node states on the tree
		node.states<-matrix(NA,nrow(tree$edge),ncol(tree$edge))
		node.states[which(tree$edge[,1]==(length(tree$tip)+1)),1]<-a
		for(j in 1:nrow(tree$edge)){
			if(tree$edge.length[j]==0){ 
				map<-vector()
				map[1]<-tree$edge.length[j]
				names(map)[1]<-
					node.states[which(tree$edge[,1]==tree$edge[j,2]),1]<-
					node.states[j,2]<-node.states[j,1]
			} else {
				time=0
				state<-node.states[j,1]
				new.state<-state
				dt<-vector()
				map<-vector()
				k<-1
				while(time<tree$edge.length[j]){
					dt[1]<-time
					dt[2]<-dt[1]+rexp(n=1,rate=-Q[state,state])
					if(dt[2]<tree$edge.length[j]) 
						new.state<-rstate(Q[,state][-match(state,rownames(Q))]/
							sum(Q[,state][-match(state,rownames(Q))]))
					dt[2]<-min(dt[2],tree$edge.length[j])
					map[k]<-dt[2]-dt[1]
					names(map)[k]<-state
					k<-k+1
					state<-new.state
					time<-dt[2]
				}
				names(map)[length(map)]->node.states[j,2]->
					node.states[which(tree$edge[,1]==tree$edge[j,2]),1]
			}
			mtree$maps[[j]]<-map
		}
		# add a couple of elements
		mtree$node.states<-node.states
		tip.states<-node.states[tree$edge[,2]<=length(tree$tip),2]
		tip.states<-tip.states[order(tree$edge[tree$edge[,2]<=length(tree$tip),2])]
		names(tip.states)<-tree$tip.label
		mtree$states<-tip.states
		# now construct the matrix "mapped.edge" (for backward compatibility
		allstates<-vector()
		for(j in 1:nrow(mtree$edge)) 
			allstates<-c(allstates,names(mtree$maps[[j]]))
		allstates<-unique(allstates)
		mtree$mapped.edge<-matrix(data=0,length(mtree$edge.length),
			length(allstates),dimnames=list(apply(mtree$edge,1,
			function(x) paste(x,collapse=",")),state=allstates))
		for(j in 1:length(mtree$maps)) 
			for(k in 1:length(mtree$maps[[j]])) 
				mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]<-
					mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]+
					mtree$maps[[j]][k]
		class(mtree)<-c("simmap",setdiff(class(mtree),"simmap"))
		mtrees[[i]]<-mtree	
	}
	if(nsim==1) mtrees<-mtrees[[1]]
	if(message) cat("Done simulation(s).\n")
	mtrees
}

## simulate DNA sequence from a tree & model parameters
## written by Liam J. Revell 2013, 2019

genSeq<-function(tree,l=1000,Q=NULL,rate=1,format="DNAbin",...){
	if(is.null(Q)){ 
		Q<-matrix(1,4,4)
		rownames(Q)<-colnames(Q)<-c("a","c","g","t")
		diag(Q)<-0
		diag(Q)<--colSums(Q)
	}
	if(length(rate)!=l){
		if(length(rate)==1) rate<-rep(rate,l)
		else {
			cat("warning: length(rate) & l should match for length(rate)>1\n")
			cat("         rate will be recycled.\n")
			rate<-rep(rate,ceiling(l/length(rate)))[1:l]
		}
	}
	cat("simulating sequences....\n")
	flush.console()
	X<-sapply(rate,function(a,b,c) sim.Mk(b,a*c),b=tree,c=Q)
	if(format=="DNAbin") return(as.DNAbin(X))
	else if(format=="phyDat") return(as.phyDat(X))
	else if(format=="matrix") return(X)
}

