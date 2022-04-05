## function for Matrix Representation Parsimony supertree estimation in R
## uses pratchet() or optim.parsimony() from the "phangorn" package
## written by Liam J. Revell 2011, 2013, 2015, 2022

compute.mr<-function(trees,type=c("phyDat","matrix")){
	type<-type[1]
	if(inherits(trees,"phylo")){ 
		trees<-list(trees)
		class(trees)<-"multiPhylo"
	}
	if(!inherits(trees,"multiPhylo")) stop("trees should be an object of class \"phylo\" or \"multiPhylo\".")
	# compute matrix representation phylogenies
	X<-list() # list of bipartitions
	characters<-0 # number of characters
	for(i in 1:length(trees)){
		temp<-prop.part(trees[[i]]) # find all bipartitions
		# create matrix representation of trees[[i]] in X[[i]]
		X[[i]]<-matrix(0,nrow=length(trees[[i]]$tip),ncol=length(temp)-1)
		for(j in 1:ncol(X[[i]])) X[[i]][c(temp[[j+1]]),j]<-1
		rownames(X[[i]])<-attr(temp,"labels") # label rows
		if(i==1) species<-trees[[i]]$tip.label
		else species<-union(species,trees[[i]]$tip.label) # accumulate labels
		characters<-characters+ncol(X[[i]]) # count characters
	}
	XX<-matrix(data="?",nrow=length(species),ncol=characters,dimnames=list(species))
	j<-1
	for(i in 1:length(X)){
		# copy each of X into supermatrix XX
		XX[rownames(X[[i]]),c(j:((j-1)+ncol(X[[i]])))]<-X[[i]][1:nrow(X[[i]]),1:ncol(X[[i]])]
		j<-j+ncol(X[[i]])
	}
	if(type=="phyDat"){
		# compute contrast matrix for phangorn
		contrast<-matrix(data=c(1,0,0,1,1,1),3,2,dimnames=list(c("0","1","?"),c("0","1")),byrow=TRUE)
		# convert XX to phyDat object
	 	XX<-phyDat(XX,type="USER",contrast=contrast)
	}
	XX
}

mrp.supertree<-function(trees,method=c("pratchet","optim.parsimony"),...){
	# set method
	method<-method[1]
	# some minor error checking
	if(!inherits(trees,"multiPhylo")) stop("trees should be an object of class \"multiPhylo\".")
	XX<-compute.mr(trees,type="phyDat")
	# estimate supertree
	if(method=="pratchet"){
		if(hasArg(start)){
			start<-list(...)$start
			if(inherits(start,"phylo")){
				supertree<-pratchet(XX,all=TRUE,...)
			} else {
				if(start=="NJ") start<-NJ(dist.hamming(XX))
				else if(start=="random") start<-rtree(n=length(XX),tip.label=names(XX))
				else {
					warning("do not recognize that option for start; using random starting tree")
					tree<-rtree(n=length(XX),tip.label=names(XX))
				}
				args<-list(...)
				args$start<-start
				args$data<-XX
				args$all<-TRUE
				supertree<-do.call(pratchet,args)
			}
		} else supertree<-pratchet(XX,all=TRUE,...)
		if(inherits(supertree,"phylo"))
			message(paste("The MRP supertree, optimized via pratchet(),\nhas a parsimony score of ",
				attr(supertree,"pscore")," (minimum ",attr(XX,"nr"),")",sep=""))
		else if(inherits(supertree,"multiPhylo"))
			message(paste("pratchet() found ",length(supertree)," supertrees\nwith a parsimony score of ",
				attr(supertree[[1]],"pscore")," (minimum ",attr(XX,"nr"),")",sep=""))
	} else if(method=="optim.parsimony"){
		if(hasArg(start)){
			start<-list(...)$start
			if(inherits(start,"phylo")){
				supertree<-optim.parsimony(tree=start,data=XX,...)
			} else {
				if(start=="NJ") start<-NJ(dist.hamming(XX))
				else if(start=="random") start<-rtree(n=length(XX),tip.label=names(XX))
				else {
					warning("do not recognize that option for tree; using random starting tree")
					start<-rtree(n=length(XX),tip.label=names(XX))
				}
				supertree<-optim.parsimony(tree=start,data=XX,...)
			}			
		} else {
			message("no input starting tree or option for optim.parsimony; using random starting tree")
			start<-rtree(n=length(XX),tip.label=names(XX))
			supertree<-optim.parsimony(tree=start,data=XX,...)
		}
		if(inherits(supertree,"phylo"))
			message(paste("The MRP supertree, optimized via optim.parsimony(),\nhas a parsimony score of ",
				attr(supertree,"pscore")," (minimum ",attr(XX,"nr"),")",sep=""))
		else if(inherits(supertree,"multiPhylo"))
			message(paste("optim.parsimony() found ",length(supertree)," supertrees\nwith a parsimony score of ",
				attr(supertree[[1]],"pscore")," (minimum ",attr(XX,"nr"),")",sep=""))
	}
	return(supertree)
}
