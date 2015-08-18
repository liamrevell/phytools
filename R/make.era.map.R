# function creates a mapped tree in which the mappings are based on eras defined by "limits"
# written by Liam J. Revell 2011, 2013, 2015

make.era.map<-function(tree,limits,...){
	## set tolerance
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-5
	# check
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	H<-nodeHeights(tree) # compute node heights
	if(limits[1]>tol){ 
		warning("first value in limits should be zero.")
		limits[1]<-0
	}
	while(limits[length(limits)]>max(H)){
		warning("last value in limits should be less than the total tree height.")
		limits<-limits[1:(length(limits)-1)]
	}
	if(is.null(names(limits))) names(limits)<-1:length(limits)
	limits[length(limits)+1]<-max(H)
	maps<-list()
	# ok, now go through the edges of the tree
	for(i in 1:nrow(tree$edge)){
		s<-1
		while((H[i,1]>=(limits[s]-tol)&&H[i,1]<limits[s+1])==FALSE) s<-s+1
		e<-1
		while((H[i,2]>limits[e]&&H[i,2]<=(limits[e+1]+tol))==FALSE) e<-e+1
		maps[[i]]<-vector()
		if(s==e){
			maps[[i]][1]<-tree$edge.length[i]
			names(maps[[i]])[1]<-names(limits)[s]
		} else {
			maps[[i]][1]<-limits[s+1]-H[i,1]
			for(j in (s+1):e)
				maps[[i]][j-s+1]<-limits[j+1]-limits[j]
			maps[[i]][e-s+1]<-H[i,2]-limits[e]
			names(maps[[i]])<-names(limits)[s:e]
		}
	}
	tree$maps<-maps
	# now construct the matrix "mapped.edge" (for backward compatibility)
	allstates<-vector()
	for(i in 1:nrow(tree$edge)) allstates<-c(allstates,names(tree$maps[[i]]))
	allstates<-unique(allstates)
	tree$mapped.edge<-matrix(data=0,length(tree$edge.length),length(allstates),
		dimnames=list(apply(tree$edge,1,function(x) paste(x,collapse=",")),
		state=allstates))
	for(i in 1:length(tree$maps)) for(j in 1:length(tree$maps[[i]])) 
		tree$mapped.edge[i,names(tree$maps[[i]])[j]]<-tree$mapped.edge[i,
		names(tree$maps[[i]])[j]]+tree$maps[[i]][j]
	class(tree)<-c("simmap",setdiff(class(tree),"simmap"))
	tree
}
