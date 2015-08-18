## paint branch(es)
## written by Liam J. Revell 2014, 2015

paintBranches<-function(tree,edge,state,anc.state="1"){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.null(tree$maps)) maps<-lapply(tree$edge.length,function(x) setNames(x,anc.state))
	else maps<-tree$maps
	ii<-sapply(edge,function(x,y) which(y==x),y=tree$edge[,2])
	for(i in 1:length(ii)) maps[[ii[i]]]<-setNames(tree$edge.length[[ii[i]]],state)
	## build mapped.edge matrix
	s<-vector()
	for(i in 1:nrow(tree$edge)) s<-c(s,names(maps[[i]]))
	s<-unique(s)
	mapped.edge<-matrix(0,length(tree$edge.length),length(s),dimnames=list(edge=apply(tree$edge,1,function(x) paste(x,collapse=",")),state=s))
	for(i in 1:length(maps)) for(j in 1:length(maps[[i]])) mapped.edge[i,names(maps[[i]])[j]]<-mapped.edge[i,names(maps[[i]])[j]]+maps[[i]][j]
	## add attributes to the tree
	tree$mapped.edge<-mapped.edge
	tree$maps<-maps
	class(tree)<-c("simmap",setdiff(class(tree),"simmap"))
	tree
}

# function paints a subtree
# written by Liam Revell 2012, 2013, 2015

paintSubTree<-function(tree,node,state,anc.state="1",stem=FALSE){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(stem==0&&node<=length(tree$tip)) stop("stem must be TRUE for node<=N")
	if(is.null(tree$edge.length)) tree<-compute.brlen(tree)
	if(is.null(tree$maps)){
		maps<-as.list(tree$edge.length)
		for(i in 1:length(maps)) names(maps[[i]])<-anc.state
		
	} else maps<-tree$maps
	if(node>length(tree$tip)){
		desc<-getDescendants(tree,node)
		z<-which(tree$edge[,2]%in%desc)
		for(i in 1:length(z)){
			maps[[z[i]]]<-sum(maps[[z[i]]])
			names(maps[[z[i]]])<-state
		}
		if(stem){
			a<-which(tree$edge[,2]==node)
			b<-bySegments(maps[[a]])/max(bySegments(maps[[a]]))
			c<-match(1,((1-stem)>=b[,1])*((1-stem)<=b[,2]))
			d<-names(maps[[a]])[1:c]
			e<-maps[[a]]
			if(c>1){
				maps[[a]]<-c(e[1:(c-1)],(1-stem)*sum(e)-sum(e[1:(c-1)]),stem*sum(e))
				names(maps[[a]])<-c(d,state)
			} else {
				maps[[a]]<-c((1-stem)*sum(e),stem*sum(e))
				names(maps[[a]])<-c(d,state)
			}
		}
	} else {
		a<-which(tree$edge[,2]==node)
		b<-bySegments(maps[[a]])/max(bySegments(maps[[a]]))
		c<-match(1,((1-stem)>=b[,1])*((1-stem)<=b[,2]))
		d<-names(maps[[a]])[1:c]
		e<-maps[[a]]
		if(c>1){
			maps[[a]]<-c(e[1:(c-1)],(1-stem)*sum(e)-sum(e[1:(c-1)]),stem*sum(e))
			names(maps[[a]])<-c(d,state)
		} else {
			maps[[a]]<-c((1-stem)*sum(e),stem*sum(e))
			names(maps[[a]])<-c(d,state)
		}
	}
	s<-vector()
	for(i in 1:nrow(tree$edge)) s<-c(s,names(maps[[i]]))
	s<-unique(s)
	mapped.edge<-matrix(0,length(tree$edge.length),length(s),dimnames=list(edge=apply(tree$edge,1,function(x) paste(x,collapse=",")),state=s))
	for(i in 1:length(maps)) for(j in 1:length(maps[[i]])) mapped.edge[i,names(maps[[i]])[j]]<-mapped.edge[i,names(maps[[i]])[j]]+maps[[i]][j]
	tree$mapped.edge<-mapped.edge
	tree$maps<-maps
	class(tree)<-c("simmap",setdiff(class(tree),"simmap"))
	return(tree)
}

# function
# written by Liam Revell 2012

bySegments<-function(z){
	XX<-matrix(0,length(z),2,dimnames=list(names(z),c("start","end")))
	XX[1,2]<-z[1]
	if(length(z)>1){
		for(j in 2:length(z)){
			XX[j,1]<-XX[j-1,2]
			XX[j,2]<-XX[j,1]+z[j]
		}
	}
	return(XX)
}
