## function drops tip or tips from an object of class "simmap"
## written by Liam J. Revell 2012, 2015, 2018

drop.tip.simmap<-function(tree,tip){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	tip<-which(tree$tip.label%in%tip)
	edges<-match(tip,tree$edge[,2])
	z<-setdiff(1:nrow(tree$edge),edges)
	tree$edge<-tree$edge[z,]
	tree$edge.length<-tree$edge.length[z]
	tree$maps<-tree$maps[z]
	z<-setdiff(tree$edge[,2],tree$edge[,1])
	z<-z[z>Ntip(tree)]
	while(length(z)>0){
		edges<-match(z,tree$edge[,2])
		y<-setdiff(1:nrow(tree$edge),edges)
		tree$edge<-tree$edge[y,]
		tree$edge.length<-tree$edge.length[y]
		tree$maps<-tree$maps[y]
		z<-setdiff(tree$edge[,2],tree$edge[,1])
		z<-z[z>Ntip(tree)]
	}
	z<-setdiff(tree$edge[,2],tree$edge[,1])
	tree$tip.label<-tree$tip.label[z]
	tree$edge[which(tree$edge[,2]%in%z),2]<-1:Ntip(tree)
	while(sum(tree$edge[1,1]==tree$edge[,1])==1){
		tree$edge<-tree$edge[2:nrow(tree$edge),]
		tree$edge.length<-tree$edge.length[2:length(tree$edge.length)]
		tree$maps<-tree$maps[2:length(tree$maps)]
	}
	i<-1
	while(i<nrow(tree$edge)){
		single<-sum(tree$edge[i,2]==tree$edge[,1])==1
		while(single){
			z<-match(tree$edge[i,2],tree$edge[,1])
			tree$edge[i,2]<-tree$edge[z,2]
			tree$edge.length[i]<-tree$edge.length[i]+tree$edge.length[z]
			if(names(tree$maps[[i]])[length(tree$maps[[i]])]!=names(tree$maps[[z]])[1]){
				tree$maps[[i]]<-c(tree$maps[[i]],tree$maps[[z]])
			} else {
				tree$maps[[i]][length(tree$maps[[i]])]<-tree$maps[[i]][length(tree$maps[[i]])]+tree$maps[[z]][1]
				if(length(tree$maps[[z]])>1) tree$maps[[i]]<-c(tree$maps[[i]],tree$maps[[z]][2:length(tree$maps[[z]])])
			}
			y<-setdiff(1:nrow(tree$edge),z)
			tree$edge<-tree$edge[y,]
			tree$edge.length<-tree$edge.length[y]
			tree$maps<-tree$maps[y]
			single<-sum(tree$edge[i,2]==tree$edge[,1])==1
		}
		i<-i+1
	}
	z<-unique(as.vector(tree$edge))
	z<-z[z>Ntip(tree)]
	y<-order(z)+Ntip(tree)
	for(i in 1:nrow(tree$edge)) for(j in 1:2) if(tree$edge[i,j]%in%z) tree$edge[i,j]<-y[which(tree$edge[i,j]==z)]
	tree$Nnode<-max(tree$edge)-Ntip(tree)
	tree$node.states<-matrix(NA,nrow(tree$edge),2)
	for(i in 1:nrow(tree$edge)) tree$node.states[i,]<-c(names(tree$maps[[i]])[1],names(tree$maps[[i]])[length(tree$maps[[i]])])
	if(!is.null(tree$states)) tree$states<-tree$states[tree$tip.label]
	allstates<-vector()
	for(i in 1:nrow(tree$edge)) allstates<-c(allstates,names(tree$maps[[i]]))
	allstates<-unique(allstates)
	tree$mapped.edge<-matrix(data=0,length(tree$edge.length),length(allstates),dimnames=list(edge=apply(tree$edge,1,function(x) paste(x,collapse=",")),state=allstates))
	for(i in 1:length(tree$maps)) for(j in 1:length(tree$maps[[i]])) tree$mapped.edge[i,names(tree$maps[[i]])[j]]<-tree$mapped.edge[i,names(tree$maps[[i]])[j]]+tree$maps[[i]][j]
	class(tree)<-c("simmap",setdiff(class(tree),"simmap"))
	return(tree)
}
