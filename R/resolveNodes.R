## function to resolve nodes in a multifurcating tree
## written by Liam J. Revell 2016

resolveNode<-function(tree,node){
	dd<-Children(tree,node)
	if(length(dd)>2){
		EL<-!is.null(tree$edge.length)
		if(!EL) tree<-compute.brlen(tree)
		n<-length(dd)
		tt<-lapply(allTrees(n,TRUE,dd),untangle,"read.tree")
		ROOT<-node==(Ntip(tree)+1)
		SPL<-if(!ROOT) splitTree(tree,split=list(node=node,
			bp=tree$edge.length[which(tree$edge[,2]==node)])) else
			list(NULL,tree)
		KIDS<-Children(SPL[[2]],SPL[[2]]$edge[1,1])
		KIDS<-setNames(KIDS,dd)[KIDS>Ntip(SPL[[2]])]
		SUBS<-list()
		if(length(KIDS)>0)
			for(i in 1:length(KIDS)) 
				SUBS[[i]]<-extract.clade(SPL[[2]],KIDS[i])
		obj<-list()
		for(i in 1:length(tt)){
			tt[[i]]$edge.length<-rep(0,nrow(tt[[i]]$edge))
			for(j in 1:Ntip(tt[[i]]))
				tt[[i]]$edge.length[which(tt[[i]]$edge[,2]==j)]<-
					tree$edge.length[which(tree$edge[,2]==
					as.numeric(tt[[i]]$tip.label[j]))]
			ind<-as.numeric(tt[[i]]$tip.label)<=Ntip(tree)
			tt[[i]]$tip.label[ind]<-
				tree$tip.label[as.numeric(tt[[i]]$tip.label[ind])]
			if(length(KIDS)>0)
				for(j in 1:length(KIDS))
					tt[[i]]<-bind.tree(tt[[i]],SUBS[[j]],
						where=which(tt[[i]]$tip.label==
						names(KIDS)[j]))	
			obj[[i]]<-if(!ROOT) bind.tree(SPL[[1]],tt[[i]],
				where=which(SPL[[1]]$tip.label=="NA")) else
				tt[[i]]
			if(!EL) obj[[i]]$edge.length<-NULL
		}
		class(obj)<-"multiPhylo"
	} else obj<-tree
	obj
}

## function to resolve all nodes in a tree with multifurcations
## written by Liam J. Revell 2016

resolveAllNodes<-function(tree){
	foo<-function(node,tree) length(Children(tree,node))
	nodes<-1:tree$Nnode+Ntip(tree) ## all nodes
	nodes<-nodes[sapply(1:tree$Nnode+Ntip(tree),foo,
		tree=tree)>2]
	for(i in 1:length(nodes)){
		if(i==1) obj<-resolveNode(tree,nodes[i])
		else {
			for(j in 1:length(obj)){
				MM<-matchNodes(tree,obj[[j]])
				NODE<-MM[which(MM[,1]==nodes[i]),2]
				if(j==1) tmp<-resolveNode(obj[[j]],NODE)
				else tmp<-c(tmp,resolveNode(obj[[j]],NODE))
			}
			obj<-tmp
		}
	}
	obj
}


