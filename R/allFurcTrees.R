# function returns all unrooted multi & bifurcating tree topologies
# written by Liam Revell 2011, 2013

allFurcTrees<-function(n,tip.label=NULL,to.plot=TRUE){
	# check to see if tip labels have been provided
	if(is.null(tip.label)) tip.label=as.character(1:n)
	# now pick three species at random to start
	new<-list(stree(n=3,tip.label=sample(tip.label,3)))
	class(new)<-"multiPhylo"
	added<-new[[1]]$tip.label; remaining<-setdiff(tip.label,added)
	# loop
	while(length(remaining)>0){
		old<-new; new<-list()
		new.tip<-sample(remaining,1)
		for(i in 1:length(old)){			
			temp<-add.to.branches.and.nodes(old[[i]],new.tip)
			new<-unlist(list(new,temp),recursive=FALSE); class(new)<-"multiPhylo"
		}
		added<-c(added,new.tip)
		remaining<-setdiff(tip.label,added)
	}
	if(to.plot) new<-read.tree(text=write.tree(new))
	return(new) # return all trees
}

# function adds a tip to all branches and nodes
# written by Liam Revell 2011

add.to.branches.and.nodes<-function(tree,tip.name){
	n.edge<-nrow(tree$edge)
	tree$edge.length<-rep(1,n.edge) # set all edge lengths to 1.0
	# create new tip
	new.tip<-list(edge=matrix(c(2L,1L),1,2),tip.label=tip.name,edge.length=1,Nnode=1L); class(new.tip)<-"phylo"
	# first add the tip to all edges
	trees<-list(); class(trees)<-"multiPhylo"
	for(i in 1:n.edge){
		trees[[i]]<-bind.tree(tree,new.tip,where=tree$edge[i,2],position=0.5) # add tip to edge
		trees[[i]]$edge.length<-NULL
	}
	# ok, now add the tip to all internal nodes
	intnodes<-unique(tree$edge[,1])
	for(i in 1:tree$Nnode){
		trees[[i+n.edge]]<-bind.tree(tree,new.tip,where=intnodes[i]) # add tip to node
		trees[[i+n.edge]]$edge.length<-NULL
	}	
	return(trees)
}

