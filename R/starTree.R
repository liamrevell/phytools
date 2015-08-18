# Function creates a polytomous tree
# written by Liam Revell 2011

starTree<-function(species,branch.lengths=NULL){
	n<-length(species)
	edge<-matrix(NA,n,2)
	if(!is.null(branch.lengths)) edge.length=branch.lengths
	edge[,1]<-n+1; edge[,2]<-1:n
	if(!is.null(branch.lengths)) phy<-list(edge=edge,edge.length=edge.length,Nnode=1,tip.label=as.vector(species))
	else phy<-list(edge=edge,Nnode=1,tip.label=as.vector(species))
	class(phy)<-"phylo"
	return(phy)
}
