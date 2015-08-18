# function does fast estimation of ML ancestral states using ace
# written by Liam J. Revell 2012, 2013, 2015

fastAnc<-function(tree,x,vars=FALSE,CI=FALSE){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	if(!is.binary.tree(tree)) btree<-multi2di(tree)
	else btree<-tree
	M<-btree$Nnode
	N<-length(btree$tip.label)
	anc<-v<-vector()
	for(i in 1:M+N){
		a<-multi2di(root(btree,node=i))
   		anc[i-N]<-ace(x,a,method="pic")$ace[1]
   		names(anc)[i-N]<-i
		if(vars||CI){
			picx<-pic(x,a,rescaled.tree=TRUE)
			b<-picx$rescaled.tree
			d<-which(b$edge[,1]==(length(b$tip.label)+1))
			v[i-N]<-(1/b$edge.length[d[1]]+1/b$edge.length[d[2]])^(-1)*mean(picx$contr^2)
			names(v)[i-N]<-names(anc)[i-N]
		}
 	}
	if(!is.binary.tree(tree)){
		ancNames<-matchNodes(tree,btree)
		anc<-anc[as.character(ancNames[,2])]
		names(anc)<-ancNames[,1]
		if(vars||CI){ 
			v<-v[as.character(ancNames[,2])]
			names(v)<-ancNames[,1]
		}
	}
	result<-list(ace=anc)
	if(vars) result$var<-v
	if(CI){ 
		result$CI95<-cbind(anc-1.96*sqrt(v),anc+1.96*sqrt(v))
		rownames(result$CI95)<-names(anc)
	}
	if(length(result)==1) return(result$ace)
	else return(result)
}
