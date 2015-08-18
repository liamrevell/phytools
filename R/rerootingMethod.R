# function to compute the marginal posterior probabilities for nodes using the rerooting method
# written by Liam J. Revell 2013, 2015

rerootingMethod<-function(tree,x,model=c("ER","SYM"),...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(!is.binary.tree(tree)){
		mt<-tree
		tree<-multi2di(tree)
		multif<-TRUE
	} else multif<-FALSE
	if(hasArg(tips)) tips<-list(...)$tips
	else tips<-NULL
	if(!is.matrix(model)) model<-model[1]
	n<-length(tree$tip.label)
	# if vector convert to binary matrix
	if(!is.matrix(x)){ 
		yy<-to.matrix(x,sort(unique(x)))
		if(is.null(tips)) tips<-FALSE
	} else { 
		if(is.null(tips)) tips<-TRUE
		yy<-x
	}
	yy<-yy[tree$tip.label,]
	yy<-yy/rowSums(yy)
	YY<-apeAce(tree,yy,model=model)
	Q<-matrix(c(0,YY$rates)[YY$index.matrix+1],length(YY$states),length(YY$states),
		dimnames=list(YY$states,YY$states))
	diag(Q)<--colSums(Q,na.rm=TRUE)
	nn<-if(tips) c(1:n,2:tree$Nnode+n) else 2:tree$Nnode+n
	ff<-function(nn){
		tt<-reroot(tree,node.number=nn,position=tree$edge.length[which(tree$edge[,2]==nn)])
		apeAce(tt,yy,model=model,fixedQ=Q)$lik.anc[1,]
	}
	XX<-t(sapply(nn,ff))
	if(tips) XX<-rbind(XX[1:n,],YY$lik.anc[1,],XX[(n+1):nrow(XX),])
	else XX<-rbind(yy,YY$lik.anc[1,],XX)
	rownames(XX)<-1:(tree$Nnode+n)
	if(tips) rownames(XX)[1:n]<-tree$tip.label
	XX<-if(tips) XX else XX[1:tree$Nnode+n,]
	if(multif){
		M<-matchNodes(mt,tree)
		ii<-sapply(M[,2],function(x,y) which(y==x),y=as.numeric(rownames(XX)))
		XX<-XX[ii,]
		rownames(XX)<-M[,1]		
	}
	return(list(loglik=YY$loglik,Q=Q,marginal.anc=XX))
}
