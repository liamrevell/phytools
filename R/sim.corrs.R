# sim.corrs
# written by Liam J. Revell 2012, 2013, 2015

sim.corrs<-function(tree,vcv,anc=NULL,internal=FALSE){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(!is.list(vcv)){
		p<-nrow(vcv)
		if(is.null(anc)) anc<-rep(0,p)
		cholvcv<-chol(vcv)
		X<-matrix(rnorm(p*nrow(tree$edge),sd=rep(sqrt(tree$edge.length),p)),nrow(tree$edge),p)
		X<-X%*%cholvcv
	} else {
		p<-nrow(vcv[[1]])
		if(is.null(anc)) anc<-rep(0,p)
		if(is.null(names(vcv))){
			names(vcv)<-colnames(tree$mapped.edge)
			message("names absent from vcv: assuming same order as $mapped.edge")
		}
		vcv<-vcv[colnames(tree$mapped.edge)]
		cholvcv<-lapply(vcv,chol)
		X<-matrix(0,nrow(tree$edge),p)
		for(i in 1:length(vcv)){
			Y<-matrix(rnorm(p*nrow(tree$edge),sd=rep(sqrt(tree$mapped.edge[,i]),p)),nrow(tree$edge),p)
			X<-X+Y%*%cholvcv[[i]]
		}
	}
	Y<-array(0,dim=c(nrow(tree$edge),ncol(tree$edge),p))
	n<-length(tree$tip)
	for(i in 1:nrow(X)){
		if(tree$edge[i,1]==(n+1))
			Y[i,1,]<-anc
		else
			Y[i,1,]<-Y[match(tree$edge[i,1],tree$edge[,2]),2,]
		Y[i,2,]<-Y[i,1,]+X[i,]
	}
	X<-matrix(data=rbind(Y[1,1,],as.matrix(Y[,2,])),length(tree$edge.length)+1,p)
	rownames(X)<-c(n+1,tree$edge[,2])
	X<-as.matrix(X[as.character(1:(n+tree$Nnode)),])
	rownames(X)[1:n]<-tree$tip.label
	if(internal==TRUE)
		return(X)
	else
		return(X[1:length(tree$tip.label),])

}
