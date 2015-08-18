# function computes the Strahler number at each node
# written by Liam J. Revell 2013, 2015

strahlerNumber<-function(tree,plot=TRUE){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	pw<-reorder(tree,"pruningwise")
	oo<-sapply(tree$edge[,2],function(x,y) which(x==y),y=pw$edge[,2])
	SS<-matrix(0,nrow(pw$edge),2)
	SS[pw$edge[,2]<=length(pw$tip.label),2]<-1
	nn<-unique(pw$edge[,1])
	for(i in 1:pw$Nnode){
		jj<-which(pw$edge[,1]==nn[i])
		s<-sort(SS[jj,2],decreasing=TRUE)
		SS[jj,1]<-if(all(sapply(s[2:length(s)],"<",s[1]))) s[1] else s[1]+1
		SS[which(pw$edge[,2]==pw$edge[jj[1],1]),2]<-SS[jj[1],1]
	}
	ss<-setNames(c(SS[oo,][1,1],SS[oo,2]),c(tree$edge[1,1],tree$edge[,2]))
	ss<-ss[order(as.numeric(names(ss)))]
	names(ss)[1:length(tree$tip.label)]<-tree$tip.label
	if(plot){ 
		plotTree(tree)
		nodelabels(ss[1:tree$Nnode+length(tree$tip.label)])
	}
	return(ss)
}

# extracts all the most inclusive clades with Strahler number i from tree
# written by Liam J. Revell 2013, 2015

extract.strahlerNumber<-function(tree,i,plot=TRUE){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	sn<-strahlerNumber(tree)
	sn<-sn[sn==i]
	# get descendant tip numbers for all clades
	ll<-lapply(as.numeric(names(sn)),getDescendants,tree=tree)
	# figure out which ones are most inclusive
	ff<-function(x,y) !all(sapply(x,"%in%",y))
	GG<-sapply(ll,function(x,y) sapply(y,ff,x=x),y=ll)
	ii<-which(colSums(GG)==(ncol(GG)-1))
	# extract these clades
	trees<-lapply(as.numeric(names(sn))[ii],extract.clade,phy=tree)
	if(plot){
		nplots<-2*ceiling(length(trees)/2)
		layout(matrix(1:nplots,ceiling(nplots/2),min(c(length(trees),2)),byrow=TRUE))
		sNN<-lapply(trees,strahlerNumber,plot=TRUE)
	}
	if(length(trees)>1) class(trees)<-"multiPhylo" else trees<-trees[[1]]
	return(trees)
}
