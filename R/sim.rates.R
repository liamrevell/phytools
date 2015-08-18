# simulates with multiple evolutionary rates in different parts of the tree
# written by Liam J. Revell 2011, 2013, 2015

sim.rates<-function(tree,sig2,anc=0,nsim=1,internal=F,plot=F){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.null(tree$mapped.edge)){
		message("tree does not contain a mapped discrete character history, using fastBM")
		X<-fastBM(tree,a=anc,sig2=sig2[1],nsim=nsim,internal=internal)
	} else {
		# first name (if necessary) and reorder sig2
		if(is.null(names(sig2))){
			message("names absent from sig2: assuming same order as $mapped.edge")
			if(length(sig2)==ncol(tree$mapped.edge)) names(sig2)<-colnames(tree$mapped.edge)
			else stop("the number of elements in sig2 should match the number of rows in mapped.edge")
		}
		sig2<-sig2[colnames(tree$mapped.edge)]
		# now create a tree for simulation
		edge.length<-rep(0,nrow(tree$edge))
		# scale the edge lengths by the rate
		for(i in 1:ncol(tree$mapped.edge))
			edge.length<-edge.length+sig2[i]*tree$mapped.edge[,i]
		names(edge.length)<-NULL
		tree<-list(Nnode=tree$Nnode,edge=tree$edge,tip.label=tree$tip.label,edge.length=edge.length)
		class(tree)<-"phylo"
		if(plot) plot(tree)
		# simulate
		X<-fastBM(tree,a=anc,nsim=nsim,internal=internal)
	}
	X
}
