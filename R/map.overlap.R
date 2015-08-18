# function computes the fractional (i.e., 0->1) overlap between the stochastic maps of two trees
# written by Liam Revell 2011, 2015

map.overlap<-function(tree1,tree2,tol=1e-6){
	if(!inherits(tree1,"phylo")||!inherits(tree2,"phylo")) 
		stop("both input trees should be objects of class \"phylo\".")
	if(!all.equal.phylo(tree1,tree2,tolerance=tol)) 
		stop("mapped trees must have the same underlying structure")
	if(is.null(tree1$maps)||is.null(tree2$maps)) 
		stop("tree should be a simmap modified \"phylo\" object")
	overlap<-0
	for(i in 1:nrow(tree1$edge)){
		XX<-matrix(0,length(tree1$maps[[i]]),2,
			dimnames=list(names(tree1$maps[[i]]),c("start","end")))
		XX[1,2]<-tree1$maps[[i]][1]
		if(length(tree1$maps[[i]])>1){
			for(j in 2:length(tree1$maps[[i]])){
				XX[j,1]<-XX[j-1,2]
				XX[j,2]<-XX[j,1]+tree1$maps[[i]][j]
			}
		}
		YY<-matrix(0,length(tree2$maps[[i]]),2,
			dimnames=list(names(tree2$maps[[i]]),c("start","end")))
		YY[1,2]<-tree2$maps[[i]][1]
		if(length(tree2$maps[[i]])>1){		
			for(j in 2:length(tree2$maps[[i]])){
				YY[j,1]<-YY[j-1,2]
				YY[j,2]<-YY[j,1]+tree2$maps[[i]][j]
			}
		}
		for(j in 1:nrow(XX)){
			lower<-which(YY[,1]<=XX[j,1]); lower<-lower[length(lower)]
			upper<-which(YY[,2]>=(XX[j,2]-tol))[1]
			for(k in lower:upper){
				if(rownames(XX)[j]==rownames(YY)[k])
					overlap<-overlap+min(YY[k,2],XX[j,2])-
						max(YY[k,1],XX[j,1])
			}
		}
	}
	return(overlap/sum(tree1$edge.length))
}
