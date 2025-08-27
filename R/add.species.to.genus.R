## function identifies the genus from a species name & attempts to add it to the tree
## written by Liam J. Revell 2013

add.species.to.genus<-function(tree,species,genus=NULL,where=c("root","random")){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(!is.ultrametric(tree))
		warning("this code has only been tested with ultrametric tree\n  your tree may be returned without edge lengths")
	where<-where[1]
	if(is.null(genus)){
		## get genus from species name
		x<-strsplit(species,"")[[1]]
		i<-1
		while(x[i]!="_"&&x[i]!=" ") i<-i+1
		genus<-paste(x[2:i-1],collapse="")
	}
	ii<-grep(paste(genus,"_",sep=""),tree$tip.label)
	if(length(ii)>1){
		if(!is.monophyletic(tree,tree$tip.label[ii]))
			warning(paste(genus,"may not be monophyletic\n  attaching to the most inclusive group containing members of this genus"))
		nn<-findMRCA(tree,tree$tip.label[ii])
		if(where=="root") tree<-bind.tip(tree,gsub(" ","_",species),where=nn)
		else if(where=="random"){
			tt<-splitTree(tree,list(node=nn,bp=tree$edge.length[which(tree$edge[,2]==nn)]))
			tt[[2]]<-add.random(tt[[2]],tips=gsub(" ","_",species))
			tree<-paste.tree(tt[[1]],tt[[2]])
		} else stop("option 'where' not recognized")
	} else if(length(ii)==1){
		nn<-ii
		if(where=="root") 
			tree<-bind.tip(tree,gsub(" ","_",species),where=nn,position=0.5*tree$edge.length[which(tree$edge[,2]==nn)])
		else if(where=="random")
			tree<-bind.tip(tree,gsub(" ","_",species),where=nn,position=runif(n=1)*tree$edge.length[which(tree$edge[,2]==nn)])
		else
			stop("option 'where' not recognized")
	} else
		warning("could not match your species to a genus\n  check spelling, including case")
	tree
}

## function take genus backbone tree & converts genus tree to species tree by simulating pure-birth subtrees
## written by Liam J. Revell 2015

genus.to.species.tree<-function(tree,species, method = c("random", "star")){
	N<-Ntip(tree)
	genera<-tree$tip.label
	species<-gsub(" ","_",species)
	method <- method[1]
	for(i in 1:N){
		jj<-grep(paste(genera[i],"_",sep=""),species)
		nn<-which(tree$tip.label==genera[i])
		if(length(jj)>1){
			if (method == "random"){
      				h <- runif(n = 1) * tree$edge.length[which(tree$edge[,2] == nn)]
      				tree$edge.length[which(tree$edge[, 2] == nn)] <- tree$edge.length[which(tree$edge[, 2] == nn)] - h
      				sub.tree <- pbtree(n = length(jj), scale = h, tip.label = species[jj])
      				} else {
        			sub.tree <- stree(length(jj), tip.label = species[jj])
      			}
			tree<-bind.tree(tree,sub.tree,where=nn)
		} else if(length(jj)==1) tree$tip.label[nn]<-species[jj]
	}
	tree
}
