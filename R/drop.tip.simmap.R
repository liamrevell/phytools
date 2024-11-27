## functions drop or keep tips for multiple phylogenetic object types
## written by Liam J. Revell 2012, 2015, 2018, 2021, 2023, 2024

drop.tip.simmap<-function(phy,tip,...){
	if(length(setdiff(phy$tip.label,tip))==0){
		phy<-NULL
		warning("drop all tips of the tree: returning NULL")
	} else if(length(setdiff(phy$tip.label,tip))==1){
		keep<-setdiff(phy$tip.label,tip)
		nn<-which(phy$tip.label==keep)
		obj<-list(
			edge=matrix(c(2,1),1,2),
			edge.length=phy$edge.length[which(phy$edge[,2]==nn)],
			Nnode=1,
			tip.label=keep,
			mapped.edge=phy$mapped.edge[which(phy$edge[,2]==nn),,drop=FALSE],
			maps=phy$maps[which(phy$edge[,2]==nn)])
		rownames(obj$mapped.edge)<-"2,1"
		class(obj)<-class(phy)
		phy<-obj
	} else {
		if(hasArg(untangle)) untangle<-list(...)$untangle
		else untangle<-FALSE
		phy<-reorder(phy)
		if(!inherits(phy,"simmap")) 
			stop("phy should be object of class \"simmap\".")
		tip<-which(phy$tip.label%in%tip)
		edges<-match(tip,phy$edge[,2])
		z<-setdiff(1:nrow(phy$edge),edges)
		phy$edge<-phy$edge[z,]
		phy$edge.length<-phy$edge.length[z]
		phy$maps<-phy$maps[z]
		z<-setdiff(phy$edge[,2],phy$edge[,1])
		z<-z[z>Ntip(phy)]
		while(length(z)>0){
			edges<-match(z,phy$edge[,2])
			y<-setdiff(1:nrow(phy$edge),edges)
			phy$edge<-phy$edge[y,]
			phy$edge.length<-phy$edge.length[y]
			phy$maps<-phy$maps[y]
			z<-setdiff(phy$edge[,2],phy$edge[,1])
			z<-z[z>Ntip(phy)]
		}
		z<-setdiff(phy$edge[,2],phy$edge[,1])
		phy$tip.label<-phy$tip.label[z]
		phy$edge[which(phy$edge[,2]%in%z),2]<-1:Ntip(phy)
		while(sum(phy$edge[1,1]==phy$edge[,1])==1){
			phy$edge<-phy$edge[2:nrow(phy$edge),]
			phy$edge.length<-phy$edge.length[2:length(phy$edge.length)]
			phy$maps<-phy$maps[2:length(phy$maps)]
		}
		i<-1
		while(i<nrow(phy$edge)){
			single<-sum(phy$edge[i,2]==phy$edge[,1])==1
			while(single){
				z<-match(phy$edge[i,2],phy$edge[,1])
				phy$edge[i,2]<-phy$edge[z,2]
				phy$edge.length[i]<-phy$edge.length[i]+phy$edge.length[z]
				if(names(phy$maps[[i]])[length(phy$maps[[i]])]!=names(phy$maps[[z]])[1]){
					phy$maps[[i]]<-c(phy$maps[[i]],phy$maps[[z]])
				} else {
					phy$maps[[i]][length(phy$maps[[i]])]<-
						phy$maps[[i]][length(phy$maps[[i]])]+phy$maps[[z]][1]
					if(length(phy$maps[[z]])>1) phy$maps[[i]]<-
						c(phy$maps[[i]],phy$maps[[z]][2:length(phy$maps[[z]])])
				}
				y<-setdiff(1:nrow(phy$edge),z)
				phy$edge<-phy$edge[y,]
				phy$edge.length<-phy$edge.length[y]
				phy$maps<-phy$maps[y]
				single<-sum(phy$edge[i,2]==phy$edge[,1])==1
			}
			i<-i+1
		}
		z<-unique(as.vector(phy$edge))
		z<-z[z>Ntip(phy)]
		y<-order(z)+Ntip(phy)
		for(i in 1:nrow(phy$edge)) 
			for(j in 1:2) 
				if(phy$edge[i,j]%in%z) phy$edge[i,j]<-y[which(phy$edge[i,j]==z)]
		phy$Nnode<-max(phy$edge)-Ntip(phy)
		phy$node.states<-matrix(NA,nrow(phy$edge),2)
		for(i in 1:nrow(phy$edge)) phy$node.states[i,]<-
			c(names(phy$maps[[i]])[1],names(phy$maps[[i]])[length(phy$maps[[i]])])
		if(!is.null(phy$states)) phy$states<-phy$states[phy$tip.label]
		allstates<-vector()
		for(i in 1:nrow(phy$edge)) allstates<-c(allstates,names(phy$maps[[i]]))
		allstates<-unique(allstates)
		phy$mapped.edge<-matrix(data=0,length(phy$edge.length),length(allstates),
			dimnames=list(edge=apply(phy$edge,1,function(x) paste(x,collapse=",")),
			state=allstates))
		for(i in 1:length(phy$maps)) for(j in 1:length(phy$maps[[i]])) 
			phy$mapped.edge[i,names(phy$maps[[i]])[j]]<-phy$mapped.edge[i,
			names(phy$maps[[i]])[j]]+phy$maps[[i]][j]
		class(phy)<-c("simmap",setdiff(class(phy),"simmap"))
		if(untangle) phy<-untangle(phy,"read.tree")
	}
	phy
}

drop.tip.multiSimmap<-function(phy,tip,...){
	if(!inherits(phy,"multiSimmap"))
        stop("phy is not an object of class \"multiSimmap\".")
	trees<-lapply(phy,drop.tip.simmap,tip=tip,...)
	class(trees)<-class(phy)
	trees
}

keep.tip.simmap<-function(phy,tip,...){
	tips<-setdiff(phy$tip.label,tip)
	drop.tip.simmap(phy,tip=tips,...)
}

keep.tip.multiSimmap<-function(phy,tip,...){
	trees<-lapply(phy,keep.tip.simmap,tip=tip,...)
	class(trees)<-class(phy)
	trees
}
