## function computes an estimate of the standing diversity in each category 
## given by x at each node
## written by Liam J. Revell 2011-2013, 2019

estDiversity<-function(tree,x,method=c("asr","simulation"),model="ER",...){
	if(!inherits(tree,"phylo")) 
		stop("tree should be object of class \"phylo\".")
	method<-matchType(method[1],c("asr","simulation"))
	if(hasArg(nsim)) nsim<-list(...)$nsim else nsim=100
	if(method=="asr"){
		tree<-reorder(tree,"cladewise") # reorder tree
		H<-nodeHeights(tree)
		if(!(model=="ER"||model=="SYM"||if(is.matrix(model)) isSymmetric(model) else FALSE)){
			cat("Warning:\n  only symmetric models allowed for method=\"asr\"\n")
			cat("            changing to model=\"ER\".\n")
			model<-"ER"
		}
		bb<-rerootingMethod(tree,x,model=model)
		aa<-bb$marginal.anc
		Q<-bb$Q
		D<-matrix(0,nrow(aa),ncol(aa),dimnames=dimnames(aa))
		# now loop through every node above the root
		message("Please wait. . . . Warning - this may take a while!")
		flush.console()
		for(i in 2:nrow(aa)){
			tt<-H[match(as.numeric(rownames(aa)[i]),tree$edge[,2]),2]
			ii<-H[,1]<tt&H[,2]>tt
			ee<-tree$edge[ii,2]
			hh<-H[ii,1]
			for(j in 1:length(ee)){
				tr<-reroot(tree,node.number=ee[j],position=tt-hh[j])
				D[i,]<-D[i,]+apeAce(tr,x[tr$tip.label],model=model,
					fixedQ=Q)$lik.anc[1,]
			}
			D[i,]<-D[i,]*aa[i,]
			if(i%%10==0){ 
				message(paste("Completed",i,"nodes"))
				flush.console()
			}
		}
		d<-rowSums(D)
	} else if(method=="simulation") {
		mtrees<-make.simmap(tree,x,nsim=nsim,model=model,message=FALSE)
		st<-sort(unique(x)); nn<-1:tree$Nnode+length(tree$tip)
		aa<-lapply(mtrees,function(x) describe.simmap(x,
			message=FALSE)$states)
		H<-nodeHeights(tree)
		D<-matrix(0,tree$Nnode,length(st),dimnames=list(nn,st))
		CC<-lapply(mtrees,function(x) lapply(x$maps,cumsum))
		# now loop through every node above the root
		message("Please wait. . . . Warning - this may take a while!")
		flush.console()
		for(i in 2:tree$Nnode){
			tt<-H[match(nn[i],tree$edge[,2]),2]
			ii<-H[,1]<tt&H[,2]>tt
			ee<-tree$edge[ii,2]
			hh<-H[ii,1]
			for(j in 1:length(ee)){
				dd<-setNames(rep(0,length(st)),st)
				for(k in 1:nsim){
					jj<-1; while((tt-hh[j])>CC[[k]][[which(tree$edge[,2]==
						ee[j])]][jj]) jj<-jj+1
					ss<-names(mtrees[[k]]$maps[[which(tree$edge[,2]==
						ee[j])]])[jj]
					dd[ss]<-dd[ss]+ if(ss==aa[[k]][i]) 1/nsim else 0
				}
				D[i,]<-D[i,]+dd
			}
			D[i,]<-D[i,]
			if(i%%10==0) message(paste("Completed",i,"nodes"))
		}
		d<-rowSums(D)
	}
	return(d)
}

