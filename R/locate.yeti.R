## code to place a missing extant taxon into a tree using ML or REML on continuous data
## written by Liam J. Revell 2014

locate.yeti<-function(tree,X,...){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	if(hasArg(method)) method<-list(...)$method
	else method<-"ML"
	if(hasArg(search)) search<-list(...)$search
	else search<-"heuristic"
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-FALSE
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(hasArg(rotate)) rotate<-list(...)$rotate
	else rotate<-if(method=="ML") TRUE else FALSE
	root.node<-length(tree$tip.label)+1
	if(hasArg(constraint)){
		if(search=="exhaustive") constraint<-list(...)$constraint
		else {
			cat("constraint only works with search==\"exhaustive\"\n")
			constraint<-c(root.node,tree$edge[,2])
		}
	} else constraint<-c(root.node,tree$edge[,2])
	if(!is.matrix(X)) X<-as.matrix(X)
	tip<-setdiff(rownames(X),tree$tip.label)
	if(method=="ML") mltree<-yetiML(tree,X,quiet,tip,root.node,constraint,plot,search,rotate)
	else if(method=="REML") mltree<-yetiREML(tree,X,quiet,tip,root.node,constraint,plot,search)
	else { 
		cat(paste("Do not recognize method ",method,".\n",sep=""))
		stop()
	}
	mltree
}

yetiML<-function(tree,X,quiet,tip,root.node,constraint,plot,search,rotate){
	if(!quiet) cat(paste("Optimizing the phylogenetic position of ",tip," using ML. Please wait....\n",sep=""))
	if(ncol(X)>1&&rotate){
		pca<-phyl.pca(tree,X[tree$tip.label,])
		obj<-phyl.vcv(X[tree$tip.label,],vcv(tree),1)
		X<-(X-matrix(rep(obj$a[,1],nrow(X)),nrow(X),ncol(X),byrow=TRUE))%*%pca$Evec
	}
	if(search=="heuristic"){
		trees<-list()
		ee<-c(root.node,tree$edge[,2])
		for(i in 1:length(ee)) trees[[i]]<-bind.tip(tree,tip,where=ee[i],position=if(ee[i]==root.node) 0 else 0.5*tree$edge.length[i-1])
		class(trees)<-"multiPhylo"
		lik.edge<-function(tree,XX,rotate){
			if(!rotate) XX<-phyl.pca(tree,XX[tree$tip.label,])$S
			obj<-phyl.vcv(as.matrix(XX[tree$tip.label,]),vcv(tree),1)
			ll<-vector()
			for(i in 1:ncol(XX)) ll[i]<-sum(dmnorm(XX[tree$tip.label,i],mean=rep(obj$a[i,1],nrow(XX)),obj$C*obj$R[i,i],log=TRUE))
			sum(ll)
		}
		logL<-sapply(trees,lik.edge,XX=X,rotate=rotate)
		if(plot){
			ll<-logL[2:length(logL)]
			ll[ll<=sort(ll,decreasing=TRUE)[ceiling(nrow(tree$edge)/2)]]<-sort(ll,decreasing=TRUE)[ceiling(nrow(tree$edge)/2)]
			layout(matrix(c(1,2),2,1),heights=c(0.95,0.05))
			plotBranchbyTrait(tree,ll,mode="edges",title="log(L)",show.tip.label=FALSE)
			edgelabels(round(logL[2:length(logL)],1),cex=0.5)
			plot.new()
			text(paste("Note: logL <=",round(min(ll),2),"set to",round(min(ll),2),"for visualization only"),x=0.5,y=0.5)
		}
		edge<-ee[which(logL==max(logL))]
	}
	lik.tree<-function(position,tip,tree,edge,XX,rt,rotate){
		if(edge==rt) tree<-bind.tip(tree,tip,edge.length=position,where=edge)
		else tree<-bind.tip(tree,tip,where=edge,position=position)
		if(!rotate) XX<-phyl.pca(tree,XX[tree$tip.label,])$S
		obj<-phyl.vcv(as.matrix(XX[tree$tip.label,]),vcv(tree),1)
		ll<-vector()
		for(i in 1:ncol(XX)) ll[i]<-sum(dmnorm(XX[tree$tip.label,i],mean=rep(obj$a[i,1],nrow(XX)),obj$C*obj$R[i,i],log=TRUE))
		sum(ll)
	}
	if(search=="heuristic"){
		ee<-edge
		if(edge!=root.node) ee<-c(ee,getAncestors(tree,node=edge,type="parent"))
		if(edge>length(tree$tip.label)) ee<-c(ee,tree$edge[which(tree$edge[,1]==edge),2])
	} else if(search=="exhaustive") ee<-c(root.node,tree$edge[,2])
	ee<-intersect(ee,constraint)
	fit<-vector(mode="list",length=length(ee))
	for(i in 1:length(ee)){
		if(ee[i]==root.node) fit[[i]]<-optimize(lik.tree,interval=c(max(nodeHeights(tree)),10*max(nodeHeights(tree))),tip=tip,tree=tree,
			edge=ee[i],XX=X,rt=root.node,rotate=rotate,maximum=TRUE)
		else fit[[i]]<-optimize(lik.tree,interval=c(0,tree$edge.length[which(tree$edge[,2]==ee[i])]),tip=tip,tree=tree,edge=ee[i],
			XX=X,rt=root.node,rotate=rotate,maximum=TRUE)
	}
	logL<-sapply(fit,function(x) x$objective)
	if(search=="exhaustive"&&plot){
		ll<-sapply(fit,function(x) x$objective)
		ll<-ll[2:length(ll)]
		ll[ll<=sort(ll,decreasing=TRUE)[ceiling(nrow(tree$edge)/2)]]<-sort(ll,decreasing=TRUE)[ceiling(nrow(tree$edge)/2)]
		layout(matrix(c(1,2),2,1),heights=c(0.95,0.05))
		plotBranchbyTrait(tree,ll,mode="edges",title="log(L)",show.tip.label=FALSE)
		edgelabels(round(ll,1),cex=0.5)
		plot.new()
		text(paste("Note: logL <=",round(min(ll),2),"set to",round(min(ll),2),"for visualization only"),x=0.5,y=0.5)
	}
	fit<-fit[[which(logL==max(logL))]]
	edge<-ee[which(logL==max(logL))]
	mltree<-if(edge==root.node) midpoint.root(bind.tip(tree,tip,where=edge,edge.length=fit$maximum)) else bind.tip(tree,tip,where=edge,position=fit$maximum)
	mltree$logL<-fit$objective
	if(!quiet) cat("Done.\n")
	mltree
}

yetiREML<-function(tree,X,quiet,tip,root.node,constraint,plot,search){
	if(!quiet){
		cat("---------------------------------------------------------------\n")
		cat("| **Warning: method=\"REML\" has not been thoroughly tested.    |\n")
		cat("|   Use with caution.**                                       |\n")
		cat("---------------------------------------------------------------\n\n")
	}
	if(!quiet) cat(paste("Optimizing the phylogenetic position of ",tip," using REML. Please wait....\n",sep=""))
	if(search=="heuristic"){
		trees<-list()
		ee<-c(root.node,tree$edge[,2])
		for(i in 1:length(ee)) trees[[i]]<-bind.tip(tree,tip,where=ee[i],position=if(ee[i]==root.node) 0 else 0.5*tree$edge.length[i-1])
		class(trees)<-"multiPhylo"
		lik.edge<-function(tree,XX){
			tree<-multi2di(tree)
			YY<-apply(XX[tree$tip.label,],2,pic,phy=tree)
			vcv<-t(YY)%*%YY/nrow(YY)
			E<-eigen(vcv)$vectors
			##a<-apply(XX,2,function(x,tree) ace(x,tree,type="continuous",method="pic")$ace[1],tree=tree)
			##S<-(X-matrix(rep(a,nrow(X)),nrow(X),ncol(X),byrow=TRUE))%*%E
			##ZZ<-apply(S,2,pic,phy=tree)
			ZZ<-YY%*%E
			vars<-diag(t(ZZ)%*%ZZ/nrow(ZZ))
			ll<-vector()
			for(i in 1:ncol(ZZ)) ll[i]<-sum(dnorm(ZZ[,i],mean=0,sd=sqrt(vars[i]),log=TRUE))
			sum(ll)
		}
		logL<-sapply(trees,lik.edge,XX=X)
		if(plot){
			ll<-logL[2:length(logL)]
			ll[ll<=sort(ll,decreasing=TRUE)[ceiling(nrow(tree$edge)/2)]]<-sort(ll,decreasing=TRUE)[ceiling(nrow(tree$edge)/2)]
			layout(matrix(c(1,2),2,1),heights=c(0.95,0.05))
			plotBranchbyTrait(tree,ll,mode="edges",title="log(L)",show.tip.label=FALSE)
			edgelabels(round(logL[2:length(logL)],1),cex=0.5)
			plot.new()
			text(paste("Note: logL <=",round(min(ll),2),"set to",round(min(ll),2),"for visualization only"),x=0.5,y=0.5)
		}
		edge<-ee[which(logL==max(logL))]
	}
	lik.tree<-function(position,tip,tree,edge,XX,rt){
		if(edge==rt) tree<-bind.tip(tree,tip,edge.length=position,where=edge)
		else tree<-bind.tip(tree,tip,where=edge,position=position)
		tree<-multi2di(tree)
		YY<-apply(XX,2,pic,phy=tree,scaled=FALSE)
		vcv<-t(YY)%*%YY/nrow(YY)
		sum(dmnorm(YY,mean=rep(0,ncol(YY)),varcov=vcv,log=TRUE))
	}
	if(search=="heuristic"){
		ee<-edge
		if(edge!=root.node) ee<-c(ee,getAncestors(tree,node=edge,type="parent"))
		if(edge>length(tree$tip.label)) ee<-c(ee,tree$edge[which(tree$edge[,1]==edge),2])
	} else if(search=="exhaustive") ee<-c(root.node,tree$edge[,2])
	ee<-intersect(ee,constraint)
	fit<-vector(mode="list",length=length(ee))
	for(i in 1:length(ee)){
		if(ee[i]==root.node) fit[[i]]<-optimize(lik.tree,interval=c(max(nodeHeights(tree)),10*max(nodeHeights(tree))),tip=tip,tree=tree,edge=ee[i],XX=X,rt=root.node,maximum=TRUE)
		else fit[[i]]<-optimize(lik.tree,interval=c(0,tree$edge.length[which(tree$edge[,2]==ee[i])]),tip=tip,tree=tree,edge=ee[i],XX=X,rt=root.node,maximum=TRUE)
	}
	logL<-sapply(fit,function(x) x$objective)
	if(search=="exhaustive"&&plot){
		ll<-sapply(fit,function(x) x$objective)
		ll<-ll[2:length(ll)]
		ll[ll<=sort(ll,decreasing=TRUE)[ceiling(nrow(tree$edge)/2)]]<-sort(ll,decreasing=TRUE)[ceiling(nrow(tree$edge)/2)]
		layout(matrix(c(1,2),2,1),heights=c(0.95,0.05))
		plotBranchbyTrait(tree,ll,mode="edges",title="log(L)",show.tip.label=FALSE)
		edgelabels(round(ll,1),cex=0.5)
		plot.new()
		text(paste("Note: logL <=",round(min(ll),2),"set to",round(min(ll),2),"for visualization only"),x=0.5,y=0.5)
	}
	fit<-fit[[which(logL==max(logL))]]
	edge<-ee[which(logL==max(logL))]
	mltree<-if(edge==root.node) midpoint.root(bind.tip(tree,tip,where=edge,edge.length=fit$maximum)) else bind.tip(tree,tip,where=edge,position=fit$maximum)
	mltree$logL<-fit$objective
	if(!quiet) cat("Done.\n")
	mltree
}

