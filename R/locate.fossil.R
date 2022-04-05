## code to place a fossil taxon into a tree using ML and continuous data
## written by Liam J. Revell 2014, 2015, 2017, 2022

locate.fossil<-function(tree,X,...){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	if(hasArg(quiet)) quiet<-list(...)$quiet
	else quiet<-FALSE
	if(hasArg(rotate)) rotate<-list(...)$rotate
	else rotate<-TRUE
	if(hasArg(edge.constraint))edge.constraint<-list(...)$edge.constraint
	else edge.constraint<-tree$edge[,2]
	if(hasArg(time.constraint)) time.constraint<-list(...)$time.constraint
	else time.constraint<-c(0,max(nodeHeights(tree)))
	if(any(time.constraint<0)) time.constraint<-max(nodeHeights(tree))+time.constraint
	if(length(time.constraint)==1) time.constraint<-rep(time.constraint,2)
	if(!is.matrix(X)) X<-as.matrix(X)
	tip<-setdiff(rownames(X),tree$tip.label)
	fossilML(tree,X,quiet,tip,edge.constraint,time.constraint,plot,search,rotate)
}

fossilML<-function(tree,X,quiet,tip,edge.constraint,time.constraint,plot,search,rotate){
	if(!quiet) cat(paste("Optimizing the phylogenetic position of ",tip," using ML. Please wait....\n",sep=""))
	if(ncol(X)>1&&rotate){
		pca<-phyl.pca(tree,X[tree$tip.label,])
		obj<-phyl.vcv(X[tree$tip.label,],vcv(tree),1)
		X<-(X-matrix(rep(obj$a[,1],nrow(X)),nrow(X),ncol(X),byrow=TRUE))%*%pca$Evec
	}
	lik.tree<-function(par,tip,tree,edge,height,XX,rotate,time.constraint){
		tip.depth<-par[1]
		tip.height<-par[2]
		if(tip.height>(tip.depth+1e-6)){
			ii<-which(tree$edge[,2]==edge)
			tree<-bind.tip(tree,tip,where=edge,position=min(height[2]-tip.depth,tree$edge.length[ii]),
				edge.length=tip.height-tip.depth)
			if(!rotate) XX<-phyl.pca(tree,XX[tree$tip.label,])$S
			obj<-phyl.vcv(as.matrix(XX[tree$tip.label,]),vcv(tree),1)
			ll<-vector()
			for(i in 1:ncol(XX)) 
				ll[i]<-sum(dmnorm(XX[tree$tip.label,i],mean=rep(obj$a[i,1],nrow(XX)),obj$C*obj$R[i,i],log=TRUE)) 
			if(plot){
				plotTree(tree,mar=c(2.1,0.1,3.1,0.1),ftype="i",fsize=0.8)
				obj<-lapply(time.constraint,function(x,tree) lines(rep(x,2),
					c(0,Ntip(tree)+1),col="red",lty="dashed"),tree=tree)
				title(paste("log(L) = ",signif(sum(ll),6),sep=""))
				axis(1)
			}
			logL<-sum(ll)
		} else logL<--.Machine$double.xmax/1e100
		logL
	}
	ee<-intersect(tree$edge[,2],edge.constraint)
	H<-nodeHeights(tree)
	hh<-unique(c(tree$edge[H[,1]<=time.constraint[2],2],tree$edge[H[,2]<=time.constraint[2],2]))
	ee<-intersect(ee,hh)
	fit<-vector(mode="list",length=length(ee))
	for(i in 1:length(ee)){
		ii<-which(tree$edge[,2]==ee[i])
		lower<-c(H[ii,1],time.constraint[1])
		upper<-c(H[ii,2],time.constraint[2])
		par<-c(mean(lower),mean(upper))
		fit[[i]]<-optim(par,lik.tree,tip=tip,tree=tree,edge=ee[i],height=H[ii,],XX=X,rotate=rotate,
			time.constraint=time.constraint,
			method="L-BFGS-B",lower=lower,upper=upper,control=list(fnscale=-1))
	}
	logL<-sapply(fit,function(x) x$value)
	ii<-which(logL==max(logL))
	if(length(ii)==1){
		fit<-fit[[ii]]
		edge<-ee[ii]
		mltree<-bind.tip(tree,tip,where=edge,
			position=H[which(tree$edge[,2]==edge),2]-fit$par[1],
			edge.length=fit$par[2]-fit$par[1])
		mltree$logL<-fit$value
	} else {
		fit<-fit[ii]
		mltree<-vector(mode="list",length=length(ii))
		for(i in 1:length(ii)){
			edge<-ee[ii[i]]
			mltree[[i]]<-bind.tip(tree,tip,where=edge,
				position=H[which(tree$edge[,2]==edge),2]-fit[[i]]$par[1],
				edge.length=fit[[i]]$par[2]-fit[[i]]$par[1])
			mltree[[i]]$logL<-fit[[i]]$value
		}
		class(mltree)<-"multiPhylo"
	}
	if(!quiet) cat("Done.\n")
	if(plot){
		if(inherits(mltree,"phylo")){
			plotTree(mltree,mar=c(2.1,0.1,3.1,0.1),ftype="i",fsize=0.8)
			obj<-lapply(time.constraint,function(x,tree) lines(rep(x,2),
				c(0,Ntip(tree)+1),col="red",lty="dashed"),tree=tree)
			title(paste("Optimized position of taxon \"",tip,"\"",sep=""))
			axis(1)
		} else {
			par(mfrow=c(1,length(mltree)))
			for(i in 1:length(mltree)){
				plotTree(mltree[[i]],mar=c(2.1,0.1,0.1,0.1),ftype="i",fsize=0.8)
				obj<-lapply(time.constraint,function(x,tree) lines(rep(x,2),
					c(0,Ntip(tree)+1),col="red",lty="dashed"),tree=tree)
				axis(1)
			}
		}					
	}
	mltree
}
