## function to compute the marginal posterior probabilities for nodes using the rerooting method
## written by Liam J. Revell 2013, 2015, 2017, 2018

rerootingMethod<-function(tree,x,model=c("ER","SYM"),...){
	if(!inherits(tree,"phylo")) 
		stop("tree should be an object of class \"phylo\".")
	if(hasArg(tips)) tips<-list(...)$tips
	else tips<-NULL
	if(!is.matrix(model)) model<-model[1]
	n<-Ntip(tree)
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
	YY<-fitMk(tree,yy,model=model,output.liks=TRUE,...)
	Q<-matrix(c(0,YY$rates)[YY$index.matrix+1],length(YY$states),
		length(YY$states),dimnames=list(YY$states,YY$states))
	diag(Q)<--colSums(Q,na.rm=TRUE)
	nn<-if(tips) c(1:n,2:tree$Nnode+n) else 2:tree$Nnode+n
	ff<-function(nn){
		tt<-if(nn>Ntip(tree)) root(tree,node=nn) else reroot(tree,nn,tree$edge.length[which(tree$edge[,2]==nn)])
		fitMk(tt,yy,model=model,fixedQ=Q,output.liks=TRUE)$lik.anc[1,]
	}
	XX<-t(sapply(nn,ff))
	if(tips) XX<-rbind(XX[1:n,],YY$lik.anc[1,],XX[(n+1):nrow(XX),])
	else XX<-rbind(yy,YY$lik.anc[1,],XX)
	rownames(XX)<-1:(tree$Nnode+n)
	if(tips) rownames(XX)[1:n]<-tree$tip.label
	XX<-if(tips) XX else XX[1:tree$Nnode+n,]
	obj<-list(loglik=YY$logLik,Q=Q,marginal.anc=XX,tree=tree,x=yy)
	class(obj)<-"rerootingMethod"
	obj
}

print.rerootingMethod<-function(x,digits=6,printlen=NULL,...){
	cat("Ancestral character estimates using re-rooting method\nof Yang et al. (1995):\n")
	if(is.null(printlen)) print(round(x$marginal.anc,digits)) else { 
		print(round(x$marginal.anc[1:printlen,],digits))
		cat("...\n")
	}
	cat("\nEstimated transition matrix,\nQ =\n")
	print(round(x$Q,digits))
	cat("\n**Note that if Q is not symmetric the marginal\nreconstructions may be invalid.\n")
	cat(paste("\nLog-likelihood =",round(x$loglik,digits),"\n\n"))
}

plot.rerootingMethod<-function(x, ...){
	args<-list(...)
	if(is.null(args$lwd)) args$lwd<-1
	if(is.null(args$ylim)) args$ylim<-c(-0.1*Ntip(x$tree),Ntip(x$tree))
	if(is.null(args$offset)) args$offset<-0.5
	if(is.null(args$ftype)) args$ftype="i"
	args$tree<-x$tree	
	do.call(plotTree,args)
	if(hasArg(piecol)) piecol<-list(...)$piecol
	else piecol<-setNames(colorRampPalette(c("blue",
		"yellow"))(ncol(x$marginal.anc)),
		colnames(x$marginal.anc))
	if(hasArg(node.cex)) node.cex<-list(...)$node.cex
	else node.cex<-0.6
	nodelabels(pie=x$marginal.anc[
		as.character(1:x$tree$Nnode+Ntip(x$tree)),],
		piecol=piecol,cex=node.cex)
	if(hasArg(tip.cex)) tip.cex<-list(...)$tip.cex
	else tip.cex<-0.4
	tiplabels(pie=x$x[x$tree$tip.label,],piecol=piecol,
		cex=tip.cex)
	legend(x=par()$usr[1],y=par()$usr[1],
		legend=colnames(x$marginal.anc),
		pch=21,pt.bg=piecol,pt.cex=2.2,bty="n")
}

logLik.rerootingMethod<-function(object,...) object$loglik
