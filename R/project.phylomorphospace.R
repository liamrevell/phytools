## morph a phylogeny to a phylomorphospace, vice versa, or there & back again
## written by Liam J. Revell 2018

project.phylomorphospace<-function(tree,X,nsteps=200,sleep=0,
	direction=c("to","from","both"),...){
	direction<-direction[1]
	tree<-minRotate(reorder(tree,"cladewise"),X[,2],print=FALSE)
	X<-X[tree$tip.label,]
	A<-cbind(fastAnc(tree,X[,1]),fastAnc(tree,X[,2]))
	cladogram<-tree
	cladogram$edge.length<-NULL
	mar<-par()$mar
	plotTree(cladogram,type="cladogram",nodes="centered",plot=FALSE)
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	obj$xx<-obj$xx*(max(c(X[,1],A[,1]))-min(c(X[,1],A[,1])))+
		min(c(X[,1],A[,1]))
	xlim<-range(obj$xx)
	xlim[2]<-xlim[2]+abs(diff(xlim))*(obj$x.lim[2]-1)
	obj$yy<-(obj$yy-min(obj$yy))/(max(obj$yy)-min(obj$yy))*
		(max(c(X[,2],A[,2]))-min(c(X[,2],A[,2])))+min(c(X[,2],A[,2]))
	ylim<-range(obj$yy)
	X0<-cbind(obj$xx[1:Ntip(tree)],obj$yy[1:Ntip(tree)])
	rownames(X0)<-tree$tip.label
	A0<-cbind(obj$xx[1:tree$Nnode+Ntip(tree)],
		obj$yy[1:tree$Nnode+Ntip(tree)])
	rownames(A0)<-1:tree$Nnode+Ntip(tree)
	par(mar=mar,new=TRUE)
	dev.hold()
	if(direction%in%c("to","both"))
		phylomorphospace(tree,X0,A0,label="horizontal",xlim=xlim,
			ylim=ylim,...)
	else if(direction=="from")
		phylomorphospace(tree,X,A,label="radial",xlim=xlim,
			ylim=ylim,...)
	dev.flush()
	if(direction=="both") nsteps<-ceiling(nsteps/2)
	if(direction%in%c("to","both")){
		for(i in 2:nsteps){
			Sys.sleep(sleep)
			dev.hold()
			phylomorphospace(tree,((nsteps-i)*X0+i*X)/nsteps,
				((nsteps-i)*A0+i*A)/nsteps,xlim=xlim,ylim=ylim,
					...)
			dev.flush()
		}
	}
	if(direction%in%c("from","both")){
		for(i in (nsteps-1):1){
			Sys.sleep(sleep)
			dev.hold()
			phylomorphospace(tree,((nsteps-i)*X0+i*X)/nsteps,
				((nsteps-i)*A0+i*A)/nsteps,xlim=xlim,ylim=ylim,
					...)
			dev.flush()
		}
		dev.hold()
		phylomorphospace(tree,X0,A0,label="horizontal",xlim=xlim,
			ylim=ylim,...)
		dev.flush()
	}
	invisible(NULL)
}
