## functions & methods to simulate & plot discrete time Brownian motion 
## on a simulated discrete-time tree

simBMphylo<-function(n,t,sig2,plot=TRUE,...){
	if(length(sig2)!=t) sig2<-rep(sig2[1],t)
	b<-exp((log(n)-log(2))/t)-1
	tree<-pbtree(b=b,d=0,t=t,n=n,type="discrete",
		tip.label=if(n<=26) LETTERS[n:1] else NULL,
		quiet=TRUE)
	H<-nodeHeights(tree)
	root<-Ntip(tree)+1
	xx<-list()
	for(i in 1:nrow(tree$edge)){
		sd<-sqrt(sig2[H[i,1]+1:tree$edge.length[i]])
		x<-rnorm(n=tree$edge.length[i],sd=sd)
		x<-c(0,cumsum(x))
		if(tree$edge[i,1]!=root){
			ii<-which(tree$edge[,2]==tree$edge[i,1])
			x<-x+xx[[ii]][length(xx[[ii]])]
		}
		xx[[i]]<-x
	}
	object<-list(tree=tree,x=xx)
	class(object)<-"simBMphylo"
	if(plot) plot(object,...)
	invisible(object)
}

plot.simBMphylo<-function(x,...){
	xx<-x$x
	tree<-x$tree
	H<-nodeHeights(tree)
	layout(mat=matrix(c(1,2),2,1),
		heights=c(1/3,2/3))
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-0.9
	if(hasArg(cex.axis)) cex.axis<-list(...)$cex.axis
	else cex.axis<-0.9
	if(hasArg(cex.lab)) cex.lab=list(...)$cex.lab
	else cex.lab<-1
	plotTree(tree,mar=c(0.1,4.1,2.1,1.1),
		xlim=c(0,1.05*max(H)),
		ylim=c(1-2*(Ntip(tree)-1)*0.04,
		Ntip(tree)+(Ntip(tree)-1)*0.04),lwd=1,
		fsize=fsize)
	## axis(1,cex.axis=cex.axis)
	mtext("a)",line=0,adj=0,cex=cex.lab)
	plot.new()
	par(mar=c(5.1,4.1,1.1,1.1))
	plot.window(xlim=c(0,1.05*max(H)),ylim=range(xx))
	axis(1,cex.axis=cex.axis)
	axis(2,cex.axis=cex.axis)
	for(i in 1:length(xx))
		lines(H[i,1]:H[i,2],xx[[i]])
	for(i in 1:Ntip(tree)){
		ii<-which(tree$edge[,2]==i)
		text(max(H),xx[[ii]][length(xx[[ii]])],
			tree$tip.label[i],pos=4,offset=0.4/3,
			cex=fsize)
	}
	mtext("b)",line=0,adj=0,cex=cex.lab)

	title(xlab="time",ylab="phenotype",
		cex.lab=cex.lab)
}

print.simBMphylo<-function(x,...){
	cat(paste("\nObject of class \"simBMphylo\" with",Ntip(x$tree),"taxa.\n"),
		sep="")
	cat("To print use plot method.\n\n")
}

