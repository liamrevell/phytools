## functions to split a tree into two (or more) plots by different methods
## written by Liam J. Revell 2012, 2015

plotTree.splits<-function(tree,splits=NULL,file=NULL,fn=NULL,...){
	if(is.null(fn)) fn<-function(){}
	ef<-0.037037037037 ## expansion factor
	if(!is.null(file)) pdf(file,width=8.5,height=11)
	if(is.null(splits)) splits<-(floor(0.5*Ntip(tree))+0.5)/Ntip(tree)
	if(hasArg(y.lim)) y.lim<-list(...)$y.lim
	else y.lim<-c(0,Ntip(tree))
	S<-matrix(c(0,splits,splits,1+1/Ntip(tree)),length(splits)+1,2)
	S<-cbind(S[,1]+ef*(S[,2]-S[,1]),S[,2]-ef*(S[,2]-S[,1]))
	S<-S*diff(y.lim)+y.lim[1]
	for(i in nrow(S):1){
		if(is.null(file)&&i<nrow(S)) par(ask=TRUE)
		plotTree(tree,ylim=S[i,],...)
		fn()
	}
	if(!is.null(file)) oo<-dev.off()
}

splitplotTree<-function(tree,fsize=1.0,ftype="reg",lwd=2,split=NULL,new.window=FALSE){
	# check font
	ftype<-which(c("off","reg","b","i","bi")==ftype)-1
	if(!ftype) fsize=0 	# check tree
	# check tree
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	# rescale tree
	tree$edge.length<-tree$edge.length/max(nodeHeights(tree))
	# check split
	n<-length(tree$tip)
	if(is.null(split)) 
		split<-0.5
	s<-(1-split)*(n+1)
	# reorder
	cw<-reorder(tree)
	pw<-reorder(tree,"pruningwise")
	# count nodes
	m<-cw$Nnode
	# y coordinates for nodes
	y<-matrix(NA,m+n,1)
	y[cw$edge[cw$edge[,2]<=length(cw$tip),2]]<-1:n
	nodes<-unique(pw$edge[,1])
	for(i in 1:m){
		desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
		y[nodes[i]]<-(min(y[desc])+max(y[desc]))/2
	}
	Y<-matrix(NA,nrow(cw$edge),ncol(cw$edge))
	for(i in 1:nrow(cw$edge)) Y[i,]<-rep(y[cw$edge[i,2]])
	# compute node heights
	H<-nodeHeights(tree)
	# now split
	a<-which(Y[,1]>=s)
	Y1<-Y[a,]; Y2<-Y[-a,]
	b<-(min(Y1)-max(Y2))/2
	Y1<-Y1-s+1; Y2<-Y2
	H1<-H[a,]; H2<-H[-a,]
	edge1<-cw$edge[a,]; edge2<-cw$edge[-a,]
	if(length(H1)==2){ H1<-matrix(H1,1,2); Y1<-matrix(Y1,1,2); edge1<-matrix(edge1,1,2) }
	if(length(H2)==2){ H2<-matrix(H2,1,2); Y2<-matrix(Y2,1,2); edge2<-matrix(edge2,1,2) }
	# label offset
	offset<-0.2*lwd/3+0.2/3
	# open plot
	par(mar=c(0.1,0.1,0.1,0.1))
	if(!new.window) layout(matrix(c(1,2),1,2))
	# first half
	plot.new();
	if(fsize*max(strwidth(cw$tip.label))<1.0){
		k<-(1-fsize*max(strwidth(cw$tip.label)))/max(H)
		H<-k*H
		H1<-k*H1
		H2<-k*H2
	} else message("Font size too large to properly rescale tree to window.")
	plot.window(xlim=c(0,max(H)+fsize*max(strwidth(cw$tip.label))),ylim=c(1-b,max(rbind(Y1,Y2))))
	for(i in 1:nrow(H1)) lines(H1[i,],Y1[i,],lwd=lwd,lend=2)
	nodes<-unique(edge1[,1])
	for(i in 1:length(nodes)) lines(H1[which(edge1[,1]==nodes[i]),1],Y1[which(edge1[,1]==nodes[i]),1],lwd=lwd)
	for(i in 1:nrow(edge1)) if(edge1[i,1]%in%edge2[,1]) lines(c(H1[i,1],H1[i,1]),c(Y1[i,1],1-b),lwd=lwd,lend=2)
	tips<-edge1[edge1[,2]<=n,2]
	for(i in 1:length(tips)) 
		if(ftype) text(H1[which(edge1[,2]==tips[i]),2],Y1[which(edge1[,2]==tips[i]),1],cw$tip.label[tips[i]],pos=4,cex=fsize,font=ftype,offset=offset)
	# second half
	if(new.window){ dev.new(); par(mar=c(0.1,0.1,0.1,0.1)) }
	plot.new()
	if(max(Y1)>max(Y2)) Y2<-Y2+max(Y1)-max(Y2)
	plot.window(xlim=c(0,max(H)+fsize*max(strwidth(cw$tip.label))),ylim=c(1,max(rbind(Y1,Y2)+b)))
	for(i in 1:nrow(H2)) lines(H2[i,],Y2[i,],lwd=lwd,lend=2)
	nodes<-unique(edge2[,1])
	for(i in 1:length(nodes)) lines(H2[which(edge2[,1]==nodes[i]),1],Y2[which(edge2[,1]==nodes[i]),1],lwd=lwd)
	for(i in 1:nrow(edge2)) if(edge2[i,1]%in%edge1[,1]) lines(c(H2[i,1],H2[i,1]),c(Y2[i,1],max(Y2)+b),lwd=lwd,lend=2)
	tips<-edge2[edge2[,2]<=n,2]
	for(i in 1:length(tips)) 
		if(ftype) text(H2[which(edge2[,2]==tips[i]),2],Y2[which(edge2[,2]==tips[i]),1],cw$tip.label[tips[i]],pos=4,cex=fsize,font=ftype,offset=offset)
	# reset margin and layout
	layout(1)
	par(mar=c(5,4,4,2)+0.1)
}
