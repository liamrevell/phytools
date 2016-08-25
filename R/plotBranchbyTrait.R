# function to plot probability or trait value by branch
# written by Liam J. Revell 2013, 2014

plotBranchbyTrait<-function(tree,x,mode=c("edges","tips","nodes"),palette="rainbow",legend=TRUE,xlims=NULL,...){
	mode<-mode[1]
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(mode=="tips"){
		x<-c(x[tree$tip.label],fastAnc(tree,x))
		names(x)[1:length(tree$tip.label)]<-1:length(tree$tip.label)
		XX<-matrix(x[tree$edge],nrow(tree$edge),2)
		x<-rowMeans(XX)
	} else if(mode=="nodes"){
		XX<-matrix(x[tree$edge],nrow(tree$edge),2)
		x<-rowMeans(XX)
	}
	# begin optional arguments
	if(hasArg(tol)) tol<-list(...)$tol
	else tol<-1e-6
	if(hasArg(prompt)) prompt<-list(...)$prompt
	else prompt<-FALSE
	if(hasArg(type)) type<-list(...)$type
	else type<-"phylogram"
	if(hasArg(show.tip.label)) show.tip.label<-list(...)$show.tip.label
	else show.tip.label<-TRUE
	if(hasArg(show.node.label)) show.node.label<-list(...)$show.node.label
	else show.node.label<-FALSE
	if(hasArg(edge.width)) edge.width<-list(...)$edge.width
	else edge.width<-4
	if(hasArg(edge.lty)) edge.lty<-list(...)$edge.lty
	else edge.lty<-1
	if(hasArg(font)) font<-list(...)$font
	else font<-3
	if(hasArg(cex)) cex<-list(...)$cex
	else cex<-par("cex")
	if(hasArg(adj)) adj<-list(...)$adj
	else adj<-NULL
	if(hasArg(srt)) srt<-list(...)$srt
	else srt<-0
	if(hasArg(no.margin)) no.margin<-list(...)$no.margin
	else no.margin<-TRUE
	if(hasArg(root.edge)) root.edge<-list(...)$root.edge
	else root.edge<-FALSE
	if(hasArg(label.offset)) label.offset<-list(...)$label.offset
	else label.offset<-0.01*max(nodeHeights(tree))
	if(hasArg(underscore)) underscore<-list(...)$underscore
	else underscore<-FALSE
	if(hasArg(x.lim)) x.lim<-list(...)$x.lim
	else x.lim<-NULL
	if(hasArg(y.lim)) y.lim<-list(...)$y.lim
	else y.lim<-if(legend&&!prompt&&type%in%c("phylogram","cladogram")) 
		c(1-0.06*length(tree$tip.label),length(tree$tip.label)) else NULL
	if(hasArg(direction)) direction<-list(...)$direction
	else direction<-"rightwards"
	if(hasArg(lab4ut)) lab4ut<-list(...)$lab4ut
	else lab4ut<-NULL
	if(hasArg(tip.color)) tip.color<-list(...)$tip.color
	else tip.color<-"black"
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	if(hasArg(rotate.tree)) rotate.tree<-list(...)$rotate.tree
	else rotate.tree<-0
	if(hasArg(open.angle)) open.angle<-list(...)$open.angle
	else open.angle<-0
	# end optional arguments
	if(palette=="heat.colors") cols<-heat.colors(n=1000)
	if(palette=="gray") cols<-gray(1000:1/1000)
	if(palette=="rainbow")	cols<-rainbow(1000,start=0.7,end=0) # blue->red
	if(is.null(xlims)) xlims<-range(x)+c(-tol,tol)
	breaks<-0:1000/1000*(xlims[2]-xlims[1])+xlims[1]
	whichColor<-function(p,cols,breaks){
		i<-1
		while(p>=breaks[i]&&p>breaks[i+1]) i<-i+1
		cols[i]
	}
	colors<-sapply(x,whichColor,cols=cols,breaks=breaks)
	par(lend=2)
	# now plot
	xx<-plot.phylo(tree,type=type,show.tip.label=show.tip.label,show.node.label=show.node.label,edge.color=colors,
	edge.width=edge.width,edge.lty=edge.lty,font=font,cex=cex,adj=adj,srt=srt,no.margin=no.margin,root.edge=root.edge,
	label.offset=label.offset,underscore=underscore,x.lim=x.lim,y.lim=y.lim,direction=direction,lab4ut=lab4ut,
	tip.color=tip.color,plot=plot,rotate.tree=rotate.tree,open.angle=open.angle,lend=2,new=FALSE)
	if(legend==TRUE&&is.logical(legend)) legend<-round(0.3*max(nodeHeights(tree)),2)
	if(legend){
		if(hasArg(title)) title<-list(...)$title
		else title<-"trait value"
		if(hasArg(digits)) digits<-list(...)$digits
		else digits<-1
		if(prompt) add.color.bar(legend,cols,title,xlims,digits,prompt=TRUE)
		else add.color.bar(legend,cols,title,xlims,digits,prompt=FALSE,x=par()$usr[1]+0.05*(par()$usr[2]-par()$usr[1]),y=par()$usr[3]+0.05*(par()$usr[4]-par()$usr[3]))
	}
	invisible(xx)
}

# function to add color bar
# written by Liam J. Revell 2013, 2015, 2016

add.color.bar<-function(leg,cols,title=NULL,lims=c(0,1),digits=1,prompt=TRUE,lwd=4,outline=TRUE,...){
	if(prompt){
		cat("Click where you want to draw the bar\n")
		flush.console()
		x<-unlist(locator(1))
		y<-x[2]
		x<-x[1]
	} else {
		if(hasArg(x)) x<-list(...)$x
		else x<-0
		if(hasArg(y)) y<-list(...)$y
		else y<-0
	}
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-1.0
	if(hasArg(subtitle)) subtitle<-list(...)$subtitle
	else subtitle<-NULL
	if(hasArg(direction)) direction<-list(...)$direction
	else direction<-"rightwards"
	if(direction%in%c("rightwards","leftwards")){
		X<-x+cbind(0:(length(cols)-1)/length(cols),1:length(cols)/length(cols))*(leg)
		if(direction=="leftwards"){ 
			X<-X[nrow(X):1,]
			if(!is.null(lims)) lims<-lims[2:1]
		}
		Y<-cbind(rep(y,length(cols)),rep(y,length(cols)))
	} else if(direction%in%c("upwards","downwards")){
		Y<-y+cbind(0:(length(cols)-1)/length(cols),1:length(cols)/length(cols))*(leg)
		if(direction=="downwards"){ 
			X<-X[nrow(X):1,]
			if(!is.null(lims)) lims<-lims[2:1]
		}
		X<-cbind(rep(x,length(cols)),rep(x,length(cols)))
	}
	if(outline) lines(c(X[1,1],X[nrow(X),2]),c(Y[1,1],Y[nrow(Y),2]),lwd=lwd+2,lend=2) 
	for(i in 1:length(cols)) lines(X[i,],Y[i,],col=cols[i],lwd=lwd,lend=2)
	if(direction%in%c("rightwards","leftwards")){
		if(!is.null(lims)) text(x=x,y=y,
			round(lims[1],digits),pos=3,cex=fsize)
		if(!is.null(lims)) text(x=x+leg,y=y+0.02*diff(par()$usr[3:4]),
			round(lims[2],digits),pos=3,cex=fsize)
		if(is.null(title)) title<-"P(state=1)"
		text(x=(2*x+leg)/2,y=y,title,pos=3,cex=fsize)
		if(is.null(subtitle)) 
			text(x=(2*x+leg)/2,y=y,paste("length=",round(leg,3),sep=""),pos=1,cex=fsize)
		else text(x=(2*x+leg)/2,y=y,subtitle,pos=1,cex=fsize)
	} else if(direction%in%c("upwards","downwards")){
		if(!is.null(lims)) text(x=x,y=y-0.02*diff(par()$usr[3:4]),round(lims[1],digits),
			pos=1,cex=fsize)
		if(!is.null(lims)) text(x=x,y=y+leg+0.02*diff(par()$usr[3:4]),
			round(lims[2],digits),
			pos=3,cex=fsize)
		if(is.null(title)) title<-"P(state=1)"
		text(x=x-0.04*diff(par()$usr[1:2]),y=(2*y+leg)/2,title,
			pos=3,cex=fsize,srt=90)
		if(is.null(subtitle)) 
			text(x=x+0.04*diff(par()$usr[1:2]),y=(2*y+leg)/2,
				paste("length=",round(leg,3),sep=""),pos=1,
				srt=90,cex=fsize)
		else text(x=x+0.04*diff(par()$usr[1:2]),y=(2*y+leg)/2,
			subtitle,pos=1,cex=fsize,srt=90)
	}
}
