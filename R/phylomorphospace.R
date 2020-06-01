## this funciton creates a phylomorphospace plot (Sidlauskas 2006)
## written by Liam J. Revell 2010-13, 2015, 2018

phylomorphospace<-function(tree,X,A=NULL,label=c("radial","horizontal","off"),control=list(),...){
	# some minor error checking
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(nrow(X)!=length(tree$tip)) stop("X must contain the same number of rows as species in tree.")
	if(is.null(rownames(X))){
		warning("X is missing row names; assuming order of tip labels.")
		rownames(X)<-tree$tip.label
	}
	if(ncol(X)!=2){ 
		warning("X has more than 2 columns.  Using only the first 2 columns.")
		X<-X[,1:2]
	}
	# get ancestral states
	if(is.null(A)) A<-apply(X,2,fastAnc,tree=tree)
	# control list
	con=list(col.edge=setNames(rep("black",nrow(tree$edge)),as.character(tree$edge[,2])),
		col.node=setNames(rep("black",max(tree$edge)),as.character(1:max(tree$edge))))
	con[(namc<-names(control))]<-control
	# get optional arguments
	if(hasArg(node.by.map)) node.by.map<-list(...)$node.by.map
	else node.by.map<-FALSE
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-"reg"
	ftype<-which(c("off","reg","b","i","bi")==ftype)-1
	if(!ftype) label<-"off"
	if(hasArg(node.size)){ 
		node.size<-list(...)$node.size
		if(length(node.size)==1) node.size<-c(node.size,node.size)
	} else node.size<-c(1,1.3)
	# set xlim & ylim
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-expand(range(c(X[,1],A[,1])),1.1)
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-expand(range(c(X[,2],A[,2])),1.1)
	# set xlab & ylab
	if(hasArg(xlab)) xlab<-list(...)$xlab
	else xlab<-colnames(X)[1]
	if(hasArg(ylab)) ylab<-list(...)$ylab
	else ylab<-colnames(X)[2]
	# set font size for tip labels
	if(hasArg(fsize)) fsize<-0.75*list(...)$fsize
	else fsize<-0.75	
	# check if colors for stochastic mapping
	if(hasArg(colors)) colors<-list(...)$colors
	else if(!is.null(tree$maps)) colors<-setNames(palette()[1:ncol(tree$mapped.edge)],
		sort(colnames(tree$mapped.edge)))
	# set lwd
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-if(is.null(tree$maps)) 1 else 2
	# other optional arguments?
	if(hasArg(axes)) axes<-list(...)$axes
	else axes<-TRUE
	if(hasArg(add)) add<-list(...)$add
	else add<-FALSE
	if(hasArg(pch)) pch<-list(...)$pch
	else pch<-21
	if(hasArg(bty)) bty<-list(...)$bty
	else bty<-"o"
	# deprecate to logical label argument
	label<-label[1]
	if(label==TRUE||label==FALSE) 
		message("options for label have changed.\nNow choose \"radial\", \"horizontal\", or \"off\".")
	if(label==TRUE) label<-"radial"
	if(label==FALSE) label<-"off"
	# do some bookkeeping
	aa<-setNames(c(X[tree$tip.label,1],A[,1]),c(1:length(tree$tip.label),rownames(A)))
	bb<-setNames(c(X[tree$tip.label,2],A[,2]),c(1:length(tree$tip.label),rownames(A)))
	XX<-matrix(aa[as.character(tree$edge)],nrow(tree$edge),2)
	YY<-matrix(bb[as.character(tree$edge)],nrow(tree$edge),2)
	# plot projection
	if(!add) plot(x=A[1,1],y=A[1,2],xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,pch=16,cex=0.1,
		col="white",axes=axes,frame.plot=TRUE,bty=bty)
	if(is.null(tree$maps)){
		for(i in 1:nrow(XX)) lines(XX[i,],YY[i,],col=con$col.edge[as.character(tree$edge[i,2])],
			lwd=lwd)
	} else {
		for(i in 1:nrow(XX)){
			xx<-tree$maps[[i]]/sum(tree$maps[[i]])*(XX[i,2]-XX[i,1])
			yy<-tree$maps[[i]]/sum(tree$maps[[i]])*(YY[i,2]-YY[i,1])
			cc<-names(tree$maps[[i]])
			x<-XX[i,1]; y<-YY[i,1]
			for(j in 1:length(xx)){
				lines(c(x,x+xx[j]),c(y,y+yy[j]),col=colors[cc[j]],lwd=lwd)
				x<-x+xx[j]; y<-y+yy[j]
			}
		}
		if(node.by.map){
			zz<-c(getStates(tree,type="tips"),getStates(tree))
			names(zz)[1:length(tree$tip.label)]<-sapply(names(zz)[1:length(tree$tip.label)],
				function(x,y) which(y==x),y=tree$tip.label)
			con$col.node<-setNames(colors[zz],names(zz))
		}
	}
	zz<-c(tree$edge[1,1],tree$edge[tree$edge[,2]>length(tree$tip.label),2])
	points(c(XX[1,1],XX[tree$edge[,2]>length(tree$tip.label),2]),c(YY[1,1],
		YY[tree$edge[,2]>length(tree$tip.label),2]),pch=pch,cex=node.size[1],
		col=if(pch%in%c(1:20)) con$col.node[as.character(zz)] else "black",
		bg=if(pch%in%c(21:25)) con$col.node[as.character(zz)] else NULL)
	zz<-tree$edge[tree$edge[,2]<=length(tree$tip.label),2]
	points(XX[tree$edge[,2]<=length(tree$tip.label),2],YY[tree$edge[,2]<=length(tree$tip.label),2],
		pch=pch,cex=node.size[2],col=if(pch%in%c(1:20)) con$col.node[as.character(zz)] else "black",
		bg=if(pch%in%c(21:25)) con$col.node[as.character(zz)] else NULL)
	zz<-sapply(1:length(tree$tip.label),function(x,y) which(x==y),y=tree$edge[,2])
	if(label!="off"){
		asp<-(par()$usr[2]-par()$usr[1])/(par()$usr[4]-par()$usr[3])
		for(i in 1:length(tree$tip.label)){
			ii<-which(tree$edge[,2]==i)
			aa<-atan(asp*(YY[ii,2]-YY[ii,1])/(XX[ii,2]-XX[ii,1]))/(2*pi)*360
			adj<-if(XX[ii,2]<XX[ii,1]) c(1,0.25) else c(0,0.25)
			tt<-if(XX[ii,2]<XX[ii,1]) paste(tree$tip.label[i],"  ",sep="") else 
				paste("  ",tree$tip.label[i],sep="")
			if(label=="radial") text(XX[ii,2],YY[ii,2],tt,cex=fsize,srt=aa,adj=adj,font=ftype)
			else if(label=="horizontal") text(XX[ii,2],YY[ii,2],tt,cex=fsize,adj=adj,font=ftype)
		}
	} else adj<-0
	xx<-setNames(c(XX[1,1],XX[,2]),c(tree$edge[1,1],tree$edge[,2]))
	xx<-xx[order(as.numeric(names(xx)))]
	yy<-setNames(c(YY[1,1],YY[,2]),c(tree$edge[1,1],tree$edge[,2]))
	yy<-yy[order(as.numeric(names(yy)))]	
	PP<-list(type="phylomorphospace",use.edge.length=TRUE,node.pos=1,
		show.tip.label=if(label!="off") TRUE else FALSE,show.node.label=FALSE,
		font=ftype,cex=fsize,adj=adj,srt=NULL,no.margin=FALSE,label.offset=0,
		x.lim=par()$usr[1:2],y.lim=par()$usr[3:4],
		direction=NULL,tip.color="black",Ntip=Ntip(tree),Nnode=tree$Nnode,
		edge=tree$edge,xx=xx,yy=yy)
	assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
}

# function to expand a range
# written by Liam J. Revell 2013
expand<-function(x,factor=1.05){
	rr<-x[2]-x[1]
	mm<-mean(x)
	x[1]<-mm-rr*factor/2
	x[2]<-mm+rr*factor/2
	return(x)
}

