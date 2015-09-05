## functions plot stochastic character mapped trees
## written by Liam Revell 2011-2015

plotSimmap<-function(tree,colors=NULL,fsize=1.0,ftype="reg",lwd=2,
	pts=FALSE,node.numbers=FALSE,mar=NULL,add=FALSE,offset=NULL,direction="rightwards",
	type="phylogram",setEnv=TRUE,part=1.0,xlim=NULL,ylim=NULL,nodes="intermediate",
	tips=NULL,maxY=NULL,hold=TRUE,split.vertical=FALSE){
	if(inherits(tree,"multiPhylo")){
		par(ask=TRUE)
		for(i in 1:length(tree)) plotSimmap(tree[[i]],colors=colors,fsize=fsize,ftype=ftype,
			lwd=lwd,pts=pts,node.numbers=node.numbers,mar,add,offset,direction,type,
			setEnv,part,xlim,ylim,nodes,tips,maxY,hold)
	} else {
		# check tree
		if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\"")
		if(is.null(tree$maps)) stop("tree should contain mapped states on edges.")
		# check font
		ftype<-which(c("off","reg","b","i","bi")==ftype)-1
		if(!ftype) fsize=0 
		# check colors
		if(is.null(colors)){
			st<-sort(unique(unlist(sapply(tree$maps,names))))
			colors<-palette()[1:length(st)]
			names(colors)<-st
			if(length(st)>1){
				cat("no colors provided. using the following legend:\n")
				print(colors)
			}
		}
		# swap out "_" character for spaces (assumes _ is a place holder)
		tree$tip.label<-gsub("_"," ",tree$tip.label)
		# get margin
		if(is.null(mar)) mar=rep(0.1,4)
		if(hold) null<-dev.hold()
		if(type=="phylogram"){
			plotPhylogram(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,add,offset,
				direction,setEnv,xlim,ylim,nodes,tips,split.vertical)
		} else if(type=="fan"){
			plotFan(tree,colors,fsize,ftype,lwd,mar,add,part,setEnv,xlim,ylim,tips,maxY)
		}
		if(hold) null<-dev.flush()
	}
}	

# function to plot simmap tree in type "phylogram"
# written by Liam J. Revell 2011-2015
plotPhylogram<-function(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,
	add,offset,direction,setEnv,xlim,ylim,placement,tips,split.vertical){
	# set offset fudge (empirically determined)
	if(split.vertical&&!setEnv){
		cat("split.vertical requires setEnv=TRUE. Setting split.vertical to FALSE.\n")
		spit.vertical<-FALSE
	}
	offsetFudge<-1.37
	# reorder
	cw<-reorderSimmap(tree)
	pw<-reorderSimmap(tree,"pruningwise")
	# count nodes and tips
	n<-Ntip(cw)
 	m<-cw$Nnode
	# Y coordinates for nodes
	Y<-matrix(NA,m+n,1)
	# first, assign y coordinates to all the tip nodes
	if(is.null(tips)) Y[cw$edge[cw$edge[,2]<=n,2]]<-1:n
	else Y[cw$edge[cw$edge[,2]<=n,2]]<-if(is.null(names(tips))) 
		tips[sapply(1:Ntip(cw),function(x,y) which(y==x),y=cw$edge[cw$edge[,2]<=n,2])]
		else tips[gsub(" ","_",cw$tip.label)]
	# get Y coordinates of the nodes
	nodes<-unique(pw$edge[,1])
	for(i in 1:m){
		if(placement=="intermediate"){ 
			desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
			Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
		} else if(placement=="centered"){
			desc<-getDescendants(tree,nodes[i])
			desc<-desc[desc<=Ntip(tree)]
			Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
		} else if(placement=="weighted"){
			desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
			n1<-desc[which(Y[desc]==min(Y[desc]))]
			n2<-desc[which(Y[desc]==max(Y[desc]))]
			v1<-tree$edge.length[which(tree$edge[,2]==n1)]
			v2<-tree$edge.length[which(tree$edge[,2]==n2)]
			Y[nodes[i]]<-((1/v1)*Y[n1]+(1/v2)*Y[n2])/(1/v1+1/v2)
		} else if(placement=="inner"){
			desc<-getDescendants(tree,nodes[i])
			desc<-desc[desc<=Ntip(tree)]
			mm<-which(abs(Y[desc]-median(Y[1:Ntip(tree)]))==min(abs(Y[desc]-
				median(Y[1:Ntip(tree)]))))
			if(length(mm>1)) mm<-mm[which(Y[desc][mm]==min(Y[desc][mm]))]
			Y[nodes[i]]<-Y[desc][mm]
		}
	}
	# compute node heights
	H<-nodeHeights(cw)
	# open plot
	par(mar=mar)
	if(is.null(offset)) offset<-0.2*lwd/3+0.2/3
	if(!add) plot.new()
	###
	if(is.null(xlim)){
		pp<-par("pin")[1]
		sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
			offsetFudge*fsize*strwidth("W",units="inches")
		alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=H,sw=sw,pp=pp,
			interval=c(0,1e6))$minimum
		xlim<-c(min(H),max(H)+sw/alp)
	}
	if(is.null(ylim)) ylim=range(Y)
	if(direction=="leftwards") plot.window(xlim=xlim[2:1],ylim=ylim)
	else plot.window(xlim=xlim,ylim=ylim)
	####
	if(!split.vertical){
		for(i in 1:m) lines(H[which(cw$edge[,1]==nodes[i]),1],
			Y[cw$edge[which(cw$edge[,1]==nodes[i]),2]],col=colors[names(cw$maps[[match(nodes[i],
			cw$edge[,1])]])[1]],lwd=lwd)
	}
	for(i in 1:nrow(cw$edge)){
		x<-H[i,1]
 		for(j in 1:length(cw$maps[[i]])){
			lines(c(x,x+cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
				col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=2)
			if(pts) points(c(x,x+cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
				pch=20,lwd=(lwd-1))
			x<-x+cw$maps[[i]][j]; j<-j+1
		}
	}
	if(node.numbers){
		symbols(0,mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),
			rectangles=matrix(c(1.2*fsize*strwidth(as.character(Ntip(cw)+1)),
			1.4*fsize*strheight(as.character(Ntip(cw)+1))),1,2),inches=FALSE,
			bg="white",add=TRUE)
		text(0,mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),Ntip(cw)+1,
			cex=fsize)
		for(i in 1:nrow(cw$edge)){
			x<-H[i,2]
			if(cw$edge[i,2]>Ntip(tree)){
				symbols(x,Y[cw$edge[i,2]],
					rectangles=matrix(c(1.2*fsize*strwidth(as.character(cw$edge[i,2])),
					1.4*fsize*strheight(as.character(cw$edge[i,2]))),1,2),inches=FALSE,
					bg="white",add=TRUE)
				text(x,Y[cw$edge[i,2]],cw$edge[i,2],cex=fsize)
			}
		}
	}
	pos<-if(direction=="leftwards") 2 else 4
	for(i in 1:n) if(ftype) text(H[which(cw$edge[,2]==i),2],Y[i],cw$tip.label[i],pos=pos,
		offset=offset,cex=fsize,font=ftype)
	if(setEnv){
		PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
			show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
			font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
			x.lim=par()$usr[1:2],y.lim=par()$usr[3:4],
			direction=direction,tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
			edge=cw$edge,xx=sapply(1:(Ntip(cw)+cw$Nnode),
			function(x,y,z) y[match(x,z)],y=H,z=cw$edge),yy=Y[,1])
		assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
	}
	if(split.vertical) splitEdgeColor(cw,colors,lwd)
}

# function to plot simmap tree in type "fan"
# written by Liam J. Revell 2013-2015
plotFan<-function(tree,colors,fsize,ftype,lwd,mar,add,part,setEnv,xlim,ylim,tips,maxY){
	# reorder
	cw<-reorder(tree)
	pw<-reorder(tree,"pruningwise")
	# count nodes and tips
	n<-Ntip(cw)
	m<-cw$Nnode 
	# get Y coordinates on uncurved space
	Y<-vector(length=m+n)
	if(is.null(tips)) tips<-1:n
	if(part<1.0) Y[cw$edge[cw$edge[,2]<=n,2]]<-0:(n-1)
	else Y[cw$edge[cw$edge[,2]<=n,2]]<-tips
	nodes<-unique(pw$edge[,1])
	for(i in 1:m){
		desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
		Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
	}
	if(is.null(maxY)) maxY<-max(Y)
	Y<-setNames(Y/maxY*2*pi,1:(n+m))
	Y<-part*cbind(Y[as.character(tree$edge[,2])],Y[as.character(tree$edge[,2])])
	R<-nodeHeights(cw)
	# now put into a circular coordinate system
	x<-R*cos(Y)
	y<-R*sin(Y)
	# optimize x & y limits
	par(mar=mar)
	offsetFudge<-1.37 # empirically determined
	offset<-0
	pp<-par("pin")[1]
 	sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
		offsetFudge*offset*fsize*strwidth("W",units="inches") 
	alp<-optimize(function(a,H,sw,pp) (2*a*1.04*max(H)+2*sw-pp)^2,H=R,sw=sw,pp=pp,
		interval=c(0,1e6))$minimum
	if(part<=0.25) x.lim<-y.lim<-c(0,max(R)+sw/alp)
	else if(part>0.25&&part<=0.5){ 
		x.lim<-c(-max(R)-sw/alp,max(R)+sw/alp)
		y.lim<-c(0,max(R)+sw/alp)
	} else x.lim<-y.lim<-c(-max(R)-sw/alp,max(R)+sw/alp)
	if(is.null(xlim)) xlim<-x.lim
	if(is.null(ylim)) ylim<-y.lim
	# plot tree
	if(!add) plot.new()
	plot.window(xlim=xlim,ylim=ylim,asp=1)
	# plot radial lines (edges)
	for(i in 1:nrow(cw$edge)){
		maps<-cumsum(cw$maps[[i]])/sum(cw$maps[[i]])
		xx<-c(x[i,1],x[i,1]+(x[i,2]-x[i,1])*maps)
		yy<-c(y[i,1],y[i,1]+(y[i,2]-y[i,1])*maps)
		for(i in 1:(length(xx)-1)) lines(xx[i+0:1],yy[i+0:1],col=colors[names(maps)[i]],
			lwd=lwd,lend=2)
	}
	# plot circular lines
	for(i in 1:m+n){
		r<-R[match(i,cw$edge)]
		a1<-min(Y[which(cw$edge==i)])
		a2<-max(Y[which(cw$edge==i)])
		draw.arc(0,0,r,a1,a2,lwd=lwd,col=colors[names(cw$maps[[match(i,cw$edge[,1])]])[1]])
	}
	# plot labels
	for(i in 1:n){
		ii<-which(cw$edge[,2]==i)
		aa<-Y[ii,2]/(2*pi)*360
		adj<-if(aa>90&&aa<270) c(1,0.25) else c(0,0.25)
		tt<-if(aa>90&&aa<270) paste(cw$tip.label[i]," ",sep="") else paste(" ",
			cw$tip.label[i],sep="")
		aa<-if(aa>90&&aa<270) 180+aa else aa
		if(ftype) text(x[ii,2],y[ii,2],tt,srt=aa,adj=adj,cex=fsize,font=ftype)
	}
	if(setEnv){
		PP<-list(type="fan",use.edge.length=TRUE,node.pos=1,
			show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
			font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
			x.lim=xlim,y.lim=ylim,direction="rightwards",tip.color="black",
			Ntip=Ntip(cw),Nnode=cw$Nnode,edge=cw$edge,
			xx=c(x[sapply(1:n,function(x,y) which(x==y)[1],y=tree$edge[,2]),2],x[1,1],
			if(m>1) x[sapply(2:m+n,function(x,y) which(x==y)[1],y=tree$edge[,2]),2] else c()),
			yy=c(y[sapply(1:n,function(x,y) which(x==y)[1],y=tree$edge[,2]),2],y[1,1],
			if(m>1) y[sapply(2:m+n,function(x,y) which(x==y)[1],y=tree$edge[,2]),2] else c()))
		assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
	}
}

# adds legend to an open stochastic map style plot
# written by Liam J. Revell 2013
add.simmap.legend<-function(leg=NULL,colors,prompt=TRUE,vertical=TRUE,...){
	if(hasArg(shape)) shape<-list(...)$shape
	else shape<-"square"
	if(prompt){
		cat("Click where you want to draw the legend\n")
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
	if(is.null(leg)) leg<-names(colors)
	h<-fsize*strheight(LETTERS[1])
	w<-h*(par()$usr[2]-par()$usr[1])/(par()$usr[4]-par()$usr[3])
	if(vertical){
		y<-y-0:(length(leg)-1)*1.5*h
		x<-rep(x+w/2,length(y))		
		text(x+w,y,leg,pos=4,cex=fsize/par()$cex)
	} else {
		sp<-fsize*max(strwidth(leg))
		x<-x-w/2+0:(length(leg)-1)*1.5*(sp+w)
		y<-rep(y+w/2,length(x))
		text(x,y,leg,pos=4,cex=fsize/par()$cex)
	}
	if(shape=="square") symbols(x,y,squares=rep(w,length(x)),bg=colors,add=TRUE,inches=FALSE)
	else if(shape=="circle") symbols(x,y,circles=rep(w,length(x)),bg=colors,add=TRUE,
		inches=FALSE)
	else stop(paste("shape=\"",shape,"\" is not a recognized option.",sep=""))
}

# function plots a tree; in the new version this is just a wrapper for plotSimmap
# written by Liam Revell 2012-2015
plotTree<-function(tree,...){
	if(hasArg(color)) color<-list(...)$color
	else color<-NULL
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-1.0
	if(hasArg(ftype)) ftype<-list(...)$ftype
	else ftype<-"reg"
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(hasArg(pts)) pts<-list(...)$pts
	else pts<-FALSE
	if(hasArg(node.numbers)) node.numbers<-list(...)$node.numbers
	else node.numbers<-FALSE
	if(hasArg(mar)) mar<-list(...)$mar
	else mar<-NULL
	if(hasArg(add)) add<-list(...)$add
	else add<-FALSE
	if(hasArg(offset)) offset<-list(...)$offset
	else offset<-NULL
	if(hasArg(type)) type<-list(...)$type
	else type<-"phylogram"
	if(hasArg(direction)) direction<-list(...)$direction
	else direction<-"rightwards"
	if(hasArg(setEnv)) setEnv<-list(...)$setEnv
	else setEnv<-TRUE
	if(hasArg(part)) part<-list(...)$part
	else part<-1.0
	if(hasArg(xlim)) xlim<-list(...)$xlim
	else xlim<-NULL
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-NULL
	if(hasArg(nodes)) nodes<-list(...)$nodes
	else nodes<-"intermediate"
	if(hasArg(tips)) tips<-list(...)$tips
	else tips<-NULL
	if(hasArg(maxY)) maxY<-list(...)$maxY
	else maxY<-NULL
	if(hasArg(hold)) hold<-list(...)$hold
	else hold<-TRUE
	if(inherits(tree,"multiPhylo")){
		par(ask=TRUE)
		if(!is.null(color)) names(color)<-"1"
		for(i in 1:length(tree)) plotTree(tree[[i]],color=color,fsize=fsize,ftype=ftype,
			lwd=lwd,pts=pts,node.numbers=node.numbers,mar=mar,add=add,offset=offset,
			direction=direction,type=type,setEnv=setEnv,part=part,xlim=xlim,ylim=ylim,
			nodes=nodes,tips=tips,maxY=maxY,hold=hold)
	} else {
		if(is.null(tree$edge.length)) tree<-compute.brlen(tree)
		tree$maps<-as.list(tree$edge.length)
		for(i in 1:length(tree$maps)) names(tree$maps[[i]])<-c("1")
		if(!is.null(color)) names(color)<-"1"
		plotSimmap(tree,colors=color,fsize=fsize,ftype=ftype,lwd=lwd,pts=pts,
			node.numbers=node.numbers,mar=mar,add=add,offset=offset,direction=direction,
			type=type,setEnv=setEnv,part=part,xlim=xlim,ylim=ylim,nodes=nodes,tips=tips,maxY=maxY,
			hold=hold)
	}
}

## S3 method for objects of class "simmap" & "multiSimmap"
## added by Liam J. Revell 2015
plot.simmap<-function(x,...) plotSimmap(x,...)
plot.multiSimmap<-function(x,...) plotSimmap(x,...)

## function to split vertical plotted lines by the states of daughter edges
## written by Liam J. Revell 2015
splitEdgeColor<-function(tree,colors,lwd=2){
    obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
    for(i in 1:tree$Nnode+Ntip(tree)){
        daughters<-tree$edge[which(tree$edge[,1]==i),2]
        for(j in 1:length(daughters)){
            jj<-which(tree$edge[,2]==daughters[j])
            color<-if(tree$maps[[jj]][1]==0) colors[names(tree$maps[[jj]])[2]] else colors[names(tree$maps[[jj]])[1]]
            lines(rep(obj$xx[i],2),obj$yy[c(i,daughters[j])],col=color,lwd=lwd)
        }
    }
}

