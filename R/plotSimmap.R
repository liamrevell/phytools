## functions plot stochastic character mapped trees
## written by Liam Revell 2011-2020

plotSimmap<-function(tree,colors=NULL,fsize=1.0,ftype="reg",lwd=2,
	pts=FALSE,node.numbers=FALSE,mar=NULL,add=FALSE,offset=NULL,direction="rightwards",
	type="phylogram",setEnv=TRUE,part=1.0,xlim=NULL,ylim=NULL,nodes="intermediate",
	tips=NULL,maxY=NULL,hold=TRUE,split.vertical=FALSE,lend=2,asp=NA,outline=FALSE,
	plot=TRUE){
	if(inherits(tree,"multiPhylo")){
		par(ask=TRUE)
		for(i in 1:length(tree)) plotSimmap(tree[[i]],colors=colors,fsize=fsize,
			ftype=ftype,lwd=lwd,pts=pts,node.numbers=node.numbers,mar,add,offset,
			direction,type,setEnv,part,xlim,ylim,nodes,tips,maxY,hold,split.vertical,
			lend,asp,outline,plot)
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
			if(direction%in%c("upwards","downwards")){
				if(outline){
					fg<-par()$fg
					par(fg="transparent")
					black<-colors
					black[]<-"black"
					updownPhylogram(tree,colors=black,fsize,ftype,lwd=lwd+2,pts,
						node.numbers,mar,add,offset,direction,setEnv,xlim,ylim,nodes,
						tips,split.vertical,lend,asp,plot)
					par(fg=fg)
				}
				updownPhylogram(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,
				add=if(outline) TRUE else add,offset,direction,setEnv,xlim,ylim,nodes,
				tips,split.vertical,lend,asp,plot)
			} else {
				if(outline){
					fg<-par()$fg
					par(fg="transparent")
					black<-colors
					black[]<-"black"
					plotPhylogram(tree,colors=black,fsize,ftype,lwd=lwd+2,pts,
						node.numbers,mar,add,offset,direction,setEnv,xlim,ylim,nodes,
						tips,split.vertical,lend,asp,plot)
					par(fg=fg)
				}
				plotPhylogram(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,
				add=if(outline) TRUE else add,offset,direction,setEnv,xlim,ylim,nodes,
				tips,split.vertical,lend,asp,plot)
			}
		} else if(type=="fan"){
			if(outline){
				fg<-par()$fg
				par(fg="transparent")
				black<-colors
				black[]<-"black"
				plotFan(tree,colors=black,fsize,ftype,lwd=lwd+2,mar,add,part,setEnv,
					xlim,ylim,tips,maxY,lend,plot,offset)
				par(fg=fg)
			}
			plotFan(tree,colors,fsize,ftype,lwd,mar,add=if(outline) TRUE else add,part,
				setEnv,xlim,ylim,tips,maxY,lend,plot,offset)
		} else if(type=="cladogram"){
			if(outline){
				fg<-par()$fg
				par(fg="transparent")
				black<-colors
				black[]<-"black"
				plotCladogram(tree,colors=black,fsize,ftype,lwd=lwd+2,mar,add,offset,
					direction,xlim,ylim,nodes,tips,lend,asp,plot)
				par(fg=fg)
			}
			plotCladogram(tree,colors,fsize,ftype,lwd,mar,add=if(outline) TRUE else add,
				offset,direction,xlim,ylim,nodes,tips,lend,asp,plot)
		}
		if(hold) null<-dev.flush()
	}
}

## function to plot simmap tree in type "phylogram"
## written by Liam J. Revell 2011-2017
updownPhylogram<-function(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,
	add,offset,direction,setEnv,xlim,ylim,placement,tips,split.vertical,lend,
	asp,plot){
	if(split.vertical&&!setEnv){
		cat("split.vertical requires setEnv=TRUE. Setting split.vertical to FALSE.\n")
		spit.vertical<-FALSE
	}
	# set offset fudge (empirically determined)
	offsetFudge<-1.37
	# reorder
	cw<-reorderSimmap(tree)
	pw<-reorderSimmap(tree,"postorder")
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
			desc<-getDescendants(pw,nodes[i])
			desc<-desc[desc<=Ntip(pw)]
			Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
		} else if(placement=="weighted"){
			desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
			n1<-desc[which(Y[desc]==min(Y[desc]))]
			n2<-desc[which(Y[desc]==max(Y[desc]))]
			v1<-pw$edge.length[which(pw$edge[,2]==n1)]
			v2<-pw$edge.length[which(pw$edge[,2]==n2)]
			Y[nodes[i]]<-((1/v1)*Y[n1]+(1/v2)*Y[n2])/(1/v1+1/v2)
		} else if(placement=="inner"){
			desc<-getDescendants(pw,nodes[i])
			desc<-desc[desc<=Ntip(pw)]
			mm<-which(abs(Y[desc]-median(Y[1:Ntip(pw)]))==min(abs(Y[desc]-
				median(Y[1:Ntip(pw)]))))
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
	if(is.null(ylim)){
		pp<-par("pin")[2]
		sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
			offsetFudge*fsize*strwidth("W",units="inches")
		alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=H,sw=sw,pp=pp,
			interval=c(0,1e6))$minimum
		ylim<-if(direction=="downwards") c(min(H)-sw/alp,max(H)) else c(min(H),max(H)+sw/alp)
	}
	if(is.null(xlim)) xlim=range(Y)
	if(direction=="downwards") H<-max(H)-H
	plot.window(xlim=xlim,ylim=ylim,asp=asp)
	####
	if(plot){
		if(!split.vertical){
			for(i in 1:m) lines(Y[cw$edge[which(cw$edge[,1]==nodes[i]),2]],
				H[which(cw$edge[,1]==nodes[i]),1],
				col=colors[names(cw$maps[[match(nodes[i],
				cw$edge[,1])]])[1]],lwd=lwd)
		}
		for(i in 1:nrow(cw$edge)){
			x<-H[i,1]
			for(j in 1:length(cw$maps[[i]])){
				if(direction=="downwards")
					lines(c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),c(x,x-cw$maps[[i]][j]),
						col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
				else lines(c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),c(x,x+cw$maps[[i]][j]),
						col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
				if(pts) points(c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),c(x,x+cw$maps[[i]][j]),
					pch=20,lwd=(lwd-1))
				x<-x+if(direction=="downwards") -cw$maps[[i]][j] else cw$maps[[i]][j]
				j<-j+1
			}
		}
		if(node.numbers){
			symbols(mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),
				if(direction=="downwards") max(H) else 0,
				rectangles=matrix(c(1.2*fsize*strwidth(as.character(Ntip(cw)+1)),
				1.4*fsize*strheight(as.character(Ntip(cw)+1))),1,2),inches=FALSE,
				bg="white",add=TRUE)
			text(mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),
				if(direction=="downwards") max(H) else 0,Ntip(cw)+1,
				cex=fsize)
			for(i in 1:nrow(cw$edge)){
				x<-H[i,2]
				if(cw$edge[i,2]>Ntip(cw)){
					symbols(Y[cw$edge[i,2]],x,
						rectangles=matrix(c(1.2*fsize*strwidth(as.character(cw$edge[i,2])),
						1.4*fsize*strheight(as.character(cw$edge[i,2]))),1,2),inches=FALSE,
						bg="white",add=TRUE)
					text(Y[cw$edge[i,2]],x,cw$edge[i,2],cex=fsize)
				}
			}
		}
		if(direction=="downwards") pos<-if(par()$usr[3]>par()$usr[4]) 2 else 4
		if(direction=="upwards") pos<-if(par()$usr[3]>par()$usr[4]) 2 else 4
		for(i in 1:n){
			shift<-offset*fsize*strwidth("W")*(diff(par()$usr[3:4])/diff(par()$usr[1:2]))
			if((direction=="downwards"&&diff(par()$usr[3:4])>0) ||
				(direction=="upwards"&&diff(par()$usr[3:4])<0)) shift<--shift
			if(ftype){
				text(labels=cw$tip.label[i],Y[i],
					H[which(cw$edge[,2]==i),2]+shift,
					pos=pos,offset=0,cex=fsize,font=ftype,
					srt=if(direction=="downwards") 270 else 90)
			}
		}
	}
	if(setEnv){
		PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
			show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
			font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
			x.lim=xlim,y.lim=ylim,
			direction=direction,tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
			edge=cw$edge,xx=Y[,1],yy=sapply(1:(Ntip(cw)+cw$Nnode),
			function(x,y,z) y[match(x,z)],y=H,z=cw$edge))
		assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
	}
	if(plot) if(split.vertical) splitEdgeColor(cw,colors,lwd)
}

# function to plot simmap tree in type "phylogram"
# written by Liam J. Revell 2011-2017
plotPhylogram<-function(tree,colors,fsize,ftype,lwd,pts,node.numbers,mar,
	add,offset,direction,setEnv,xlim,ylim,placement,tips,split.vertical,lend,
	asp,plot){
	if(split.vertical&&!setEnv){
		cat("split.vertical requires setEnv=TRUE. Setting split.vertical to FALSE.\n")
		spit.vertical<-FALSE
	}
	# set offset fudge (empirically determined)
	offsetFudge<-1.37
	# reorder
	cw<-reorderSimmap(tree)
	pw<-reorderSimmap(tree,"postorder")
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
			desc<-getDescendants(pw,nodes[i])
			desc<-desc[desc<=Ntip(pw)]
			Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
		} else if(placement=="weighted"){
			desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
			n1<-desc[which(Y[desc]==min(Y[desc]))]
			n2<-desc[which(Y[desc]==max(Y[desc]))]
			v1<-pw$edge.length[which(pw$edge[,2]==n1)]
			v2<-pw$edge.length[which(pw$edge[,2]==n2)]
			Y[nodes[i]]<-((1/v1)*Y[n1]+(1/v2)*Y[n2])/(1/v1+1/v2)
		} else if(placement=="inner"){
			desc<-getDescendants(pw,nodes[i])
			desc<-desc[desc<=Ntip(pw)]
			mm<-which(abs(Y[desc]-median(Y[1:Ntip(pw)]))==min(abs(Y[desc]-
				median(Y[1:Ntip(pw)]))))
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
		xlim<-if(direction=="leftwards") c(min(H)-sw/alp,max(H)) else c(min(H),max(H)+sw/alp)
	}
	if(is.null(ylim)) ylim=range(Y)
	if(direction=="leftwards") H<-max(H)-H
	plot.window(xlim=xlim,ylim=ylim,asp=asp)
	if(plot){
		####
		if(!split.vertical){
			for(i in 1:m) lines(H[which(cw$edge[,1]==nodes[i]),1],
				Y[cw$edge[which(cw$edge[,1]==nodes[i]),2]],col=colors[names(cw$maps[[match(nodes[i],
				cw$edge[,1])]])[1]],lwd=lwd)
		}
		for(i in 1:nrow(cw$edge)){
			x<-H[i,1]
			for(j in 1:length(cw$maps[[i]])){
				if(direction=="leftwards")
					lines(c(x,x-cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
						col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
				else lines(c(x,x+cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
						col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
				if(pts) points(c(x,x+cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
					pch=20,lwd=(lwd-1))
				x<-x+if(direction=="leftwards") -cw$maps[[i]][j] else cw$maps[[i]][j]
				j<-j+1
			}
		}
		if(node.numbers){
			symbols(if(direction=="leftwards") max(H) else 0,
				mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),
				rectangles=matrix(c(1.2*fsize*strwidth(as.character(Ntip(cw)+1)),
				1.4*fsize*strheight(as.character(Ntip(cw)+1))),1,2),inches=FALSE,
				bg="white",add=TRUE)
			text(if(direction=="leftwards") max(H) else 0,
				mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),Ntip(cw)+1,
				cex=fsize)
			for(i in 1:nrow(cw$edge)){
				x<-H[i,2]
				if(cw$edge[i,2]>Ntip(cw)){
					symbols(x,Y[cw$edge[i,2]],
						rectangles=matrix(c(1.2*fsize*strwidth(as.character(cw$edge[i,2])),
						1.4*fsize*strheight(as.character(cw$edge[i,2]))),1,2),inches=FALSE,
						bg="white",add=TRUE)
					text(x,Y[cw$edge[i,2]],cw$edge[i,2],cex=fsize)
				}
			}
		}
		if(direction=="leftwards") pos<-if(par()$usr[1]>par()$usr[2]) 4 else 2
		if(direction=="rightwards") pos<-if(par()$usr[1]>par()$usr[2]) 2 else 4
		for(i in 1:n) if(ftype) text(H[which(cw$edge[,2]==i),2],Y[i],cw$tip.label[i],pos=pos,
			offset=offset,cex=fsize,font=ftype)
	}
	if(setEnv){
		PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
			show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
			font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
			x.lim=xlim,y.lim=ylim,
			direction=direction,tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
			edge=cw$edge,xx=sapply(1:(Ntip(cw)+cw$Nnode),
			function(x,y,z) y[match(x,z)],y=H,z=cw$edge),yy=Y[,1])
		assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
	}
	if(plot) if(split.vertical) splitEdgeColor(cw,colors,lwd)
}

# function to plot simmap tree in type "fan"
# written by Liam J. Revell 2013-2017
plotFan<-function(tree,colors,fsize,ftype,lwd,mar,add,part,setEnv,xlim,ylim,tips,maxY,lend,plot,offset){
	if(!plot) cat("plot=FALSE option is not permitted for type=\"fan\". Tree will be plotted.\n")
	if(is.null(offset)) offset<-1
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
	Y<-part*cbind(Y[as.character(cw$edge[,2])],Y[as.character(cw$edge[,2])])
	R<-nodeHeights(cw)
	# now put into a circular coordinate system
	x<-R*cos(Y)
	y<-R*sin(Y)
	# optimize x & y limits
	par(mar=mar)
	offsetFudge<-1.37 # empirically determined
	OFFSET<-0
	pp<-par("pin")[1]
 	sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
		offsetFudge*OFFSET*fsize*strwidth("W",units="inches") 
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
	## first, the lines emerging from the root (if there are only two):
	jj<-which(cw$edge[,1]==(Ntip(cw)+1))
	if(length(jj)==2){
		m.left<-cumsum(cw$maps[[jj[1]]])/sum(cw$maps[[jj[1]]])
		xx.left<-c(x[jj[1],1],x[jj[1],1]+(x[jj[1],2]-x[jj[1],1])*m.left)
		yy.left<-c(y[jj[1],1],y[jj[1],1]+(y[jj[1],2]-y[jj[1],1])*m.left)
		m.right<-cumsum(cw$maps[[jj[2]]])/sum(cw$maps[[jj[2]]])
		xx.right<-c(x[jj[2],1],x[jj[2],1]+(x[jj[2],2]-x[jj[2],1])*m.right)
		yy.right<-c(y[jj[2],1],y[jj[2],1]+(y[jj[2],2]-y[jj[2],1])*m.right)
		xx<-c(xx.left[length(xx.left):1],xx.right[2:length(xx.right)])
		yy<-c(yy.left[length(yy.left):1],yy.right[2:length(yy.right)])
		col<-colors[c(names(m.left)[length(m.left):1],names(m.right))]
		segments(xx[2:length(xx)-1],yy[2:length(yy)-1],xx[2:length(xx)],yy[2:length(yy)],
			col=col,lwd=lwd,lend=lend)
	} else jj<-NULL
	for(i in 1:nrow(cw$edge)){
		if(i%in%jj==FALSE){
			maps<-cumsum(cw$maps[[i]])/sum(cw$maps[[i]])
			xx<-c(x[i,1],x[i,1]+(x[i,2]-x[i,1])*maps)
			yy<-c(y[i,1],y[i,1]+(y[i,2]-y[i,1])*maps)
			for(i in 1:(length(xx)-1)) lines(xx[i+0:1],yy[i+0:1],col=colors[names(maps)[i]],
				lwd=lwd,lend=lend)
		}
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
		tt<-if(aa>90&&aa<270) paste(cw$tip.label[i],paste(rep(" ",offset),
			collapse=""),sep="") else paste(paste(rep(" ",offset),collapse=""),
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
			xx=c(x[sapply(1:n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2],x[1,1],
			if(m>1) x[sapply(2:m+n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2] else c()),
			yy=c(y[sapply(1:n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2],y[1,1],
			if(m>1) y[sapply(2:m+n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2] else c()))
		assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
	}
}

## internal function for slanted cladogram
## written by Liam J. Revell 2017
plotCladogram<-function(tree,colors=NULL,fsize=1.0,ftype="reg",lwd=2,mar=NULL,
	add=FALSE,offset=NULL,direction="rightwards",xlim=NULL,ylim=NULL,
	nodes="intermediate",tips=NULL,lend=2,asp=NA,plot=TRUE){
	placement<-nodes
	# set offset fudge (empirically determined)
	offsetFudge<-1.37
	# reorder
	cw<-reorderSimmap(tree)
	pw<-reorderSimmap(tree,"postorder")
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
			desc<-getDescendants(pw,nodes[i])
			desc<-desc[desc<=Ntip(pw)]
			Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
		} else if(placement=="weighted"){
			desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
			n1<-desc[which(Y[desc]==min(Y[desc]))]
			n2<-desc[which(Y[desc]==max(Y[desc]))]
			v1<-pw$edge.length[which(pw$edge[,2]==n1)]
			v2<-pw$edge.length[which(pw$edge[,2]==n2)]
			Y[nodes[i]]<-((1/v1)*Y[n1]+(1/v2)*Y[n2])/(1/v1+1/v2)
		} else if(placement=="inner"){
			desc<-getDescendants(pw,nodes[i])
			desc<-desc[desc<=Ntip(pw)]
			mm<-which(abs(Y[desc]-median(Y[1:Ntip(pw)]))==min(abs(Y[desc]-
				median(Y[1:Ntip(pw)]))))
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
		xlim<-if(direction=="leftwards") c(min(H)-sw/alp,max(H)) else c(min(H),max(H)+sw/alp)
	}
	if(is.null(ylim)) ylim=range(Y)
	if(direction=="leftwards") H<-max(H)-H
	plot.window(xlim=xlim,ylim=ylim,asp=asp)
	if(plot){
		####
		for(i in 1:nrow(cw$edge)){
			x<-H[i,1]
			y<-Y[cw$edge[i,1]]
			m<-(Y[cw$edge[i,2]]-Y[cw$edge[i,1]])/(H[i,2]-H[i,1])
			if(is.finite(m)){
				for(j in 1:length(cw$maps[[i]])){
					if(direction=="leftwards")
						lines(c(x,x-cw$maps[[i]][j]),c(y,y-cw$maps[[i]][j]*m),
							col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
					else lines(c(x,x+cw$maps[[i]][j]),c(y,y+cw$maps[[i]][j]*m),
						col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
					x<-x+if(direction=="leftwards") -cw$maps[[i]][j] else cw$maps[[i]][j]
					y<-y+if(direction=="leftwards") -m*cw$maps[[i]][j] else m*cw$maps[[i]][j]
					j<-j+1
				}
			} else {
				lines(rep(x,2),Y[cw$edge[i,]],col=colors[names(cw$maps[[i]])[1]],lwd=lwd,
					lend=lend)
			}
		}
		if(direction=="leftwards") pos<-if(par()$usr[1]>par()$usr[2]) 4 else 2
		if(direction=="rightwards") pos<-if(par()$usr[1]>par()$usr[2]) 2 else 4
		for(i in 1:n) if(ftype) text(H[which(cw$edge[,2]==i),2],Y[i],cw$tip.label[i],pos=pos,
			offset=offset,cex=fsize,font=ftype)
	}
	PP<-list(type="phylogram",use.edge.length=TRUE,node.pos=1,
		show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
		font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
		x.lim=xlim,y.lim=ylim,
		direction=direction,tip.color="black",Ntip=Ntip(cw),Nnode=cw$Nnode,
		edge=cw$edge,xx=sapply(1:(Ntip(cw)+cw$Nnode),
		function(x,y,z) y[match(x,z)],y=H,z=cw$edge),yy=Y[,1])
	assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
}


## adds legend to an open stochastic map style plot
## written by Liam J. Revell 2013, 2016, 2017
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
	w<-par()$mfcol[2]*h*abs(diff(par()$usr[1:2])/diff(par()$usr[3:4]))
	flipped<-par()$usr[1]>par()$usr[2]
	if(vertical){
		y<-y-0:(length(leg)-1)*1.5*h
		x<-rep(x+w/2,length(y))
		text(x + if(flipped) -w else w,y,leg,pos=4,cex=fsize/par()$cex)
	} else {
		sp<-abs(fsize*max(strwidth(leg)))
		x<-x + if(flipped) w/2-0:(length(leg)-1)*1.5*(sp+w) else -w/2+0:(length(leg)-1)*1.5*(sp+w)
		y<-rep(y+w/2,length(x))
		text(x,y,leg,pos=4,cex=fsize/par()$cex)
	}
	if(shape=="square") symbols(x,y,squares=rep(w,length(x)),bg=colors,add=TRUE,inches=FALSE)
	else if(shape=="circle") nulo<-mapply(draw.circle,x=x,y=y,col=colors,
		MoreArgs=list(nv=200,radius=w/2))
	else stop(paste("shape=\"",shape,"\" is not a recognized option.",sep=""))
}

# function plots a tree; in the new version this is just a wrapper for plotSimmap
# written by Liam Revell 2012-2017
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
	if(hasArg(lend)) lend<-list(...)$lend
	else lend<-2
	if(hasArg(asp)) asp<-list(...)$asp
	else asp<-NA
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-TRUE
	if(inherits(tree,"multiPhylo")){
		par(ask=TRUE)
		if(!is.null(color)) names(color)<-"1"
		for(i in 1:length(tree)) plotTree(tree[[i]],color=color,fsize=fsize,ftype=ftype,
			lwd=lwd,pts=pts,node.numbers=node.numbers,mar=mar,add=add,offset=offset,
			direction=direction,type=type,setEnv=setEnv,part=part,xlim=xlim,ylim=ylim,
			nodes=nodes,tips=tips,maxY=maxY,hold=hold,lend=lend,asp=asp,plot=plot)
	} else {
		if(is.null(tree$edge.length)) tree<-compute.brlen(tree)
		tree$maps<-as.list(tree$edge.length)
		for(i in 1:length(tree$maps)) names(tree$maps[[i]])<-c("1")
		if(!is.null(color)) names(color)<-"1"
		plotSimmap(tree,colors=color,fsize=fsize,ftype=ftype,lwd=lwd,pts=pts,
			node.numbers=node.numbers,mar=mar,add=add,offset=offset,direction=direction,
			type=type,setEnv=setEnv,part=part,xlim=xlim,ylim=ylim,nodes=nodes,tips=tips,maxY=maxY,
			hold=hold,lend=lend,asp=asp,plot=plot)
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
		cols<-vector()
		for(j in 1:length(daughters)){
			jj<-which(tree$edge[,2]==daughters[j])
			cols[j]<-if(tree$maps[[jj]][1]==0&&length(tree$maps[[jj]])>1) colors[names(tree$maps[[jj]])[2]] 
				else colors[names(tree$maps[[jj]])[1]]
		}
		ii<-order(obj$yy[c(i,daughters)])
		jj<-order(obj$yy[daughters])
		x0<-x1<-rep(obj$xx[i],length(daughters))
		y0<-obj$yy[c(i,daughters)][ii][1:length(daughters)]
		y1<-obj$yy[c(i,daughters)][ii][2:(length(daughters)+1)]
		cols<-cols[jj]
		for(j in 1:length(x0)) segments(x0[j],y0[j],x1[j],y1[j],col=cols[j],lwd=lwd,lend=2)
	}
}


