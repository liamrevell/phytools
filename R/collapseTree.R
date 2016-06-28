## function to interactively expand and contract subtrees on a phylogeny
## inspired by the phylogeny interface of sharksrays.org by Gavin Naylor
## written by Liam J. Revell 2015

collapseTree<-function(tree,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(hasArg(nodes)) nodes<-list(...)$nodes
	else nodes<-TRUE
	if(hasArg(hold)) hold<-list(...)$hold
	else hold<-TRUE
	if(hasArg(drop.extinct)) drop.extinct<-list(...)$drop.extinct
	else drop.extinct=FALSE
	cat("Click on the nodes that you would like to collapse...\n")
	## turn off locator bell (it's annoying)
	options(locatorBell=FALSE)
	## check for node labels
	if(is.null(tree$node.label)) tree$node.label<-as.character(Ntip(tree)+1:tree$Nnode)
	else if(any(tree$node.label=="")){
		tree$node.label[which(tree$node.label)==""]<-
			which(tree$node.label=="")+Ntip(tree)
	}
	## remove any spaces
	tree$node.label<-sapply(tree$node.label,gsub,pattern=" ",replacement="_")
	tree$tip.label<-sapply(tree$tip.label,gsub,pattern=" ",replacement="_")
	## copy original tree:
	otree<-tree<-reorder(tree)
	## plot initial tree:
	if(hold) dev.hold()
	fan(tree,...)
	lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	circles(x=lastPP$xx[1:Ntip(tree)],y=lastPP$yy[1:Ntip(tree)],
		r=0.004*(par()$usr[2]-par()$usr[1]))
	circles(x=lastPP$xx[1:tree$Nnode+Ntip(tree)],
		y=lastPP$yy[1:tree$Nnode+Ntip(tree)],
		r=0.007*(par()$usr[2]-par()$usr[1]))
	check<-textbox(x=c(par()$usr[1],par()$usr[1]+
		0.1*(par()$usr[2]-par()$usr[1])),
		y=par()$usr[4],c("click to stop"),justify="c")
	dev.flush()
	x<-unlist(locator(1))
	y<-x[2]
	x<-x[1]
	d<-sqrt((x-lastPP$xx)^2+(y-lastPP$yy)^2)
	nn<-which(d==min(d,na.rm=TRUE))
	## collapse tree & replot:
	while(!(x>par()$usr[1]&&x<par()$usr[1]+0.1*(par()$usr[2]-par()$usr[1])&&
		y>check&&y<par()$usr[4])){
		obj<-list(tree)
		if(nn>(Ntip(tree)+1)){
			obj<-splitTree(tree,list(node=nn,
				bp=tree$edge.length[which(tree$edge[,2]==nn)]))
			obj[[1]]$tip.label[which(obj[[1]]$tip.label=="NA")]<-
				tree$node.label[nn-Ntip(tree)]
			tips<-which(tree$tip.label%in%obj[[1]]$tip.label)
			theta<-atan(lastPP$yy[nn]/lastPP$xx[nn])
			if(lastPP$yy[nn]>0&&lastPP$xx[nn]<0) theta<-pi+theta
			else if(lastPP$yy[nn]<0&&lastPP$xx[nn]<0) theta<-pi+theta
			else if(lastPP$yy[nn]<0&&lastPP$xx[nn]>0) theta<-2*pi+theta
			ii<-which((c(tips,Ntip(tree)+1)-c(0,tips))>1)
			if(ii>1&&ii<=length(tips))
				tips<-c(tips[1:(ii-1)],theta/(2*pi)*Ntip(tree),tips[ii:length(tips)])
			else if(ii==1) tips<-c(theta/(2*pi)*Ntip(tree),tips)
			else if(ii>length(tips)) tips<-c(tips,theta/(2*pi)*Ntip(tree))
			tree<-read.tree(text=write.tree(obj[[1]]))
			M<-matrix(NA,min(c(max(4,Ntip(obj[[2]])),10)),length(tips))
			for(i in 1:ncol(M)) M[,i]<-seq(from=tips[i],to=i,length.out=nrow(M))
			colnames(M)<-tree$tip.label
			maxY<-seq(from=sum(sapply(obj,Ntip))-length(obj)+1,to=Ntip(tree),length.out=nrow(M))
			pw<-reorder(tree,"pruningwise")
			H<-nodeHeights(tree)
			for(i in 1:nrow(M)){
				if(hold) dev.hold()
				fan(tree,pw,H,xlim=lastPP$x.lim,ylim=lastPP$y.lim,
					tips=M[i,],maxY=maxY[i],...)
				if(nodes||i==nrow(M)){
					lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
					circles(x=lastPP$xx[1:Ntip(tree)],y=lastPP$yy[1:Ntip(tree)],
						r=0.004*(par()$usr[2]-par()$usr[1]))
					circles(x=lastPP$xx[1:tree$Nnode+Ntip(tree)],
						y=lastPP$yy[1:tree$Nnode+Ntip(tree)],
						r=0.007*(par()$usr[2]-par()$usr[1]))
				}
				dev.flush()
			}
		} else if(nn<=Ntip(tree)) {
			if(tree$tip.label[nn]%in%otree$node.label){
				on<-which(otree$node.label==tree$tip.label[nn])+Ntip(otree)
				obj<-splitTree(otree,list(node=on,
					bp=otree$edge.length[which(otree$edge[,2]==on)]))
				nlabel<-tree$tip.label[nn]
				tree$tip.label[nn]<-"NA"
				if(nn==1) tips<-c(rep(nn,Ntip(obj[[2]])),(nn+1):Ntip(tree))
				else if(nn>1&&nn<Ntip(tree)) tips<-c(1:(nn-1),rep(nn,Ntip(obj[[2]])),(nn+1):Ntip(tree))
				else if(nn==Ntip(tree)) tips<-c(1:(nn-1),rep(nn,Ntip(obj[[2]])))
				tree<-read.tree(text=write.tree(paste.tree(tree,obj[[2]])))
				M<-matrix(NA,min(c(max(4,Ntip(obj[[2]])),10)),length(tips))
				for(i in 1:ncol(M)) M[,i]<-seq(from=tips[i],to=i,length.out=nrow(M))
				colnames(M)<-tree$tip.label
				pw<-reorder(tree,"pruningwise")
				H<-nodeHeights(tree)
				for(i in 1:nrow(M)){
					if(hold) dev.hold()
					fan(tree,pw,H,xlim=lastPP$x.lim,ylim=lastPP$y.lim,
						tips=M[i,],maxY=NULL,...)
					if(nodes||i==nrow(M)){
						lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
						circles(x=lastPP$xx[1:Ntip(tree)],y=lastPP$yy[1:Ntip(tree)],
							r=0.004*(par()$usr[2]-par()$usr[1]))
						circles(x=lastPP$xx[1:tree$Nnode+Ntip(tree)],
							y=lastPP$yy[1:tree$Nnode+Ntip(tree)],
							r=0.007*(par()$usr[2]-par()$usr[1]))
					}
					dev.flush()
				}
			}
		} else { 
			cat("Collapsing to the root is not permitted. Choose another node.\n")
			flush.console()
		}
		check<-textbox(x=c(par()$usr[1],par()$usr[1]+
			0.1*(par()$usr[2]-par()$usr[1])),y=par()$usr[4],
			c("click to stop"),justify="c")
		x<-unlist(locator(1))
		y<-x[2]
		x<-x[1]
		lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		d<-sqrt((x-lastPP$xx)^2+(y-lastPP$yy)^2)
		nn<-which(d==min(d,na.rm=TRUE))
	}
	## turn locator bell back on
	options(locatorBell=TRUE)
	if(drop.extinct){ 
		if(!is.ultrametric(otree)) cat("Input tree was not ultrametric. Ignoring argument drop.extinct.\n")
		else { 
			th<-setNames(sapply(1:Ntip(tree),nodeheight,tree=tree),tree$tip.label)
			tips<-names(th)[which(th<(mean(sapply(1:Ntip(otree),
				nodeheight,tree=otree))-.Machine$double.eps^0.5))]
			tree<-drop.tip(tree,tips)
		}
	}
	tree
}

circles<-function(x,y,r,col="blue")
	nulo<-mapply(draw.circle,x=x,y=y,radius=r,MoreArgs=list(border=col,
		col="white",nv=20))

# simplified function to plot tree in type "fan"
# written by Liam J. Revell 2015
fan<-function(tree,pw=NULL,nH=NULL,colors="black",fsize=0.6,ftype="i",lwd=1,mar=rep(0.1,4),
	add=FALSE,part=1,xlim=NULL,ylim=NULL,tips=NULL,maxY=NULL,...){
	## set colors
	if(length(colors)==1) colors<-rep(colors,nrow(tree$edge))
	## set font
	ftype<-which(c("off","reg","b","i","bi")==ftype)-1
	# reorder
	cw<-tree
	if(is.null(pw)) pw<-reorder(tree,"pruningwise")
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
	R<-if(is.null(nH)) nodeHeights(cw) else nH
	# now put into a circular coordinate system
	x<-R*cos(Y)
	y<-R*sin(Y)
	# optimize x & y limits
	par(mar=mar)
	offsetFudge<-1.37 # empirically determined
	offset<-0
	if(is.null(xlim)||is.null(ylim)){
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
		xlim<-x.lim
		ylim<-y.lim
	}
	# plot tree
	if(!add) plot.new()
	plot.window(xlim=xlim,ylim=ylim,asp=1)
	# plot radial lines (edges)
	segments(x[,1],y[,1],x[,2],y[,2],col=colors,lwd=lwd,lend=2)
	# plot circular lines
	ii<-sapply(1:m+n,function(x,y) which(y==x),y=cw$edge)
	r<-sapply(1:m+n,function(x,y,R) R[match(x,y)],y=cw$edge,R=R)
	a1<-sapply(ii,function(x,Y) min(Y[x]),Y=Y)
	a2<-sapply(ii,function(x,Y) max(Y[x]),Y=Y)
	draw.arc(rep(0,tree$Nnode),rep(0,tree$Nnode),r,a1,a2,lwd=lwd,col=colors)
	# plot labels
	cw$tip.label<-gsub("_"," ",cw$tip.label)
	for(i in 1:n){
		ii<-which(cw$edge[,2]==i)
		aa<-Y[ii,2]/(2*pi)*360
		adj<-if(aa>90&&aa<270) c(1,0.25) else c(0,0.25)
		tt<-if(aa>90&&aa<270) paste(cw$tip.label[i]," ",sep="") else paste(" ",
			cw$tip.label[i],sep="")
		aa<-if(aa>90&&aa<270) 180+aa else aa
		if(ftype) text(x[ii,2],y[ii,2],tt,srt=aa,adj=adj,cex=fsize,font=ftype)
	}
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
