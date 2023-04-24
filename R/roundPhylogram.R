## function plots a sigmoid or spline phylogram or cladogram
## written by Liam J. Revell 2022, 2023

splinePhylogram<-function(tree,...){
	args<-list(...)
	args$tree<-tree
	args$plot<-FALSE
	args$direction<-"rightwards"
	do.call(plotTree,args)
	pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	if(hasArg(df)) df<-list(...)$df
	else df<-50
	if(hasArg(res)) res<-list(...)$res
	else res<-4*df
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-par()$lwd
	if(hasArg(color)) color<-list(...)$color
	else color<-par()$fg
	if(is.null(tree$edge.length)) tree<-compute.brlen(tree)
	h<-max(nodeHeights(tree))/res
	tree<-make.era.map(tree,seq(0,max(nodeHeights(tree)),by=h))
	tree<-map.to.singleton(tree)
	args<-list(...)
	args$tree<-tree
	args$color<-"transparent"
	args$direction<-"rightwards"
	args$add<-TRUE
	dev.hold()
	do.call(plotTree,args)
	obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	phy<-list(edge=obj$edge,Nnode=obj$Nnode,
		tip.label=1:obj$Ntip)
	attr(phy,"class")<-"phylo"
	for(i in 1:Ntip(phy)){
		aa<-c(i,Ancestors(phy,i))
		xx<-obj$xx[aa]
		yy<-obj$yy[aa]
		DF<-min(min(length(xx)-1,df),df)
		tmp<-predict(smooth.spline(xx,yy,df=DF))
		tmp$x<-c(xx[length(xx)],tmp$x)
		tmp$y<-c(yy[length(yy)],tmp$y)
		lines(tmp,lwd=lwd,col=color)
	}
	nulo<-dev.flush()
	args$tree<-collapse.singles(tree)
	args$add<-TRUE
	ffg<-par()$fg
	par(fg="transparent")
	do.call(plotTree,args)
	par(fg=ffg)
	assign("last_plot.phylo",pp,envir=.PlotPhyloEnv)
}

compute.brlen.simmap<-function(phy,method="Grafen",power=1,...){
	if(inherits(phy,"simmap")){
		tt<-as.phylo(phy)
		tt<-compute.brlen(tt,method=method,power=power,...)
		ss<-tt$edge.length/phy$edge.length
		phy$maps<-mapply("*",ss,phy$maps,SIMPLIFY=FALSE)
		phy$mapped.edge<-phy$mapped.edge*
			matrix(rep(ss,ncol(phy$mapped.edge)),
			nrow(phy$mapped.edge),ncol(phy$mapped.edge))
		phy$edge.length<-tt$edge.length
	} else if(inherits(phy,"phylo")) 
		phy<-compute.brlen(phy,method=method,power=power,...)
	phy
}		

sigmoidPhylogram<-function(tree,...){
	if(hasArg(outline)) outline<-list(...)$outline
	else outline<-FALSE
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-2
	if(outline){
		args<-list(...)
		args$tree<-as.phylo(tree)
		args$lwd<-lwd+2
		args$colors<-par()$fg
		args$outline<-FALSE
		do.call(sigmoidPhylogram,args)
	}
	## b=5,m1=0.01,m2=0.5,v=1
	b<-if(hasArg(b)) list(...)$b else 5
	m1<-if(hasArg(m1)) list(...)$m1 else 0.01
	m2<-if(hasArg(m2)) list(...)$m2 else 0.5
	v<-if(hasArg(v)) list(...)$v else 1
	direction<-if(hasArg(direction)) list(...)$direction else "rightwards"
	show.hidden<-if(hasArg(show.hidden)) list(...)$show.hidden else FALSE
	if(inherits(tree,"simmap")){
		if(hasArg(colors)) colors<-list(...)$colors
		else {
			ss<-sort(unique(c(getStates(tree,"nodes"),
				getStates(tree,"tips"))))
			mm<-length(ss)
			colors<-setNames(
				colorRampPalette(palette()[1:min(8,mm)])(mm),
				ss)
		}
	} else if(inherits(tree,"phylo")) {
		if(hasArg(color)){ 
			color<-setNames(list(...)$color,"1")
			colors<-color
		} else colors<-setNames(par()$fg,"1")
		tree<-paintSubTree(tree,Ntip(tree)+1,"1")
	}
	if(hasArg(res)) res<-list(...)$res
	else res<-199
	if(hasArg(use.edge.length)) 
		use.edge.length<-list(...)$use.edge.length
	else use.edge.length<-TRUE
	if(!use.edge.length){
		if(hasArg(power)) power<-list(...)$power
		else power<-1
		tree<-compute.brlen.simmap(tree,power=power)
	}
	h<-max(nodeHeights(tree))
	args<-list(...)
	args$power<-NULL
	args$res<-NULL
	args$colors<-NULL
	args$b<-NULL
	args$m1<-NULL
	args$m2<-NULL
	args$v<-NULL
	args$tree<-tree
	args$color<-if(show.hidden) make.transparent("red",0.25) else
		"transparent"
	if(outline) args$add<-TRUE
	args$outline<-FALSE
	dev.hold()
	par_fg<-par()$fg
	par(fg="transparent")
	do.call(plotTree,args)
	par(fg=par_fg)
	pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	par_usr<-par()$usr
	args$add<-TRUE
	if(direction=="downwards"){
		args$xlim<-pp$x.lim
		args$ylim<-pp$y.lim[2:1]-pp$y.lim[1]
		args$direction<-"upwards"
	} else if(direction=="leftwards"){
		args$xlim<-pp$x.lim[2:1]-pp$x.lim[1]
		args$ylim<-pp$y.lim
		args$direction<-"rightwards"
	}
	if(outline){
		par_fg<-par()$fg
		par(fg="transparent")
	}
	do.call(plotTree,args)
	if(outline) par(fg=par_fg)
	pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	xx<-if(direction%in%c("rightwards","leftwards")) pp$xx else pp$yy
	yy<-if(direction%in%c("rightwards","leftwards")) pp$yy else pp$xx
	## Yt<-A+(K-A)/((C+exp(-B*(t-M)))^(1/v))
	sigmoid<-function(x,.A=A,.K=K,.C=C,.B=B,.M=M,.v=v)
		return(.A+(.K-.A)/((.C+exp(-.B*(x-.M)))^(1/.v)))
	for(i in 1:nrow(tree$edge)){
		A<-yy[tree$edge[i,1]]
		K<-yy[tree$edge[i,2]]
		if(i==1) dy<-abs(A-K)
		B<-b*Ntip(tree)/h
		t<-seq(xx[tree$edge[i,1]],xx[tree$edge[i,2]],
			length.out=res)
		t<-sort(c(t,t[1]+cumsum(tree$maps[[i]])))
		dd<-diff(range(t))
		M<-t[1] + if(m1*h>(m2*dd)) m2*dd else m1*h
		C<-1
		Yt<-c(A,sigmoid(t),K)
		t<-c(t[1],t,t[length(t)])
		COLS<-vector()
		bb<-setNames(t[1]+cumsum(tree$maps[[i]]),names(tree$maps[[i]]))
		for(j in 1:length(t)) 
			COLS[j]<-colors[names(bb[which(t[j]<=bb)])[1]]
		nn<-length(t)
		if(direction%in%c("rightwards","leftwards"))
			segments(t[1:(nn-1)],Yt[1:(nn-1)],x1=t[2:nn],y1=Yt[2:nn],
				col=COLS,lwd=lwd)
		else if(direction%in%c("upwards","downwards"))
			segments(Yt[1:(nn-1)],t[1:(nn-1)],x1=Yt[2:nn],y1=t[2:nn],
				col=COLS,lwd=lwd)
	}
	nulo<-dev.flush()
	pp$edge<-tree$edge
	assign("last_plot.phylo",pp,envir=.PlotPhyloEnv)
}

## function plots a round phylogram
## written by Liam J. Revell 2014, 2015, 2016, 2023

roundPhylogram<-function(tree,fsize=1.0,ftype="reg",lwd=2,mar=NULL,offset=NULL,
	direction="rightwards",type="phylogram",xlim=NULL,ylim=NULL,...){
	if(inherits(tree,"multiPhylo")){
		par(ask=TRUE)
		tt<-lapply(tree,roundPhylogram,fsize=fsize,ftype=ftype,lwd=lwd,mar=mar,offset=offset, 
			direction=direction,type=type,xlim=xlim,ylim=ylim)
	} else {
		if(hasArg(lty)) lty<-list(...)$lty
		else lty<-"solid"
		if(length(lty)!=nrow(tree$edge)) lty<-rep(lty,ceiling(nrow(tree$edge)/length(lty)))
		if(length(lwd)!=nrow(tree$edge)) lwd<-rep(lwd,ceiling(nrow(tree$edge)/length(lwd)))
		if(type=="cladogram"||is.null(tree$edge.length)) tree<-compute.brlen(tree)
		ftype<-which(c("off","reg","b","i","bi")==ftype)-1
		if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
		# swap out "_" character for spaces (assumes _ is a place holder)
		tree$tip.label<-gsub("_"," ",tree$tip.label)
		if(is.null(mar)) mar=rep(0.1,4)
		n<-length(tree$tip.label)
		# set offset fudge (empirically determined)
		offsetFudge<-1.37
		# reorder cladewise to assign tip positions
		cw<-reorder(tree,"cladewise")
		io<-reorder(tree,index.only=TRUE)
		lwd<-lwd[io]
		lty<-lty[io]
		y<-vector(length=n+cw$Nnode)
		y[cw$edge[cw$edge[,2]<=n,2]]<-1:n
		# reorder pruningwise for post-order traversal
		pw<-reorder(tree,"pruningwise")
		nn<-unique(pw$edge[,1])
		# compute vertical position of each edge
		for(i in 1:length(nn)){
			yy<-y[pw$edge[which(pw$edge[,1]==nn[i]),2]]
			y[nn[i]]<-mean(range(yy))
		}
		# compute start & end points of each edge
		X<-nodeHeights(cw)
		if(is.null(xlim)){
			pp<-par("pin")[1]
			sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+offsetFudge*fsize*strwidth("W",units="inches")
			alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=X,sw=sw,pp=pp,interval=c(0,1e6))$minimum
			xlim<-c(min(X),max(X)+sw/alp)
		}
		if(is.null(ylim)) ylim=c(1,max(y))
		## end preliminaries
		# open & size a new plot
		plot.new()
		par(mar=mar)
		if(is.null(offset)) offset<-0.2*lwd/3+0.2/3
		plot.window(xlim=xlim,ylim=ylim)
		# plot edges
		for(i in 1:nrow(X)){
			x<-NULL
			b<-y[cw$edge[i,1]]
			c<-X[i,1]
			d<-if(y[cw$edge[i,2]]>y[cw$edge[i,1]]) 1 else -1
			xx<-X[i,2]
			yy<-y[cw$edge[i,2]]
			a<-(xx-c)/(yy-b)^2
			curve(d*sqrt((x-c)/a)+b,from=X[i,1],to=X[i,2],add=TRUE,lwd=lwd[i],lty=lty[i])
		}
		# plot tip labels
		for(i in 1:n)
			if(ftype) text(X[which(cw$edge[,2]==i),2],y[i],tree$tip.label[i],pos=4,offset=offset,font=ftype,cex=fsize)
		PP<-list(type=type,use.edge.length=if(type=="phylogram") TRUE else FALSE,node.pos=1,
			show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,font=ftype,cex=fsize,
			adj=0,srt=0,no.margin=FALSE,label.offset=offset,x.lim=par()$usr[1:2],y.lim=par()$usr[3:4],
			direction=direction,tip.color="black",Ntip=length(cw$tip.label),Nnode=cw$Nnode,edge=tree$edge,
			xx=sapply(1:(length(cw$tip.label)+cw$Nnode),function(x,y,z) y[match(x,z)],y=X,z=cw$edge),
			yy=y)
		assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
	}
}
