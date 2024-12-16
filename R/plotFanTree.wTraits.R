## plot fan or arc tree with trait data

plotFanTree.wTraits<-function(tree,X,type=c("arc","fan"),...){
	X<-if(is.vector(X)) as.matrix(X[tree$tip.label]) else 
		X[tree$tip.label,,drop=FALSE]
	h<-max(nodeHeights(tree))
	d<-min(ncol(X)*0.07*h,h)
	type<-type[1]
	if(!(type%in%c("arc","fan"))) type<-"fan"
	ftype<-if(hasArg(ftype)) list(...)$ftype else "i"
	fsize<-if(hasArg(fsize)) list(...)$fsize else 0.5
	part<-if(hasArg(part)) list(...)$part else 
		min(0.99,(Ntip(tree)-2)/Ntip(tree))
	arc_height<-if(hasArg(arc_height)) list(...)$arc_height else 
		0.7
	spacer<-if(hasArg(spacer)) list(...)$spacer else 0.025
	spacer<-spacer*(2*pi*part/(Ntip(tree)-1))/2
	xlim<-if(hasArg(xlim)) list(...)$xlim else NULL
	ylim<-if(hasArg(ylim)) list(...)$ylim else NULL
	if(hasArg(colors)) colors<-list(...)$colors
	else {
		colors<-list()
		for(i in 1:ncol(X)){
			if(is.numeric(X[,i])){
				colors[[i]]<-setNames(hcl.colors(n=100),1:100)
			} else {
				if(!is.factor(X[,i])) X[,i]<-as.factor(X[,i])
				colors[[i]]<-setNames(
					palette.colors(n=length(levels(X[,i]))),
					levels(X[,i]))
			}
		}
	}
	tt<-uu<-force.ultrametric(tree,method="extend",message=FALSE)
	tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]<-
		tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]+d
	plotTree(tt,type=type,ftype=ftype,fsize=fsize,
		part=part,color="transparent",
		arc_height=arc_height*h/max(nodeHeights(tt)),
		xlim=xlim,ylim=ylim)
	pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	outer_rad<-max(pp$xx)
	plotTree(uu,type=type,ftype="off",part=part,
		lwd=1,add=TRUE,xlim=pp$x.lim,ylim=pp$y.lim,
		arc_height=arc_height,ftype="off",color="#D3D3D3")
	pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
	inner_rad<-max(pp$xx)
	par(lty="solid")
	plotTree(tree,type=type,ftype="off",part=part,
		lwd=1,add=TRUE,xlim=pp$x.lim,ylim=pp$y.lim,
		arc_height=arc_height,ftype="off",color=par()$fg)
	par(lend=3)
	for(i in 1:ncol(X)){
		if(is.numeric(X[,i])){
			x_seq<-seq(min(X[,i]),max(X[,i]),length.out=100)
			x_ind<-sapply(X[,i],function(x,y) which.min((x-y)^2),
				y=x_seq)
			colors[[i]]<-colorRampPalette(colors[[i]])(n=100)
			cols<-colors[[i]][x_ind]
		} else {
			cols<-colors[[i]][X[tree$tip.label,i]]	
		}
		for(j in 1:Ntip(tree)){
			start<-if(pp$xx[j]>0) 
				(i-1)*(d/ncol(X))+(2/7)*(d/ncol(X)) else 
					-((i-1)*(d/ncol(X))+(2/7)*(d/ncol(X)))
			end<-if(pp$xx[j]>0) i*d/ncol(X) else -i*d/ncol(X)
			th<-atan(pp$yy[j]/pp$xx[j])
			theta<-(2*pi*part/(Ntip(tree)-1))/2-spacer
			sign<-if(pp$xx[j]>0) 1 else -1
			H1<-(sign*inner_rad+start)/cos(theta)
			H2<-(sign*inner_rad+end)/cos(theta)
			th_up<-th+theta
			th_down<-th-theta
			x<-c(H1*cos(th_down),H2*cos(th_down),
				H2*cos(th_up),H1*cos(th_up))
			y<-c(H1*sin(th_down),H2*sin(th_down),
				H2*sin(th_up),H1*sin(th_up))
			polygon(x,y,col=cols[j],border=FALSE)
		}
	}
	invisible(colors)
}