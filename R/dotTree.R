## function to plot a tree with dots/circles for a plotted phenotype
## written by Liam J. Revell 2016, 2017, 2018, 2023

dotTree<-function(tree,x,legend=TRUE,method="plotTree",standardize=FALSE,...){
	if(is.data.frame(x)) x<-as.matrix(x)
	if(hasArg(data.type)) data.type<-list(...)$data.type
	else {
		## try to detect type
		if(is.numeric(x)) data.type<-"continuous"
		else data.type<-"discrete"
	}
	if(data.type=="continuous"){ 
		if(hasArg(colors)) colors<-list(...)$colors
		else colors<-"blue"
		dotTree.continuous(tree,x,colors[1],legend,method,standardize,...)
	} else if(data.type=="discrete"){ 
		if(hasArg(colors)) colors<-list(...)$colors
		else { 
			ss<-unique(as.vector(x))
			colors<-setNames(palette()[1:length(ss)],ss)
		}
		dotTree.discrete(tree,x,colors,legend,method,...)
	}
}
	
dotTree.continuous<-function(tree,x,color,legend,method,standardize,...){
	if(hasArg(border)) border<-list(...)$border
	else border<-par()$fg
	if(hasArg(cex.dot)) cex.dot<-list(...)$cex.dot
	else cex.dot<-1.0
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)&&method=="plotTree"){
		if(ncol(x)>1) method<-"phylogram"
		else x<-x[,1]
	}
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-c(1,0.8)
	if(length(fsize)==1) fsize<-rep(fsize,2)
	if(hasArg(x.space)) x.space<-list(...)$x.space
	else x.space<-0.1
	if(hasArg(k)) k<-list(...)$k
	else k<-0.8
	## reorder tree
	tree<-reorder(tree,"cladewise")
	## if standardize==TRUE
	if(standardize){
		if(is.matrix(x)){
			sd<-apply(x,2,function(x) sqrt(var(x)))
			x<-(x-matrix(rep(1,Ntip(tree)),Ntip(tree),1)%*%colMeans(x))/
				(matrix(rep(1,Ntip(tree)),Ntip(tree),1)%*%sd)
		} else if(is.vector(x)) x<-(x-mean(x))/sqrt(var(x))
	}
	## in case any x<0
	min.x<-min(x)
	max.x<-max(x)
	if(any(x<0)) x<-x-min(x)
	if(method=="plotTree"){
		fsize<-fsize[1]
		x<-x[tree$tip.label]
		## plot tree
		plotTree(tree,offset=1.7,ylim=c(-Ntip(tree)/25,
			Ntip(tree)),...)
		## get last phylo plot parameters
		obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		x.tip<-obj$xx[1:obj$Ntip]
		y.tip<-obj$yy[1:obj$Ntip]
		## plot points
		rr<-(k*x/max(x)+x.space)/2*diff(par()$usr[1:2])/
			diff(par()$usr[3:4])
		if(k<=0.8&&any(rr>(strwidth("W")*fsize/2))) 
			rr<-rr/max(rr)*strwidth("W")*fsize/2
		nulo<-mapply(draw.circle,x=x.tip+1.2*strwidth("W"),y=y.tip,
			radius=cex.dot*rr,MoreArgs=list(nv=200,col=color,border=border))
		## add legend
		if(legend){ 
			h<-dot.legend(x=par()$usr[1]+0.1*max(nodeHeights(tree)),
				y=0,min.x,max.x,Ntip=Ntip(tree),
				method="plotTree",...)
			if(standardize) text(h,0.1*(1+par()$usr[3]),"(SD units)",pos=4)
		}
	} else if(method=="phylogram"){
		if(is.vector(x)) x<-as.matrix(x)
		x[]<-x[tree$tip.label,]
		if(hasArg(mar)) mar<-list(...)$mar
		else mar<-rep(0.1,4)
		if(hasArg(xlim)) xlim<-list(...)$xlim
		else xlim<-c(-0.5,0.55+x.space*ncol(x)+x.space/2)
		if(hasArg(labels)) labels<-list(...)$labels
		else labels<-FALSE
		if(hasArg(ylim)) ylim<-list(...)$ylim
		else ylim<-c(if(legend) -0.1 else 0,if(labels) 1.1 else 1)
		## plot tree
		plot.new()
		par(mar=mar)
		plot.window(xlim=xlim,ylim=ylim)
		h<-phylogram(tree,...)
		## get last phylo plot parameters
		obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		x.tip<-rep(h,obj$Ntip)
		y.tip<-obj$yy[1:obj$Ntip]
		## plot points
		rr<-(k*x/max(x)+x.space)/2*diff(par()$usr[1:2])/diff(par()$usr[3:4])/
			(Ntip(tree)-1)
		if(k<=0.8&&any(rr>(strwidth("W")*fsize[1]/2))) 
			rr<-rr/max(rr)*strwidth("W")*fsize[1]/2
		for(i in 1:ncol(x)){
			nulo<-mapply(draw.circle,x=x.tip+1.2*strwidth("W")+x.space*(i-1),
				y=y.tip,radius=cex.dot*rr[,i],MoreArgs=list(nv=200,col=color,
				border=border))
		}
		## add legend
		if(legend){ 
			h<-dot.legend(x=-0.45,y=-0.04,min.x,max.x,Ntip=Ntip(tree),
				method="phylogram",...)
			if(standardize) text(h,-0.04,"(SD units)",pos=4)
		}
		if(labels){
			text(x=seq(max(x.tip)+1.2*strwidth("W"),
				max(x.tip)+1.2*strwidth("W")+x.space*(ncol(x)-1),by=x.space),
				y=rep(1.02,ncol(x)),colnames(x),srt=70,adj=c(0,0.5),cex=fsize[2])
		}
	}
}

dotTree.discrete<-function(tree,x,color,legend,method,...){
	if(hasArg(border)) border<-list(...)$border
	else border<-par()$fg
	if(hasArg(cex.dot)) cex.dot<-list(...)$cex.dot
	else cex.dot<-1.0
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)&&method=="plotTree"){
		if(ncol(x)>1) method<-"phylogram"
		else x<-x[,1]
	}
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-c(1,0.8)
	if(length(fsize)==1) fsize<-rep(fsize,2)
	if(hasArg(x.space)) x.space<-list(...)$x.space
	else x.space<-0.1
	## reorder tree
	tree<-reorder(tree,"cladewise")
	if(method=="plotTree"){
		x<-x[tree$tip.label]
		## plot tree
		plotTree(tree,offset=1.7,ylim=c(-1/25*Ntip(tree),
			Ntip(tree)),...)
		## get last phylo plot parameters
		obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		x.tip<-obj$xx[1:obj$Ntip]
		y.tip<-obj$yy[1:obj$Ntip]
		## plot points
		r<-min(0.8/2*diff(par()$usr[1:2])/diff(par()$usr[3:4]),
			strwidth("W")*fsize/2)
		nulo<-mapply(draw.circle,x=x.tip+1.2*strwidth("W"),y=y.tip,
			col=color[as.character(x)],MoreArgs=list(nv=200,
			radius=cex.dot*r,border=border))
		if(legend){ 
			add.simmap.legend(colors=color,prompt=FALSE,
				vertical=FALSE,shape="circle",
				x=par()$usr[1]+0.1*max(nodeHeights(tree)),
				y=-1/25*Ntip(tree),border=border)
		}
	} else if(method=="phylogram"){
		if(is.vector(x)) x<-as.matrix(x)
		x[]<-x[tree$tip.label,]
		if(hasArg(mar)) mar<-list(...)$mar
		else mar<-rep(0.1,4)
		if(hasArg(xlim)) xlim<-list(...)$xlim
		else xlim<-c(-0.5,0.55+x.space*ncol(x)+x.space/2)
		if(hasArg(labels)) labels<-list(...)$labels
		else labels<-FALSE
		if(hasArg(ylim)) ylim<-list(...)$ylim
		else ylim<-c(if(legend) -0.1 else 0,if(labels) 1.1 else 1)
		## plot tree
		plot.new()
		par(mar=mar)
		plot.window(xlim=xlim,ylim=ylim)
		h<-phylogram(tree,...)
		## get last phylo plot parameters
		obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		x.tip<-rep(h,obj$Ntip)
		y.tip<-obj$yy[1:obj$Ntip]
		## plot points
		r<-min(0.8/2*diff(par()$usr[1:2])/diff(par()$usr[3:4])/(Ntip(tree)-1),
			strwidth("W")*fsize/2)
		for(i in 1:ncol(x)){
			nulo<-mapply(draw.circle,x=x.tip+1.2*strwidth("W")+x.space*(i-1),
				y=y.tip,col=color[as.character(x[,i])],MoreArgs=list(nv=20,
				radius=cex.dot*r,border=border))
		}
		## add legend
		if(legend){ 
			add.simmap.legend(colors=color,prompt=FALSE,
				vertical=FALSE,shape="circle",x=-0.45,y=-0.06,
				border=border)
		}
		if(labels){
			text(x=seq(max(x.tip)+1.2*strwidth("W"),
				max(x.tip)+1.2*strwidth("W")+x.space*(ncol(x)-1),by=x.space),
				y=rep(1.02,ncol(x)),colnames(x),srt=70,adj=c(0,0.5),cex=fsize[2])
		}
	}
}

## dot legend
## written by Liam J. Revell 2016, 2023

dot.legend<-function(x,y,min,max,Ntip,length=5,prompt=FALSE,
	method="plotTree",...){
	if(hasArg(border)) border<-list(...)$border
	else border<-par()$fg
	if(hasArg(cex.dot)) cex.dot<-list(...)$cex.dot
	else cex.dot<-1.0
	if(hasArg(cex)) cex<-list(...)$cex
	else cex<-1
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-1
	if(hasArg(colors)) colors<-list(...)$colors
	else colors<-"blue"
	if(hasArg(leg.space)) leg.space<-list(...)$leg.space
	else leg.space<-0.2
	if(hasArg(k)) k<-list(...)$k
	else k<-0.8
	if(prompt){
		obj<-locator(1)
		x<-obj$x
		y<-obj$y
	}
	if(method=="plotTree"){
		text(x,y-0.5*(Ntip/25),round(min,2),pos=1,cex=cex)
		s<-(k*max(min,0)/min(max,max+min)+0.1)/2*diff(par()$usr[1:2])/
			diff(par()$usr[3:4])
		e<-(k+0.1)/2*diff(par()$usr[1:2])/diff(par()$usr[3:4])
		rr<-seq(s,e,length.out=length)
		if(k<=0.8&&any(rr>(strwidth("W")*fsize/2))) 
			rr<-rr/max(rr)*strwidth("W")*fsize/2
		temp<-c(0,cumsum((1+leg.space)*rep(2*max(rr),length-1)))
		nulo<-mapply(draw.circle,x=x+temp,y=rep(y,length),radius=cex.dot*rr,
			MoreArgs=list(nv=200,col=colors,border=border))
		text(max(x+temp),y-0.5*(Ntip/25),round(max,2),pos=1,cex=cex)
		y1<-0.1/25*Ntip
		lines(c(x,max(x+temp)),rep(y-0.5*(Ntip/25)-y1,2))
		lines(c(x,x),y-c(y1+0.5*(Ntip/25),2*y1+0.5*(Ntip/25)))
		lines(c(max(x+temp),max(x+temp)),y-c(y1+0.5*(Ntip/25),2*y1+0.5*(Ntip/25)))
	} else if(method=="phylogram"){
		text(x,y-0.04,round(min,2),pos=1,cex=cex)		
		s<-(k*max(min,0)/(min(max,max+min))+0.1)/2*
			diff(par()$usr[1:2])/diff(par()$usr[3:4])/(Ntip-1)
		e<-(k+0.1)/2*diff(par()$usr[1:2])/diff(par()$usr[3:4])/(Ntip-1)
		rr<-seq(s,e,length.out=length)
		if(k<=0.8&&any(rr>(strwidth("W")*fsize/2))) 
			rr<-rr/max(rr)*strwidth("W")*fsize/2
		temp<-c(0,cumsum((1+leg.space)*rep(2*max(rr),length-1)))
		nulo<-mapply(draw.circle,x=x+temp,y=rep(y,length),radius=cex.dot*rr,
			MoreArgs=list(nv=200,col=colors,border=border))
		text(max(x+temp),y-0.04,round(max,2),pos=1,cex=cex)
		y1<-0.01
		lines(c(x,max(x+temp)),rep(y-0.02-y1,2))
		lines(c(x,x),y-c(y1+0.02,2*y1+0.02))
		lines(c(max(x+temp),max(x+temp)),y-c(y1+0.02,2*y1+0.02))
	}
	invisible(max(x+temp)+0.5*max(rr))
}