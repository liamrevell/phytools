## function to plot a tree with dots/circles for a plotted phenotype
## written by Liam J. Revell 2016

dotTree<-function(tree,x,legend=TRUE,method="plotTree",standardize=FALSE,...){
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
			ss<-sort(unique(x))
			colors<-setNames(palette()[1:length(ss)],ss)
		}
		dotTree.discrete(tree,x,colors,legend,method,...)
	}
}
	
dotTree.continuous<-function(tree,x,color,legend,method,standardize,...){	
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)&&method=="plotTree"){
		if(ncol(x)>1) method<-"phylogram"
		else x<-x[,1]
	}
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
		x<-x[tree$tip.label]
		## plot tree
		plotTree(tree,offset=1.7,ylim=c(-1/25*Ntip(tree),
			Ntip(tree)),...)
		## get last phylo plot parameters
		obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		x.tip<-obj$xx[1:obj$Ntip]
		y.tip<-obj$yy[1:obj$Ntip]
		## plot points
		draw.circle(x.tip+1.2*strwidth("W"),y.tip,nv=200,
			radius=(0.8*x/max(x)+0.1)/2*diff(par()$usr[1:2])/
			diff(par()$usr[3:4]),col=color)
		## add legend
		if(legend){ 
			h<-dot.legend(x=par()$usr[1]+0.1*max(nodeHeights(tree)),
				y=0.1*(1+par()$usr[3]),min.x,max.x,Ntip=Ntip(tree),
				method="plotTree",...)
			if(standardize) text(h,0.1*(1+par()$usr[3]),"(SD units)",pos=4)
		}
	} else if(method=="phylogram"){
		if(is.vector(x)) x<-as.matrix(x)
		x<-x[tree$tip.label,]
		## plot tree
		plot.new()
		par(mar=rep(0.1,4))
		plot.window(xlim=c(-0.5,0.55+0.1*ncol(x)),
			ylim=c(if(legend) -0.1 else 0,1))
		h<-phylogram(tree,...)
		## get last phylo plot parameters
		obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		x.tip<-rep(h,obj$Ntip)
		y.tip<-obj$yy[1:obj$Ntip]
		## plot points
		for(i in 1:ncol(x)){
			draw.circle(x.tip+1.2*strwidth("W")+0.1*(i-1),y.tip,
				nv=200,radius=(0.8*x[,i]/max(x)+0.1)/Ntip(tree)/
					2*diff(par()$usr[1:2])/diff(par()$usr[3:4]),col=color)
		}
		## add legend
		if(legend){ 
			h<-dot.legend(x=-0.45,y=-0.04,min.x,max.x,Ntip=Ntip(tree),
				method="phylogram",...)
			if(standardize) text(h,-0.04,"(SD units)",pos=4)
		}
	}
}

dotTree.discrete<-function(tree,x,color,legend,method,...){
	if(is.data.frame(x)) x<-as.matrix(x)
	if(is.matrix(x)&&method=="plotTree"){
		if(ncol(x)>1) method<-"phylogram"
		else x<-x[,1]
	}
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
		draw.circle(x.tip+1.2*strwidth("W"),y.tip,nv=200,
			radius=0.8/2*diff(par()$usr[1:2])/diff(par()$usr[3:4]),
			col=color[as.character(x)])
		if(legend){ 
			add.simmap.legend(colors=color,prompt=FALSE,
				vertical=FALSE,shape="circle",
				x=par()$usr[1]+0.1*max(nodeHeights(tree)),
				y=0.1*(1+par()$usr[3]))
		}
	} else if(method=="phylogram"){
		if(is.vector(x)) x<-as.matrix(x)
		x<-x[tree$tip.label,]
		## plot tree
		plot.new()
		par(mar=rep(0.1,4))
		plot.window(xlim=c(-0.5,0.55+0.1*ncol(x)),
			ylim=c(if(legend) -0.06 else 0,1))
		h<-phylogram(tree,...)
		## get last phylo plot parameters
		obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
		x.tip<-rep(h,obj$Ntip)
		y.tip<-obj$yy[1:obj$Ntip]
		## plot points
		for(i in 1:ncol(x)){
			draw.circle(x.tip+1.2*strwidth("W")+0.1*(i-1),y.tip,
				nv=200,radius=0.8/Ntip(tree)/2*diff(par()$usr[1:2])/diff(par()$usr[3:4]),
				col=color[as.character(x[,i])])
		}
		## add legend
		if(legend){ 
			add.simmap.legend(colors=color,prompt=FALSE,
				vertical=FALSE,shape="circle",x=-0.45,y=-0.06)
		}
	}
}

## dot legend
## written by Liam J. Revell 2016

dot.legend<-function(x,y,min,max,Ntip,length=5,prompt=FALSE,method="plotTree",...){
	if(hasArg(cex)) cex<-list(...)$cex
	else cex<-1
	if(hasArg(colors)) colors<-list(...)$colors
	else colors<-"blue"
	if(prompt){
		obj<-locator(1)
		x<-obj$x
		y<-obj$y
	}
	if(method=="plotTree"){
		text(x,y-0.5,round(min,2),pos=1,cex=cex)
		s<-(0.8*max(min,0)/min(max,max+min)+0.1)/2*diff(par()$usr[1:2])/
			diff(par()$usr[3:4])
		e<-(0.8/2*diff(par()$usr[1:2])+0.1)/diff(par()$usr[3:4])
		rr<-seq(s,e,length.out=length)
		temp<-c(0,cumsum(1.1*rep(2*max(rr),length-1)))
		draw.circle(x+temp,rep(y,length),nv=200,radius=rr,col=colors)
		text(max(x+temp),y-0.5,round(max,2),pos=1,cex=cex)
		y1<-0.1/25*Ntip
		lines(c(x,max(x+temp)),rep(y-0.5-y1,2))
		lines(c(x,x),y-c(y1+0.5,2*y1+0.5))
		lines(c(max(x+temp),max(x+temp)),y-c(y1+0.5,2*y1+0.5))
	} else if(method=="phylogram"){
		text(x,y-0.04,round(min,2),pos=1,cex=cex)
		s<-(0.8*max(min,0)/(min(max,max+min))+0.1)/(2*Ntip)*
			diff(par()$usr[1:2])/diff(par()$usr[3:4])
		e<-(0.8*diff(par()$usr[1:2])+0.1)/(2*Ntip*diff(par()$usr[3:4]))
		rr<-seq(s,e,length.out=length)
		temp<-c(0,cumsum(1.1*rep(2*max(rr),length-1)))
		draw.circle(x+temp,rep(y,length),nv=200,radius=rr,col=colors)
		text(max(x+temp),y-0.04,round(max,2),pos=1,cex=cex)
		y1<-0.01
		lines(c(x,max(x+temp)),rep(y-0.02-y1,2))
		lines(c(x,x),y-c(y1+0.02,2*y1+0.02))
		lines(c(max(x+temp),max(x+temp)),y-c(y1+0.02,2*y1+0.02))
	}
	invisible(max(x+temp)+0.5*max(rr))
}