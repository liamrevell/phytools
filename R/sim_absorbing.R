## simulate stochastic diffusion on a tree with reflective or absorbing boundaries

sim.reflective<-function(tree,x0=0,sig2=1,
	bounds=c(-Inf,Inf),nsteps=1000,...){
	args<-list(...)
	args$tree<-tree
	args$x0<-x0
	args$sig2<-sig2
	args$bounds<-bounds
	args$nsteps<-nsteps
	if(is.null(args$reflective)) args$reflective<-TRUE
	do.call(sim.absorbing,args)
}

sim.absorbing<-function(tree,x0=0,sig2=1,
	bounds=c(-Inf,Inf),nsteps=1000,...){
	if(hasArg(plot)) plot<-list(...)$plot
	else plot=TRUE
	if(hasArg(reflective)) reflective<-list(...)$reflective
	else reflective<-FALSE
	if(!reflective){
		if(hasArg(semithreshold)) semithreshold<-list(...)$semithreshold
		else semithreshold<-FALSE
	} else semithreshold<-FALSE
	tree<-reorder(tree,"cladewise")
	ll<-max(nodeHeights(tree))
	tt<-map.to.singleton(make.era.map(tree,
		limits=seq(0,ll,length.out=nsteps+1)))
	delta_x<-rnorm(n=nrow(tt$edge),sd=sqrt(sig2*tt$edge.length))
	X<-matrix(NA,nrow(tt$edge),2)
	ROOT<-Ntip(tt)+1
	X[which(tt$edge[,1]==ROOT),1]<-x0
	for(i in 1:nrow(tt$edge)){
		X[i,2]<-X[i,1]+delta_x[i]
		if(!reflective){
			if(!semithreshold){
				if(X[i,2]<=bounds[1]) X[i,2]<--Inf
				else if(X[i,2]>=bounds[2]) X[i,2]<-Inf
			}
		} else {
			while(X[i,2]<=bounds[1]||X[i,2]>=bounds[2]){
				if(X[i,2]<=bounds[1]) X[i,2]<-2*bounds[1]-X[i,2]
				if(X[i,2]>=bounds[2]) X[i,2]<-2*bounds[2]-X[i,2]
			}
		}
		jj<-which(tt$edge[,1]==tt$edge[i,2])
		if(length(jj)>0) X[jj,1]<-X[i,2]
	}
	ii<-which(X==-Inf)
	if(length(ii)>0)	X[ii]<-bounds[1]
	ii<-which(X==Inf)
	if(length(ii)>0) X[ii]<-bounds[2]
	if(plot){
		if(!semithreshold){
			LIMS<-vector()
			LIMS[1]<-if(bounds[1]==-Inf) min(X) else bounds[1]
			LIMS[2]<-if(bounds[2]==Inf) max(X) else bounds[2]
		} else LIMS<-bounds
		hh<-hcl.colors(n=201)
		if(!semithreshold){
			edge_col<-hh[floor(200*((X[,1]-LIMS[1])/diff(LIMS)))+1]
			lwd<-rep(3,length(edge_col))
		} else {
			ii<-floor(200*((X[,1]-LIMS[1])/diff(LIMS)))+1
			lwd<-rep(3,length(ii))
			lwd[ii<1]<-1
			lwd[ii>201]<-1
			ii[ii<1]<-1
			ii[ii>201]<-201
			edge_col<-hh[ii]
		}	
		par(mfrow=c(2,1))
		par(mar=c(0.1,4.1,1.1,1.1))
		plot(tt,edge.color=edge_col,show.tip.label=FALSE,
			edge.width=3,y.lim=c(-0.15*Ntip(tt),Ntip(tt)))
		add.color.bar(0.5*max(nodeHeights(tree)),
			cols=hcl.colors(n=100),
			title="phenotype",lims=LIMS,digits=2,
			prompt=FALSE,lwd=10,outline=FALSE,
			x=0,y=-0.075*Ntip(tree),fsize=0.8,
			subtitle="")
		par(mar=c(5.1,4.1,1.1,1.1))
		plot(NA,xlim=c(0,ll),ylim=range(X),
			xlab="time",ylab="phenotype",bty="n",
			las=1)
		H<-nodeHeights(tt)
		for(i in 1:nrow(tt$edge))
			lines(H[i,],X[i,],col=edge_col[i],lwd=lwd[i])
		bcols<-hcl.colors(2)
		if(is.finite(bounds[1])){
			lines(c(0,max(nodeHeights(tree))),
				rep(bounds[1],2),lwd=8,
				col=make.transparent(bcols[1],0.2))
			text(0.1*max(nodeHeights(tree)),bounds[1],
				"lower bound",pos=3)
		}
		if(is.finite(bounds[2])){
			lines(c(0,max(nodeHeights(tree))),
				rep(bounds[2],2),lwd=8,
				col=make.transparent(bcols[2],0.2))
			text(0.1*max(nodeHeights(tree)),bounds[2],
				"upper bound",pos=1)
		}
	}
	if(semithreshold){
		X[X<bounds[1]]<-bounds[1]
		X[X>bounds[2]]<-bounds[2]
	}
	setNames(sapply(X=1:Ntip(tt),
		FUN=function(i,x,e) X[which(e[,2]==i),2],
		x=X,e=tt$edge),tt$tip.label)
}