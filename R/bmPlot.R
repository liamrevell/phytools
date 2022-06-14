## visualize discrete time Brownian simulation on a tree
## written by Liam J. Revell 2012, 2013, 2015, 2020, 2022

bmPlot<-function(tree,type="BM",anc=0,sig2=1/1000,ngen=1000,...){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(hasArg(return.tree)) return.tree<-list(...)$return.tree
	else return.tree<-TRUE

	tr<-reorder(tree)
	H<-nodeHeights(tr)
	tr$edge.length<-round(tr$edge.length/max(H)*ngen)
	tr<-di2multi(tr)
	h<-nodeHeights(tr)[,1]

	bmSim<-function(start,n)
		cumsum(c(start,rnorm(n,sd=sqrt(sig2))))

	X<-T<-list()
	N<-length(tr$tip)
	for(i in 1:nrow(tr$edge)){
		if(tr$edge[i,1]==(N+1)) X[[i]]<-bmSim(anc,tr$edge.length[i])
		else {
			parent<-match(tr$edge[i,1],tr$edge[,2])
			X[[i]]<-bmSim(X[[parent]][length(X[[parent]])],tr$edge.length[i])
		}
		T[[i]]<-h[i]+0:tr$edge.length[i]
	}
	
	if(type=="BM") cols<-bm(X,T,...)
	else if(type=="threshold") cols<-th(X,T,...)

	if(return.tree){
		if(type=="BM") tr$X<-X
		else if(type=="threshold"){
			if(hasArg(thresholds)) thresholds<-list(...)$thresholds
			else stop("no thresholds provided for type=\"threshold\"")
			thresholds<-setNames(c(thresholds,Inf),letters[1:(length(thresholds)+1)])
			xx<-lapply(X,function(x,y) sapply(x[2:length(x)],threshState,thresholds=y),y=thresholds)
			maps<-lapply(tr$edge.length,function(x) rep(1,x))	
			maps<-mapply(function(x,y) setNames(x,y),x=maps,y=xx)
			maps<-lapply(maps,mergeAdjacent)
			tr$maps<-maps
			class(tr)<-c("simmap",setdiff(class(tr),"simmap"))
		}
	}

	Y<-matrix(NA,nrow(tr$edge),2)
	Y[,1]<-sapply(X,function(x) x[1])
	Y[,2]<-sapply(X,function(x) x[length(x)])

	x<-Y[tr$edge[,2]%in%1:N,2]
	names(x)<-tr$tip.label[tr$edge[tr$edge[,2]%in%1:N,2]]
	a<-c(anc,Y[tr$edge[,2]%in%(N+2:tr$Nnode),2])
	names(a)<-c(N+1,tr$edge[tr$edge[,2]%in%(N+2:tr$Nnode),2])
	
	xx<-c(x[tr$tip],a[as.character(N+1:tr$Nnode)])
	if(!return.tree) return(xx)
	else return(list(x=xx,tree=tr,colors=cols))
}

# plots type="BM"
bm<-function(X,T,...){
	if(hasArg(bty)) bty<-list(...)$bty
	else bty<-"o"
	if(hasArg(lwd)) lwd<-list(...)$lwd
	else lwd<-1
	minX<-min(sapply(X,min))
	maxX<-max(sapply(X,max))
	plot(X[[1]][1],T[[1]][1],ylim=c(0,max(sapply(T,max))),
		xlim=c(minX,maxX),ylab="time",xlab="phenotype",
		bty=bty)
	for(i in 1:length(X)) lines(X[[i]],T[[i]],lwd=lwd)
	if(hasArg(colors)) cols<-list(...)$colors
	else cols<-"black"
	return(cols)
}

# plots type="threshold"
th<-function(X,T,...){
	if(hasArg(bty)) bty<-list(...)$bty
	else bty<-"o"
	minX<-min(sapply(X,min))
	maxX<-max(sapply(X,max))
	if(hasArg(thresholds)) thresholds<-list(...)$thresholds
	else stop("no thresholds provided for type='threshold'")
	if(hasArg(colors)) cols<-list(...)$colors
	else {
		cols<-c("black","red","blue","green")
		while(length(thresholds)>(length(cols)-1)) cols<-c(cols,c("black","red","blue","green"))
	}
	thresholds<-thresholds+1e-16
	plot(X[[1]][1],T[[1]][1],ylim=c(0,max(sapply(T,max))),xlim=c(minX,maxX),ylab="time",xlab="liability",
		bty=bty)
	for(i in 1:length(X))
		if(length(X[[i]])>1)
			for(j in 1:(length(X[[i]])-1))
				if(threshCol(X[[i]][j],thresholds,cols)==threshCol(X[[i]][j+1],thresholds,cols))
					lines(X[[i]][c(j,j+1)],T[[i]][c(j,j+1)],
						col=threshCol(X[[i]][j+1],thresholds,cols),lwd=lwd)
				else {
					t<-thresholds[which(sort(c(thresholds,min(X[[i]][c(j,j+1)])))==min(X[[i]][c(j,j+1)]))]
					lines(c(X[[i]][j],t),c(T[[i]][j],
						T[[i]][j]+abs(diff(c(X[[i]][j],t))/diff(X[[i]][c(j,j+1)]))),
						col=threshCol(X[[i]][j],thresholds,cols),lwd=lwd)
					lines(c(X[[i]][j+1],t),c(T[[i]][j+1],
						T[[i]][j+1]-abs(diff(c(X[[i]][j+1],t))/diff(X[[i]][c(j,j+1)]))),
						col=threshCol(X[[i]][j+1],thresholds,cols),lwd=lwd)
				}
	for(i in 1:length(thresholds)) lines(rep(thresholds[i],2),c(0,max(sapply(T,max))),lty="dashed")
	return(setNames(cols[1:(length(thresholds)+1)],letters[1:(length(thresholds)+1)]))
}

# returns a color based on position relative to thresholds
threshCol<-function(x,thresholds,cols){
	t<-c(-Inf,thresholds,Inf); i<-1
	while(x>t[i]) i<-i+1
	return(cols[i-1])
}

# merge adjacent map in the same state
mergeAdjacent<-function(x){
	ii<-1
	y<-x[1]
	if(length(x)>1){
		for(i in 2:length(x)){
			if(names(x)[i]==names(y)[ii]) y[ii]<-y[ii]+x[i]
			else{
				y<-c(y,x[i])
				ii<-ii+1
			}
		}
	}
	return(y)
}

