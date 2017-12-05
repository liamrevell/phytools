## computing the mean number of character changes through time from a set of stochastic map trees
## written by Liam J. Revell 2017

ctt<-function(trees,segments=20,...){
	if(!(inherits(trees,"multiSimmap")))
		stop("trees should be an object of class \"multiSimmap\".")
	tree<-as.phylo(trees[[1]])
	changes<-sapply(trees,getChanges)
	h<-max(nodeHeights(tree))
	b<-segments
	segs<-cbind(seq(0,h-h/b,h/b),
		seq(1/b*h,h,h/b))
	nchanges<-rep(0,b)
	for(i in 1:length(changes)){
		for(j in 1:length(changes[[i]])){
			ind<-which((changes[[i]][j]>segs[,1])+
				(changes[[i]][j]<=segs[,2])==2)
			nchanges[ind]<-nchanges[ind]+1/length(changes)
		}
	}
	LTT<-ltt(tree,plot=FALSE)
	LTT<-cbind(LTT$ltt[2:(length(LTT$ltt)-1)],
		LTT$times[2:(length(LTT$ltt)-1)],
		LTT$times[3:length(LTT$ltt)])
	ii<-1
	edge.length<-rep(0,b)
	for(i in 1:nrow(segs)){
		done.seg<-FALSE
		while(LTT[ii,2]<=segs[i,2]&&done.seg==FALSE){
			edge.length[i]<-edge.length[i]+
				LTT[ii,1]*(min(segs[i,2],LTT[ii,3])-
				max(segs[i,1],LTT[ii,2]))
			if(LTT[ii,3]>=segs[i,2]) done.seg<-TRUE
			if(LTT[ii,3]<=segs[i,2]) ii<-if(ii<nrow(LTT)) ii+1 else ii
		}
	}
	object<-list(segments=segs,nchanges=nchanges,edge.length=edge.length,tree=tree)
	class(object)<-"ctt"
	object
}	
	
plot.ctt<-function(x,...){
	h<-max(nodeHeights(x$tree))
	args<-list(...)
	if(!is.null(args$type)){ 
		type<-args$type
		args$type<-NULL
	} else type<-"rate"
	if(!is.null(args$show.tree)){
		show.tree<-args$show.tree
		args$show.tree<-NULL
	} else show.tree<-FALSE
	if(!is.null(args$add)){ 
		add<-args$add
		args$add<-NULL
	} else add<-FALSE
	if(is.null(args$ylim)) 
		args$ylim<-if(type=="number")c(0,max(x$nchanges)) else 
			c(0,max(x$nchanges/x$edge.length))
	if(is.null(args$xlim))
		args$xlim<-c(max(x$segments),min(x$segments))
	if(is.null(args$lwd)) args$lwd<-2
	if(is.null(args$xlab)) args$xlab<-"time since the present"
	if(is.null(args$ylab)) 
		args$ylab<-if(type=="number") "mean number of changes" else
			"mean number of changes / total edge length"
	args$type<-"l"
	args$x<-h-as.vector(t(x$segments))
	args$y<-if(type=="number") rbind(x$nchanges,x$nchanges) else 
		rbind(x$nchanges/x$edge.length,x$nchanges/x$edge.length)
	if(!add) do.call(plot,args)
	else do.call(lines,args)
	if(show.tree) plotTree(x$tree,add=TRUE,ftype="off",lwd=1,
		color=make.transparent("blue",0.1),mar=par()$mar,
		direction="leftwards",xlim=args$xlim)
}

sim.ctt<-function(tree,Q,anc=NULL,nmaps=100,...){
	x<-as.factor(sim.history(tree,Q,anc=anc,message=FALSE)$states)
	while(length(levels(x))!=ncol(Q)) 
		x<-as.factor(sim.history(tree,Q,anc=anc,message=FALSE)$states)
	flush.console()
	cat("Starting stochastic mapping with simulated data vector.... ")
	flush.console()
	trees<-make.simmap(tree,x,Q=Q,nsim=nmaps,message=FALSE)
	cat("Done.\n")
	flush.console()
	ctt(trees)
}

sim.multiCtt<-function(tree,Q,anc=NULL,nmaps=100,nsim=100,...){
	object<-replicate(nsim,sim.ctt(tree,Q,anc=anc,nmaps=nmaps,...),simplify=FALSE)
	class(object)<-"multiCtt"
	object
}

getChanges<-function(tree){
	states<-sort(unique(getStates(tree)))
	nc<-sapply(tree$maps,length)-1
	b<-which(nc>0)
	nc<-nc[b]
	xx<-vector()
	H<-nodeHeights(tree)
	for(i in 1:length(b)){
		for(j in 1:nc[i]){
			ss<-names(tree$maps[[b[i]]])[j+1]
			x<-rep(H[b[i],1]+cumsum(tree$maps[[b[i]]])[j],2)
			xx<-c(xx,setNames(x[1],
				paste(names(tree$maps[[b[i]]])[j:(j+1)],
				collapse="->")))
		}
	}
	xx
}

print.ctt<-function(x,...){
	cat("Object of class \"ctt\" consisting of:\n")
	cat("   (1) a matrix (segments) with the beginning & ending time of each segment.\n")
	cat("   (2) a vector (nchanges) with the mean number of changes in each segment.\n")
	cat("   (3) a vector (edge.length) containing the total edge length of each segement.\n")
	cat("   (4) an object of class \"phylo\".\n\n")
}

print.multiCtt<-function(x,...){
	cat(paste(length(x),"objects of class \"ctt\" in a list.\n\n"))
}

plot.multiCtt<-function(x,...){
	if(hasArg(alpha)) alpha<-list(...)$alpha
	else alpha<-0.05
	segments<-x[[1]]$segments
	nchanges<-sapply(x,function(x) x$nchanges)
	if(hasArg(type)) type<-list(...)$type
	else type<-"rate"
	edge.length<-sapply(x,function(x) x$edge.length)
	obj<-list(segments=segments,nchanges=rowMeans(nchanges),
		edge.length=rowMeans(edge.length),tree=x[[1]]$tree)
	class(obj)<-"ctt"
	lower<-max(floor(alpha/2*length(x)),1)
	upper<-min(ceiling((1-alpha/2)*length(x)),ncol(nchanges))
	xx<-max(nodeHeights(x[[1]]$tree))-as.vector(t(segments))
	xx<-c(xx,xx[length(xx):1])
	y.lower<-if(type=="number") apply(nchanges,1,sort)[lower,] else
		if(type=="rate") apply(nchanges/edge.length,1,sort)[lower,]
	y.upper<-if(type=="number") apply(nchanges,1,sort)[upper,] else
		if(type=="rate") apply(nchanges/edge.length,1,sort)[upper,]
	y.lower<-as.vector(rbind(y.lower,y.lower))
	y.upper<-as.vector(rbind(y.upper,y.upper))
	yy<-c(y.lower,y.upper[length(y.upper):1])
	args<-list(...)
	if(!is.null(args$alpha)) args$alpha<-NULL
	if(is.null(args$col)) args$col<-"blue"
	if(is.null(args$ylim)) args$ylim<-range(yy)
	args$x<-obj
	do.call(plot,args)
	polygon(xx,yy,col=make.transparent("grey",0.4),border=0)
}
