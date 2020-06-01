## This function computes the data for a lineage-through-time plot and 
## (optionally) creates this plot the function does not require a tree 
## that is ultrametric.  Optionally, the function can remove extinct
## species from the phylogeny. If the input tree is an object of class 
## "multiPhylo" then the function will simultaneously plot all ltts.
## written by Liam J. Revell 2010-2015

ltt<-function(tree,plot=TRUE,drop.extinct=FALSE,log.lineages=TRUE,gamma=TRUE,...){
	# set tolerance
	tol<-1e-6
	# check "phylo" object
	if(!inherits(tree,"phylo")&&!inherits(tree,"multiPhylo"))
		stop("tree must be object of class \"phylo\" or \"multiPhylo\".")
	if(inherits(tree,"multiPhylo")){
		obj<-lapply(tree,ltt,plot=FALSE,drop.extinct=drop.extinct,
			log.lineages=log.lineages,gamma=gamma)
		class(obj)<-"multiLtt"
	} else {
		# reorder the tree
		tree<-reorder.phylo(tree,order="cladewise")
		if(!is.null(tree$node.label)){ 
			node.names<-setNames(tree$node.label,1:tree$Nnode+Ntip(tree))
			tree$node.label<-NULL
		} else node.names<-NULL
		## check if tree is ultrametric & if yes, then make *precisely* ultrametric
		if(is.ultrametric(tree)){
			h<-max(nodeHeights(tree))
			time<-c(0,h-sort(branching.times(tree),decreasing=TRUE),h)
			nodes<-as.numeric(names(time)[2:(length(time)-1)])
			ltt<-c(cumsum(c(1,sapply(nodes,function(x,y) sum(y==x)-1,y=tree$edge[,1]))),
				length(tree$tip.label))
			names(ltt)<-names(time)
		} else {
			# internal functions
			# drop extinct tips
			drop.extinct.tips<-function(phy){
				temp<-diag(vcv(phy))
				if(length(temp[temp<(max(temp)-tol)])>0)
					pruned.phy<-drop.tip(phy,names(temp[temp<(max(temp)-tol)]))
				else
					pruned.phy<-phy
				return(pruned.phy)
			}
			# first, if drop.extinct==TRUE
			if(drop.extinct==TRUE)
				tree<-drop.extinct.tips(tree)
			# compute node heights
			root<-length(tree$tip)+1
			node.height<-matrix(NA,nrow(tree$edge),2)
			for(i in 1:nrow(tree$edge)){
				if(tree$edge[i,1]==root){
					node.height[i,1]<-0.0
					node.height[i,2]<-tree$edge.length[i]
				} else {
					node.height[i,1]<-node.height[match(tree$edge[i,1],tree$edge[,2]),2]
					node.height[i,2]<-node.height[i,1]+tree$edge.length[i]
				}
			}
			ltt<-vector()
			tree.length<-max(node.height) # tree length
			n.extinct<-sum(node.height[tree$edge[,2]<=length(tree$tip),2]<(tree.length-tol))
			# fudge things a little bit
			node.height[tree$edge[,2]<=length(tree$tip),2]<-
				node.height[tree$edge[,2]<=length(tree$tip),2]+1.1*tol
			time<-c(0,node.height[,2]); names(time)<-as.character(c(root,tree$edge[,2]))
			temp<-vector()
			time<-time[order(time)]
			time<-time[1:(tree$Nnode+n.extinct+1)] # times
			# get numbers of lineages
			for(i in 1:(length(time)-1)){
				ltt[i]<-0
				for(j in 1:nrow(node.height))
					ltt[i]<-ltt[i]+(time[i]>=(node.height[j,1]-
						tol)&&time[i]<=(node.height[j,2]-tol))
			}
			ltt[i+1]<-0
			for(j in 1:nrow(node.height))
				ltt[i+1]<-ltt[i+1]+(time[i+1]<=(node.height[j,2]+tol))
			# set names (these are the node indices)
			names(ltt)<-names(time)
			# append 0,1 for time 0
			ltt<-c(1,ltt)
			time<-c(0,time)
			# subtract fudge factor
			time[length(time)]<-time[length(time)]-1.1*tol
		}
		if(!is.null(node.names)){ nn<-sapply(names(time),function(x,y) 
			if(any(names(y)==x)) y[which(names(y)==x)] else "",y=node.names)
			names(ltt)<-names(time)<-nn
		}
		if(gamma==FALSE){
			obj<-list(ltt=ltt,times=time,tree=tree)
			class(obj)<-"ltt"
		} else {
			gam<-gammatest(list(ltt=ltt,times=time))
			obj<-list(ltt=ltt,times=time,gamma=gam$gamma,p=gam$p,tree=tree)
			class(obj)<-"ltt"
		}
	}
	if(plot) plot(obj,log.lineages=log.lineages,...)
	obj
}

## function computes the gamma-statistic & a two-tailed P-value
## written by Liam J. Revell 2011, 2019

gammatest<-function(x){
	n<-max(x$ltt)
	g<-vector()
	for(i in 2:(length(x$times))) g[i-1]<-x$times[i]-x$times[i-1]
	T<-sum((2:n)*g[2:n])
	doublesum<-0
	for(i in 1:(n-1)) for(k in 1:i) doublesum<-doublesum+k*g[k]
	gamma<-(1/(n-2)*doublesum-T/2)/(T*sqrt(1/(12*(n-2))))
	p<-2*pnorm(abs(gamma),lower.tail=F)
	object<-list(gamma=gamma,p=p)
	class(object)<-"gammatest"
	object
}

## print method for object class "gammatest"
## written by Liam J. Revell 2019

print.gammatest<-function(x,...){
	if(hasArg(digits)) digits<-list(...)$digits
	else digits<-4
	cat("\nAn object of class \"gammatest\" with:\n")
	cat(paste("(1) Pybus & Harvey's gamma = ",
		round(x$gamma,digits),sep=""))
	cat(paste("\n(2) p-value = ",round(x$p,digits),
		"\n\n",sep=""))
}

## S3 print method for object of class "ltt"
## written by Liam J. Revell 2015

print.ltt<-function(x,digits=4,...){
	cat("Object of class \"ltt\" containing:\n\n")
	cat(paste("(1) A phylogenetic tree with ",Ntip(x$tree), 
		" tips and ",x$tree$Nnode," internal\n",sep= ""))
	cat("    nodes.\n\n")
	cat(paste("(2) Vectors containing the number of lineages (ltt) and\n",
		"    branching times (times) on the tree.\n\n",sep=""))
	if(!is.null(x$gamma))
		cat(paste("(3) A value for Pybus & Harvey's \"gamma\"",
			" statistic of\n    gamma = ",round(x$gamma,digits),
			", p-value = ",
			round(x$p,digits),".\n\n",sep=""))
}

## S3 print method for object of class "multiLtt"
## written by Liam J. Revell 2015

print.multiLtt<-function(x,...){
	cat(paste(length(x),"objects of class \"ltt\" in a list\n"))
}

## S3 plot method for object of class "ltt"
## written by Liam J. Revell 2015

plot.ltt<-function(x,...){
	args<-list(...)
	args$x<-x$time
	if(!is.null(args$log.lineages)){ 
		logl<-args$log.lineages
		args$log.lineages<-NULL
	} else logl<-TRUE
	if(!is.null(args$show.tree)){
		show.tree<-args$show.tree
		args$show.tree<-NULL
	} else show.tree<-FALSE
	if(!is.null(args$transparency)){
		transparency<-args$transparency
		args$transparency<-NULL
	} else transparency<-0.3
	args$y<-if(logl) log(x$ltt) else x$ltt
	if(is.null(args$xlab)) args$xlab<-"time"
	if(is.null(args$ylab)) 
		args$ylab<-if(logl) "log(lineages)" else "lineages"
	if(is.null(args$type)) args$type<-"s"
	if(hasArg(add)){ 
		add<-list(...)$add
		args$add<-NULL
	} else add<-FALSE
	if(!add) do.call(plot,args)
	else do.call(lines,args)
	if(show.tree){
		tips<-if(par()$ylog) setNames(exp(1:Ntip(x$tree)),x$tree$tip.label) 
			else setNames(1:Ntip(x$tree),x$tree$tip.label)
		plotTree(x$tree,color=rgb(0,0,1,transparency),
			ftype="off",add=TRUE,mar=par()$mar,tips=tips)
	}
}

## S3 plot method for object of class "multiLtt"
## written by Liam J. Revell 2015, 2020

plot.multiLtt<-function(x,...){
	max.lineages<-max(sapply(x,function(x) max(x$ltt)))
	max.time<-max(sapply(x,function(x) max(x$time)))
	args<-list(...)
	if(!is.null(args$log.lineages)) logl<-list(...)$log.lineages
	else logl<-TRUE	
	if(is.null(args$xlim)) args$xlim<-c(0,max.time)
	if(is.null(args$ylim)) 
		args$ylim<-if(logl) c(0,log(max.lineages)) else c(1,max.lineages)
	args$x<-x[[1]]
	do.call(plot,args)
	args$add<-TRUE
	if(!is.null(args$log)) args$log<-NULL
	if(length(x)>2){
		for(i in 2:length(x)){
			args$x<-x[[i]]
			do.call(plot,args)
		}
	} else {
		args$x<-x[[2]]
		do.call(plot,args)
	}
}

## compute the gamma statistic through time by slicing the tree n times
## written by Liam J. Revell 2017

gtt<-function(tree,n=100,...){
	if(!inherits(tree,"phylo")) 
		stop("tree must be object of class \"phylo\".")
	if(inherits(tree,"simmap")) tree<-as.phylo(tree)
	if(hasArg(plot)) plot<-list(...)$plot
	else plot<-FALSE
	obj<-ltt(tree,plot=FALSE)
	t<-obj$times[which(obj$ltt==3)[1]]
	h<-max(nodeHeights(tree))
	x<-seq(t,h,by=(h-t)/(n-1))
	trees<-lapply(x,treeSlice,tree=tree,orientation="rootwards")
	gamma<-sapply(trees,function(x,plot){ 
		obj<-unlist(gammatest(ltt<-ltt(x,plot=FALSE)));
		if(plot) plot(ltt,xlim=c(0,h),ylim=c(1,Ntip(tree)),
			log.lineages=FALSE,log="y");
		Sys.sleep(0.01);
		obj},plot=plot)
	object<-list(t=x,gamma=gamma[1,],p=gamma[2,],tree=tree)
	class(object)<-"gtt"
	object
}

## plot method for "gtt" object class
plot.gtt<-function(x,...){
	args<-list(...)
	args$x<-x$t
	args$y<-x$gamma
	if(!is.null(args$show.tree)){ 
		show.tree<-args$show.tree
		args$show.tree<-NULL
	} else show.tree<-TRUE
	if(is.null(args$xlim)) args$xlim<-c(0,max(x$t))
	if(is.null(args$xlab)) args$xlab<-"time"
	if(is.null(args$ylab)) args$ylab<-expression(gamma)
	if(is.null(args$lwd)) args$lwd<-3
	if(is.null(args$type)) args$type<-"s"
	if(is.null(args$bty)) args$bty<-"l"
	if(is.null(args$main)) args$main<-expression(paste(gamma," through time plot"))
	do.call(plot,args)
	if(show.tree) plotTree(x$tree,add=TRUE,ftype="off",mar=par()$mar,
		xlim=args$xlim,color=make.transparent("blue",0.1))
}

## print method for "gtt" object class
print.gtt<-function(x,...)
	cat("Object of class \"gtt\".\n\n")
	
## perform the MCCR test of Pybus & Harvey (2000)
## written by Liam J. Revell 2018

mccr<-function(obj,rho=1,nsim=100,...){
	N<-round(Ntip(obj$tree)/rho)
	tt<-pbtree(n=N,nsim=nsim)
	foo<-function(x,m) drop.tip(x,sample(x$tip.label,m))
	tt<-lapply(tt,foo,m=N-Ntip(obj$tree))
	g<-sapply(tt,function(x) ltt(x,plot=FALSE)$gamma)
	P<-if(obj$gamma>median(g)) 2*mean(g>=obj$gamma) else 2*mean(g<=obj$gamma)
	result<-list(gamma=obj$gamma,"P(two-tailed)"=P,null.gamma=g)
	class(result)<-"mccr"
	result
}

## print method for "mccr" object class

print.mccr<-function(x,digits=4,...){
	cat("Object of class \"mccr\" consisting of:\n\n")
	cat(paste("(1) A value for Pybus & Harvey's \"gamma\"",
		" statistic of \n    gamma = ",round(x$gamma,digits),
		".\n\n",sep=""))
	cat(paste("(2) A two-tailed p-value from the MCCR test of ",
		round(x$'P(two-tailed)',digits),".\n\n", sep = ""))
	cat(paste("(3) A simulated null-distribution of gamma from ",
		length(x$null.gamma),"\n    simulations.\n\n",sep=""))
}

## plot method for "mccr" object class

plot.mccr<-function(x,...){
	if(hasArg(main)) main<-list(...)$main
	else main=expression(paste("null distribution of ",
		gamma))
	hh<-hist(x$null.gamma,breaks=min(c(max(12,
		round(length(x$null.gamma)/10)),20)),
		plot=FALSE)
	if(hasArg(ylim)) ylim<-list(...)$ylim
	else ylim<-c(0,1.15*max(hh$counts))
	plot(hh,xlim=range(c(x$gamma,x$null.gamma)),
		main=main,xlab=expression(gamma),col="lightgrey",
		ylim=ylim)
	arrows(x0=x$gamma,y0=par()$usr[4],y1=0,length=0.12,
		col=make.transparent("blue",0.5),lwd=2)
	text(x$gamma,0.96*par()$usr[4],
		expression(paste("observed value of ",gamma)),
		pos=if(x$gamma>mean(x$null.gamma)) 2 else 4,
		cex=0.9)
}