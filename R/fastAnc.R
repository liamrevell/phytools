## function does fast estimation of ML ancestral states using ace
## written by Liam J. Revell 2012, 2013, 2015, 2019, 2020, 2021, 2023, 2026

fastAnc<-function(tree,x,vars=FALSE,CI=FALSE,...){
	if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\".")
	if(is.null(names(x))){
		warn<-paste("x should be a vector with names corresponding to the taxon labels of the tree.\n",
			" Assuming x is in the order of tree$tip.label (this is seldom true).")
		warning(warn)
	}
	if(length(class(tree)>1)) class(tree)<-"phylo"
	if(hasArg(anc.states)) anc.states<-list(...)$anc.states
	else anc.states<-NULL
	if(!is.null(anc.states)){
		nodes<-as.numeric(names(anc.states))
		tt<-tree
		for(i in 1:length(nodes)){
			M<-matchNodes(tt,tree,method="distances",quiet=TRUE)
			ii<-M[which(M[,2]==nodes[i]),1]
			tt<-bind.tip(tt,nodes[i],edge.length=0,where=ii)
		}
		x<-c(x,anc.states)
	} else tt<-tree
	if(!is.binary(tt)) btree<-multi2di(tt,random=FALSE)
	else btree<-tt
	M<-btree$Nnode
	N<-length(btree$tip.label)
	anc<-v<-vector()
	for(i in 1:M+N){
		a<-collapse.singles(multi2di(ape::root.phylo(btree,node=i),random=FALSE))
   		anc[i-N]<-ace(x,a,method="pic")$ace[1]
   		names(anc)[i-N]<-i
		if(vars||CI){
			picx<-pic(x,a,rescaled.tree=TRUE)
			b<-picx$rescaled.tree
			d<-which(b$edge[,1]==(length(b$tip.label)+1))
			v[i-N]<-(1/b$edge.length[d[1]]+1/b$edge.length[d[2]])^(-1)*mean(picx$contr^2)
			names(v)[i-N]<-names(anc)[i-N]
		}
 	}
	if(!is.binary(tree)||!is.null(anc.states)){
		ancNames<-matchNodes(tree,btree,method="distances",quiet=TRUE)
		anc<-anc[as.character(ancNames[,2])]
		names(anc)<-ancNames[,1]
		if(vars||CI){ 
			v<-v[as.character(ancNames[,2])]
			names(v)<-ancNames[,1]
		}
	}
	obj<-list(ace=anc)
	if(vars) obj$var<-v
	if(CI){ 
		obj$CI95<-cbind(anc-1.96*sqrt(v),anc+1.96*sqrt(v))
		rownames(obj$CI95)<-names(anc)
	}
	if(length(obj)==1) obj<-obj$ace
	class(obj)<-"fastAnc"
	attr(obj,"tree")<-tree
	attr(obj,"data")<-x
	obj
}

## print method for "fastAnc"
## written by Liam J. Revell 2015, 2026
print.fastAnc<-function(x,digits=6,printlen=NULL,...){
	cat("Ancestral character estimates using fastAnc:\n")
	if(!is.list(x)){ 
		if(is.null(printlen)||printlen>=length(x)) 
			print(round(unclass(x[1:length(x)]),digits)) 
		else printDotDot(unclass(x[1:length(x)]),digits,
			printlen)
	} else {
		Nnode<-length(x$ace)
		if(is.null(printlen)||printlen>=Nnode) print(round(x$ace,digits))
		else printDotDot(x$ace,digits,printlen)
		if(!is.null(x$var)){
			cat("\nVariances on ancestral states:\n")
			if(is.null(printlen)||printlen>=Nnode) print(round(x$var,digits))
			else printDotDot(x$var,digits,printlen)
		}
		if(!is.null(x$CI95)){
			cat("\nLower & upper 95% CIs:\n")
			colnames(x$CI95)<-c("lower","upper")
			if(is.null(printlen)||printlen>=Nnode) print(round(x$CI95,digits))
			else printDotDot(x$CI95,digits,printlen)
		}
	}
	cat("\n")
}

## internal function
## written by Liam J. Revell 2015
printDotDot<-function(x,digits,printlen){
	if(is.vector(x)){
		x<-as.data.frame(t(as.matrix(unclass(round(x[1:printlen],digits)))))
		x<-cbind(x,"....")
		rownames(x)<-""
		colnames(x)[printlen+1]<-""
		print(x)
	} else if(is.matrix(x)){
		x<-as.data.frame(rbind(round(x[1:printlen,],digits),c("....","....")))
		print(x)
	}
}

## S3 plot method for fastAnc object class
plot.fastAnc<-function(x,...){
  plotTree(attr(x,"tree"),...)
  if(hasArg(tip.cex)) tip.cex<-list(...)$tip.cex
  else tip.cex<-1
  if(hasArg(node.cex)) node.cex<-list(...)$node.cex
  else node.cex<-1.5
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  gradient<-hcl.colors(n=101)
  ace <- if(!is.list(x)) x[1:length(x)] else x$ace[1:length(x$ace)]
  x_range<-range(c(attr(x,"data"),ace))
  tip.cols<-gradient[round((attr(x,"data")[attr(x,"tree")$tip.label]-
      x_range[1])/diff(x_range)*100)]
  node.cols<-gradient[round((ace-x_range[1])/diff(x_range)*100)]
  points(pp$xx,pp$yy,col=c(tip.cols,node.cols),pch=16,
    cex=c(rep(tip.cex,Ntip(attr(x,"tree"))),
      rep(node.cex,Nnode(attr(x,"tree")))))
  if(hasArg(legend)) legend<-list(...)$legend
  else legend<-"topleft"
  if(legend!=FALSE){
    h<-diff(par()$usr[1:2])*0.4
    if(legend=="topleft"){
      x.pos<-par()$usr[1]+0.1*h
      y.pos<-par()$usr[4]-0.05*diff(par()$usr[3:4])
    } else if(legend=="bottomleft"){
      x.pos<-par()$usr[1]+0.1*h
      y.pos<-par()$usr[3]+0.05*diff(par()$usr[3:4])
    } else if(legend=="topright"){
      x.pos<-par()$usr[2]-h
      y.pos<-par()$usr[4]-0.05*diff(par()$usr[3:4])
    } else if(legend=="bottomright"){
      x.pos<-par()$usr[2]-h
      y.pos<-par()$usr[3]+0.05*diff(par()$usr[3:4])
    }
    if(hasArg(title)) title<-list(...)$title
    else title<-"trait value"
    add.color.bar(0.75*h,gradient,
      title=title,lims=x_range,digits=3,prompt=FALSE,
      lwd=6,outline=FALSE,x=x.pos,y=y.pos,subtitle="",
	  fsize=0.8)
  }
}
