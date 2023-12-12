## functions to compute the parsimony score (mostly for diagnostic purposes)

Fitch<-function(x,pw,nn,nm,return="score"){
	xx<-vector("list",Ntip(pw)+Nnode(pw))
	xx[1:Ntip(pw)]<-setNames(x,nm)[pw$tip.label]
	pp<-0
	for(i in 1:length(nn)){
		ee<-which(pw$edge[,1]==nn[i])
		Intersection<-Reduce(intersect,xx[pw$edge[ee,2]])
		if(length(Intersection)>0){ 
			xx[[nn[i]]]<-Intersection
		} else {
			xx[[nn[i]]]<-Reduce(union,xx[pw$edge[ee,2]])
			pp<-pp+1
		}
	}
	if(return=="score") pp else if(return=="nodes") xx
}

pscore<-function(tree,x,...){
	pw<-if(!is.null(attr(tree,"order"))&&
			attr(tree,"order")=="postorder") tree else 
		reorder(tree,"postorder")
	nn<-unique(pw$edge[,1])
	if(is.matrix(x)||is.data.frame(x)){
		nm<-rownames(x)
		apply(x,2,Fitch,pw=pw,nn=nn,nm=nm,...)
	} else {
		nm<-names(x)
		Fitch(x,pw,nn,nm,...)
	}
}

acctran<-function(tree,x,...){
	pw<-if(!is.null(attr(tree,"order"))&&
			attr(tree,"order")=="postorder") tree else 
		reorder(tree,"postorder")
	cw<-if(!is.null(attr(tree,"order"))&&
			attr(tree,"order")=="cladewise") tree else
		reorder(tree,"cladewise")
	nn<-unique(pw$edge[,1])
	if(is.matrix(x)||is.data.frame(x)){
		nm<-rownames(x)
		Nodes<-apply(x,2,AccDelTran,pw=pw,cw=cw,nn=nn,nm=nm,
			method="acctran")
	} else {
		nm<-names(x)
		AccDelTran(x,pw,cw,nn,nm,method="acctran")
	}
}

AccDelTran<-function(x,pw,cw,nn,nm,method="acctran"){
	Sets<-Fitch(x,pw,nn,nm,return="nodes")
  ## preorder traversal
	for(i in length(nn):1){
	  ee<-which(cw$edge[,1]==nn[i])
	  
	}

#	ee<-which(pw$edge[,1]==nn[i])
#	Intersection<-Reduce(intersect,xx[pw$edge[ee,2]])
#	if(length(Intersection)>0){ 
#	  xx[[nn[i]]]<-Intersection
#	} else {
#	  xx[[nn[i]]]<-Reduce(union,xx[pw$edge[ee,2]])
#	  pp<-pp+1
#	}
	
	
}