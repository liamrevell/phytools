# function writes a "phylo" object to a Newick string with ancestor state estimates
# written by Liam J. Revell 2013, 2015

writeAncestors<-function(tree,Anc=NULL,file="",digits=6,format=c("phylip","nexus"),...){
	format=format[1]
	if(hasArg(CI)) CI<-list(...)$CI
	else CI<-TRUE
	if(is.null(Anc)){
		if(hasArg(x)){ 
			x<-list(...)$x
			if(inherits(tree,"multiPhylo")){
				if(is.list(x)) Anc<-mapply(fastAnc,tree,x,MoreArgs=list(CI=CI),SIMPLIFY=FALSE)
				else Anc<-lapply(tree,fastAnc,x=x,CI=CI)
			} else if(inherits(tree,"phylo")){
				if(is.list(x)){ 
					Anc<-lapply(x,fastAnc,tree=tree,CI=CI)
					tree<-repPhylo(tree,length(x))
					class(tree)<-"multiPhylo"
				} else Anc<-fastAnc(tree,x,CI=CI)
			}
		} else stop("must have argument 'Anc' or 'x'")
	}
	if(format=="phylip"){
		if(class(tree)=="multiPhylo")
			XX<-mapply(writeAnc,tree,Anc,MoreArgs=list(digits=digits))
		else if(class(tree)=="phylo")
			XX<-writeAnc(tree,Anc,digits)
		else stop("tree should be an object of class 'phylo' or 'multiPhylo'")
		write(XX,file)
		invisible(XX)
	} else if(format=="nexus"){
		writeNex(tree,Anc,file,digits)
	}
}

# internal function to create a Nexus style output file
# written by Liam J. Revell 2013
writeNex<-function(tree,Anc,file="",digits){
	if(inherits(tree,"multiPhylo")) N<-length(tree)
	else { 
		N<-1
		tree<-list(tree)
		Anc<-list(Anc)
	}
	n<-length(tree[[1]]$tip.label)
	write("#NEXUS",file)
	write(paste("[R-package PHYTOOLS, ",date(),"]\n",sep=""),file,append=TRUE)
	write("BEGIN TAXA;",file,append=TRUE)
	write(paste("\tDIMENSIONS NTAX = ",n,";",sep=""),file,append=TRUE)
	write("\tTAXLABELS",file,append=TRUE)
	trans<-tree[[1]]$tip.label; trans<-sort(trans)
	for(i in 1:n) write(paste("\t\t",trans[i],sep=""),file,append=TRUE)
	write("\t;",file,append=TRUE)
	write("END;",file,append=TRUE)
	write("BEGIN TREES;\n\tTRANSLATE",file,append=TRUE)
	for(i in 1:(n-1)) write(paste("\t\t",i,"\t",trans[i],",",sep=""),file,append=TRUE)
	write(paste("\t\t",i+1,"\t",trans[i+1],sep=""),file,append=TRUE)
	write("\t;",file,append=TRUE)
	for(i in 1:N){
		tree[[i]]$tip.label<-sapply(tree[[i]]$tip.label,function(x) which(x==trans))
		write(paste("\tTREE * UNTITLED = [&R] ",writeAnc(tree[[i]],Anc[[i]],digits),sep=""),file,append=TRUE)
	}
	write("END;",file,append=TRUE)
}

# internal function to write the Newick string with ancestor states
# written by Liam J. Revell
writeAnc<-function(tree,Anc,digits){
	tree<-reorder.phylo(tree,"cladewise")
	n<-length(tree$tip.label)
	if(!is.list(Anc)) Anc<-list(ace=Anc)
	Anc$ace<-Anc$ace[order(names(Anc$ace))]
	Anc$ace<-round(Anc$ace,digits)
	if(!is.null(Anc$CI95)){
		Anc$CI95<-round(Anc$CI95,digits)
		tree$node.label<-paste("[&CI={",Anc$CI95[,1],",",Anc$CI95[,2],"},ancstate={",Anc$ace,"}]",sep="")
	} else {
		Anc$CI95<-Anc$CI95[names(Anc$ace),]
		tree$node.label<-paste("[&ancstate={",Anc$ace,"}]",sep="")
	}
	tree$edge.length<-round(tree$edge.length,digits)
	string<-vector(); string[1]<-"("; j<-2
	for(i in 1:nrow(tree$edge)){
		if(tree$edge[i,2]<=n){ 
			string[j]<-tree$tip.label[tree$edge[i,2]]; j<-j+1
			if(!is.null(tree$edge.length)){
				string[j]<-paste(c(":",tree$edge.length[i]),collapse="") 
				j<-j+1
			}
			v<-which(tree$edge[,1]==tree$edge[i,1]); k<-i
			while(length(v)>0&&k==v[length(v)]){
				string[j]<-")"; j<-j+1
				string[j]<-tree$node.label[tree$edge[k,1]-n]; j<-j+1
				w<-which(tree$edge[,2]==tree$edge[k,1])
				if(!is.null(tree$edge.length)){
					string[j]<-paste(c(":",tree$edge.length[w]),collapse="")
					j<-j+1
				}
				v<-which(tree$edge[,1]==tree$edge[w,1]); k<-w 
			} 
			string[j]<-","; j<-j+1
		} else if(tree$edge[i,2]>=n){
			string[j]<-"("; j<-j+1
		}
	}
	if(is.null(tree$edge.length)) string<-c(string[1:(length(string)-1)],";")
	else string<-c(string[1:(length(string)-2)],";")
	string<-paste(string,collapse="")
	return(string)
}

