## makes xml data and tree file for SIMMAP
## written by Liam J. Revell 2012, 2015, 2022
 
export.as.xml<-function(file,trees,X){
	if(is.vector(X)) X<-data.frame(X)
	if(inherits(X,"DNAbin")) { 
		X<-as.character(X)
		datatype="nucleotide"
	} else datatype="standard"
	if(inherits(trees,"phylo")){ 
		trees<-list(trees)
		class(trees)<-"multiPhylo"
	}
	ntaxa<-length(trees[[1]]$tip)
	nchars<-ncol(X)
	write("<?xml version=\"1.0\"?>",file)
	write("<simmap>",file,append=TRUE)
	write(paste("\t<data ntaxa=\"",ntaxa,"\" nchars=\"",nchars,"\" datatype=\"",datatype,"\">",sep=""),file,append=TRUE)
	for(i in 1:ntaxa) write(paste("\t\t<seq name=\"",rownames(X)[i],"\">",paste(X[i,],collapse=""),"</seq>",sep=""),file,append=TRUE)
	write("\t</data>",file,append=TRUE)
	write("\t<trees>",file,append=TRUE)
	trans<-trees[[1]]$tip.label
	for(i in 1:ntaxa) write(paste("\t\t<translate id=\"",i,"\">",trans[i],"</translate>",sep=""),file,append=TRUE)
	for(i in 1:length(trees)){
		trees[[i]]$tip.label<-match(trans,trees[[i]]$tip.label)
		temp<-unlist(strsplit(write.tree(trees[[i]]),NULL))
		temp<-paste(temp[1:(length(temp)-1)],collapse="")
		write(paste("\t\t<tree>",temp,"</tree>",sep=""),file,append=TRUE)
	}
	write("\t</trees>",file,append=TRUE)
	write("</simmap>",file,append=TRUE)
}
