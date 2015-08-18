# function
# written by Liam J. Revell 2012, 2015

writeNexus<-function(tree,file=""){
	if(inherits(tree,"multiPhylo")) N<-length(tree)
	else { 
		N<-1
		tree<-list(tree)
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
		write(paste("\tTREE * UNTITLED = [&R] ",write.tree(tree[[i]]),sep=""),file,append=TRUE)
	}
	write("END;",file,append=TRUE)
}
