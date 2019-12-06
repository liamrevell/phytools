## function writes a modified "phylo" object to a simmap Newick string
## written by Liam Revell 2011, 2013, 2015, 2017, 2019

write.simmap<-function(tree,file=NULL,append=FALSE,map.order=NULL,quiet=FALSE,format="phylip",version=1.0){
	if(!(inherits(tree,"simmap")||inherits(tree,"multiSimmap"))) 
		stop("tree should be an object of class \"simmap\" or \"multiSimmap\".")
	if(inherits(tree,"multiPhylo")) N<-length(tree)
	else { 
		N<-1
		tree<-list(tree)
	}
	n<-sapply(tree,Ntip)
	if(format=="nexus"){
		trans<-unique(sort(sapply(tree,function(x) x$tip.label)))
		nn<-length(trans)
		write("#NEXUS",file)
		write(paste("[R-package PHYTOOLS, ",date(),"]\n",sep=""),file,append=TRUE)
		write("BEGIN TAXA;",file,append=TRUE)
		write(paste("\tDIMENSIONS NTAX = ",nn,";",sep=""),file,append=TRUE)
		write("\tTAXLABELS",file,append=TRUE)
		for(i in 1:length(trans)) write(paste("\t\t",trans[i],sep=""),file,append=TRUE)
		write("\t;",file,append=TRUE)
		write("END;",file,append=TRUE)
		if(version==1.0) write("BEGIN SMPTREES;\n\tTRANSLATE",file,append=TRUE)
		else write("BEGIN TREES;\n\tTRANSLATE",file,append=TRUE)
		for(i in 1:(nn-1)) write(paste("\t\t",i,"\t",trans[i],",",sep=""),file,
			append=TRUE)
		write(paste("\t\t",i+1,"\t",trans[i+1],sep=""),file,append=TRUE)
		write("\t;",file,append=TRUE)
		for(i in 1:N){
			tree[[i]]$tip.label<-sapply(tree[[i]]$tip.label,function(x) 
				which(x==trans))
			tree.string<-if(version==1.0) write.v1(tree[[i]],map.order=map.order,
				quiet=quiet) else write.v2(tree[[i]],map.order=map.order,
				quiet=quiet)
			write(paste("\tTREE * UNTITLED = [&R] ",tree.string,sep=""),file,
				append=TRUE)
		}
		write("END;",file,append=TRUE)
	} else {
		for(i in 1:N){ 
			if(version==1.0)
				write(write.v1(tree[[i]],map.order=map.order,quiet=quiet),file,append=append)
			else
				write(write.v2(tree[[i]],map.order=map.order,quiet=quiet),file,append=append)
		}
	}
}	
	
write.v1<-function(tree,map.order,quiet){
	if(is.null(map.order)){
		if(!is.null(attr(tree,"map.order")))
			map.order<-attr(tree,"map.order")
		else {
			if(!quiet) 
				message("map order should be specified in function call or by tree attribute \"map.order\".\nAssuming right-to-left order.")
			map.order<-"R"
		}
	}
	map.order<-toupper(unlist(strsplit(map.order,NULL))[1])
	if(map.order!="R"&&map.order!="L"){
		if(!quiet) message("do not recognize map order. Assuming right-to-left order.")
		map.order<-"R"
	}
	tree<-reorderSimmap(tree,"cladewise")
	n<-Ntip(tree)
	string<-vector()
	string[1]<-"("
	j<-2
	for(i in 1:nrow(tree$edge)){
		if(tree$edge[i,2]<=n){ 
			string[j]<-tree$tip.label[tree$edge[i,2]]
			j<-j+1
			string[j]<-":{"
			j<-j+1
			if(map.order=="L"){
				for(l in 1:length(tree$maps[[i]])){
					string[j]<-paste(c(names(tree$maps[[i]])[l],",",
						round(tree$maps[[i]][l],8)),collapse="")
					string[j+1]<-":"
					j<-j+2
				}
			} else {
				for(l in length(tree$maps[[i]]):1){
					string[j]<-paste(c(names(tree$maps[[i]])[l],",",
						round(tree$maps[[i]][l],8)),collapse="")
					string[j+1]<-":"
					j<-j+2
				}
			}					
			string[j-1]<-"}"
			v<-which(tree$edge[,1]==tree$edge[i,1])
			k<-i
			while(length(v)>0&&k==v[length(v)]){
				string[j]<-")"
				j<-j+1
				w<-which(tree$edge[,2]==tree$edge[k,1])
				if(length(w)>0){
					string[j]<-":{"
					j<-j+1
					if(map.order=="L"){
						for(l in 1:length(tree$maps[[w]])){
							string[j]<-paste(c(names(tree$maps[[w]])[l],",",
								round(tree$maps[[w]][l],8)),collapse="")
							string[j+1]<-":"
							j<-j+2
						}
					} else {
						for(l in length(tree$maps[[w]]):1){
							string[j]<-paste(c(names(tree$maps[[w]])[l],",",
								round(tree$maps[[w]][l],8)),collapse="")
							string[j+1]<-":"
							j<-j+2
						}
					}					
					string[j-1]<-"}"
				}
				v<-which(tree$edge[,1]==tree$edge[w,1])
				k<-w
			} 
			string[j]<-","
			j<-j+1
		} else if(tree$edge[i,2]>=n){
			string[j]<-"("
			j<-j+1
		}
	}
	string<-c(string[1:(length(string)-1)],";")
	string<-paste(string,collapse="")
	string
}

write.v2<-function(tree,map.order,quiet){
	if(is.null(map.order)){
		if(!is.null(attr(tree,"map.order")))
			map.order<-attr(tree,"map.order")
		else {
			if(!quiet) 
				message("map order should be specified in function call or by tree attribute \"map.order\".\nAssuming right-to-left order.")
			map.order<-"R"
		}
	}
	map.order<-toupper(unlist(strsplit(map.order,NULL))[1])
	if(map.order!="R"&&map.order!="L"){
		if(!quiet) message("do not recognize map order. Assuming right-to-left order.")
		map.order<-"R"
	}
	tree<-reorderSimmap(tree,"cladewise")
	n<-Ntip(tree)
	string<-vector()
	string[1]<-"("
	j<-2
	for(i in 1:nrow(tree$edge)){
		if(tree$edge[i,2]<=n){
			string[j]<-tree$tip.label[tree$edge[i,2]]
			j<-j+1
			string[j]<-":[&map={"
			j<-j+1
			nn<-length(tree$maps[[i]])
			if(nn==1){ 
				string[j]<-names(tree$maps[[i]])[1]
				j<-j+1
			} else {
				if(map.order=="L"){
					for(l in 1:(nn-1)){
						string[j]<-paste(c(names(tree$maps[[i]])[l],",",
							round(tree$maps[[i]][l],8)),collapse="")
						string[j+1]<-","
						j<-j+2
					}
					string[j]<-names(tree$maps[[i]])[nn]
					j<-j+1
				} else {
					for(l in nn:2){
						string[j]<-paste(c(names(tree$maps[[i]])[l],",",
							round(tree$maps[[i]][l],8)),collapse="")
						string[j+1]<-","
						j<-j+2
					}
					string[j]<-names(tree$maps[[i]])[1]
					j<-j+1
				}					
			}
			string[j]<-paste(c("}]",round(tree$edge.length[i],8)),collapse="")
			j<-j+1
			v<-which(tree$edge[,1]==tree$edge[i,1])
			k<-i
			while(length(v)>0&&k==v[length(v)]){
				string[j]<-")"
				j<-j+1
				w<-which(tree$edge[,2]==tree$edge[k,1])
				if(length(w)>0){
					nn<-length(tree$maps[[w]])
					string[j]<-":[&map={"
					j<-j+1
					if(nn==1){ 
						string[j]<-names(tree$maps[[w]])[1]
						j<-j+1
					} else {	
						if(map.order=="L"){
							for(l in 1:(nn-1)){
								string[j]<-paste(c(names(tree$maps[[w]])[l],",",
									round(tree$maps[[w]][l],8)),collapse="")
								string[j+1]<-","
								j<-j+2
							}
							string[j]<-names(tree$maps[[w]])[nn]
							j<-j+1
						} else {
							for(l in nn:2){
								string[j]<-paste(c(names(tree$maps[[w]])[l],",",
									round(tree$maps[[w]][l],8)),collapse="")
								string[j+1]<-","
								j<-j+2
							}
							string[j]<-names(tree$maps[[w]])[1]
							j<-j+1
						}
					}
					string[j]<-paste(c("}]",round(tree$edge.length[w],8)),collapse="")
					j<-j+1
				}
				v<-which(tree$edge[,1]==tree$edge[w,1])
				k<-w
			}
			string[j]<-","
			j<-j+1

		} else if(tree$edge[i,2]>=n){
			string[j]<-"("
			j<-j+1
		}
	}
	string<-c(string[1:(length(string)-1)],";")
	string<-paste(string,collapse="")
	string
}
