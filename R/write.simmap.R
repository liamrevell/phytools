# function writes a modified "phylo" object to a simmap Newick string
# written by Liam Revell 2011, 2013, 2015

write.simmap<-function(tree,file=NULL,append=FALSE,map.order=NULL,quiet=FALSE){
	if(inherits(tree,"multiPhylo")) for(i in 1:length(tree)) write.simmap(tree,file,if(i==1) append else TRUE,map.order)
	else {
		if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\" or \"multiPhylo\".")
		if(is.null(tree$maps)) stop("tree is does not contain a stochastic character map.")
		if(is.null(map.order)){
			if(!is.null(attr(tree,"map.order")))
				map.order<-attr(tree,"map.order")
			else {
				if(!quiet) message("map order should be specified in function call or by tree attribute \"map.order\".\nAssuming right-to-left order.")
				map.order<-"R"
			}
		}
		map.order<-toupper(unlist(strsplit(map.order,NULL))[1])
		if(map.order!="R"&&map.order!="L"){
			if(!quiet) message("do not recognize map order. Assuming right-to-left order.")
			map.order<-"R"
		}
		tree<-reorderSimmap(tree,"cladewise")
		n<-length(tree$tip)
		string<-vector(); string[1]<-"("; j<-2
		for(i in 1:nrow(tree$edge)){
			if(tree$edge[i,2]<=n){ 
				string[j]<-tree$tip.label[tree$edge[i,2]]; j<-j+1
				string[j]<-":{"; j<-j+1
				if(map.order=="L"){
					for(l in 1:length(tree$maps[[i]])){
						string[j]<-paste(c(names(tree$maps[[i]])[l],",",round(tree$maps[[i]][l],8)),collapse="")
						string[j+1]<-":"; j<-j+2
					}
				} else {
					for(l in length(tree$maps[[i]]):1){
						string[j]<-paste(c(names(tree$maps[[i]])[l],",",round(tree$maps[[i]][l],8)),collapse="")
						string[j+1]<-":"; j<-j+2
					}
				}					
				string[j-1]<-"}"
				v<-which(tree$edge[,1]==tree$edge[i,1]); k<-i
				while(length(v)>0&&k==v[length(v)]){
					string[j]<-")"; j<-j+1
					w<-which(tree$edge[,2]==tree$edge[k,1])
					if(length(w)>0){
						string[j]<-":{"; j<-j+1
						if(map.order=="L"){
							for(l in 1:length(tree$maps[[w]])){
								string[j]<-paste(c(names(tree$maps[[w]])[l],",",round(tree$maps[[w]][l],8)),collapse="")
								string[j+1]<-":"; j<-j+2
							}
						} else {
							for(l in length(tree$maps[[w]]):1){
								string[j]<-paste(c(names(tree$maps[[w]])[l],",",round(tree$maps[[w]][l],8)),collapse="")
								string[j+1]<-":"; j<-j+2
							}
						}					
						string[j-1]<-"}"
					}
					v<-which(tree$edge[,1]==tree$edge[w,1]); k<-w
				} 
				string[j]<-","; j<-j+1
			} else if(tree$edge[i,2]>=n){
				string[j]<-"("; j<-j+1
			}
		}
		string<-c(string[1:(length(string)-1)],";")
		string<-paste(string,collapse="")
		if(is.null(file)) return(string)
		else write(string,file=file,append=append)
	}
}

