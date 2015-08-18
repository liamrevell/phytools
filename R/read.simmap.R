# function reads SIMMAP v1.0 & v1.5 style trees into file
# written by Liam J. Revell 2011-2015

read.simmap<-function(file="",text,format="nexus",rev.order=TRUE,version=1.0){
	if(file!=""){
		trans<-NULL
		if(format=="nexus"){ 
			XX<-readNexusData(file,version)
			text<-XX$text
			trans<-XX$trans
			Ntree<-XX$Ntree
		} else if(format=="phylip"){
			text<-scan(file,sep="\n",what="character")
			Ntree<-length(text)
			
		}
	}  else {
		trans<-NULL
		Ntree<-length(text)
	}
	tree<-lapply(text,text_to_tree,version,rev.order,trans)
	class(tree)<-c("multiSimmap","multiPhylo")
	if(length(tree)==1) tree<-tree[[1]]
	return(tree)
}

# function gets branch length
# written by Liam J. Revell 2011-2013

getBranch<-function(text,start,version){
	if(version<=1.0){
		i<-start+1
		l<-1
		if(text[i]=="{"){
			maps<-vector()
			i<-i+1
			while(text[i]!="}"){
				temp<-vector()
				m<-1
				while(text[i]!=","){
					temp[m]<-text[i]
					i<-i+1
					m<-m+1
				}
				state<-paste(temp,collapse="")
				i<-i+1
				temp<-vector()
				m<-1
				while(text[i]!=":"&&text[i]!="}"){
					temp[m]<-text[i]
					i<-i+1
					m<-m+1
				}
				length<-as.numeric(paste(temp,collapse=""))
				maps[l]<-length
				names(maps)[l]<-as.character(state)
				l<-l+1
				if(text[i]==":") i<-i+1
			}
		}
		edge.length<-sum(maps) # create branch length
		i<-i+1
	} else if(version>1.0){
		i<-start+1
		l<-1
		if(text[i]=="["){
			maps<-vector()
			i<-i+7
			while(text[i]!="}"){
				temp<-vector()
				m<-1
				while(text[i]!=","&&text[i]!="}"){
					temp[m]<-text[i]
					i<-i+1
					m<-m+1
				}
				state<-paste(temp,collapse="")
				if(text[i]!="}"){
					i<-i+1
					temp<-vector()
					m<-1
					while(text[i]!=","&&text[i]!="}"){
						temp[m]<-text[i]
						i<-i+1
						m<-m+1
					}
					ll<-as.numeric(paste(temp,collapse=""))
					maps[l]<-ll
					names(maps)[l]<-as.character(state)
					l<-l+1
					if(text[i]==",") i<-i+1
				}
			}
			temp<-vector()
			m<-1	
			i<-i+2
			while(is.na(match(text[i],c(",",")")))){
				temp[m]<-text[i]
				i<-i+1
				m<-m+1
			}
			maps[l]<-as.numeric(paste(temp,collapse=""))-sum(maps)
			names(maps)[l]<-as.character(state)
		}
		edge.length<-sum(maps)
	}
	return(list(maps=maps,edge.length=edge.length,end=i))
}

# make a mapped edge matrix
# written by Liam J. Revell 2013

makeMappedEdge<-function(edge,maps){
	st<-sort(unique(unlist(sapply(maps,function(x) names(x)))))
	mapped.edge<-matrix(0,nrow(edge),length(st))
	rownames(mapped.edge)<-apply(edge,1,function(x) paste(x,collapse=","))
	colnames(mapped.edge)<-st
	for(i in 1:length(maps)) 
		for(j in 1:length(maps[[i]])) 
			mapped.edge[i,names(maps[[i]])[j]]<-mapped.edge[i,names(maps[[i]])[j]]+maps[[i]][j]
	return(mapped.edge)
}


# function translates text-string to tree
# written by Liam J. Revell 2011-2015

text_to_tree<-function(text,version,rev.order,trans){
	text<-unlist(strsplit(text, NULL))
	tip.label<-vector(mode="character") 
	edge<-matrix(c(1,NA),1,2) 
	edge.length<-vector()
	maps<-list()
	currnode<-1
	Nnode<-currnode
	i<-j<-k<-1
	while(text[i]!=";"){
		if(text[i]=="("){
			if(j>nrow(edge)) edge<-rbind(edge,c(NA,NA))
			edge[j,1]<-currnode
			i<-i+1
			# is the next element a label?
			if(is.na(match(text[i],c("(",")",",",":",";")))){
				temp<-getLabel(text,i)
				tip.label[k]<-temp$label
				i<-temp$end
				edge[j,2]<--k
				k<-k+1
				# is there a branch length?
				if(text[i]==":"){
					temp<-getBranch(text,i,version)
					maps[[j]]<-temp$maps
					edge.length[j]<-temp$edge.length
					i<-temp$end
				}	
			} else if(text[i]=="("){
				Nnode<-Nnode+1 # creating a new internal node
				currnode<-Nnode
				edge[j,2]<-currnode # move to new internal node
			}
			j<-j+1
		} else if(text[i]==")"){
			i<-i+1
			# is there a branch length?
			if(text[i]==":"){
				temp<-getBranch(text,i,version)
				ii<-match(currnode,edge[,2])
				maps[[ii]]<-temp$maps
				edge.length[ii]<-temp$edge.length
				i<-temp$end
			}	
			currnode<-edge[match(currnode,edge[,2]),1] # move down the tree
		} else if(text[i]==","){
			if(j>nrow(edge)) edge<-rbind(edge,c(NA,NA))
			edge[j,1]<-currnode
			i<-i+1
			# is the next element a label?
			if(is.na(match(text[i],c("(",")",",",":",";")))){
				temp<-getLabel(text,i)
				tip.label[k]<-temp$label
				i<-temp$end
				edge[j,2]<--k
				k<-k+1
				# is there a branch length?
				if(text[i]==":"){
					temp<-getBranch(text,i,version)
					maps[[j]]<-temp$maps
					edge.length[j]<-temp$edge.length
					i<-temp$end
				}
			} else if(text[i]=="("){
				Nnode<-Nnode+1 # creating a new internal node
				currnode<-Nnode
				edge[j,2]<-currnode # move to internal node
			}
			j<-j+1
		}
	}
	Ntip<-k-1
	edge[edge>0]<-edge[edge>0]+Ntip
	edge[edge<0]<--edge[edge<0]
	mapped.edge<-makeMappedEdge(edge,maps)
	maps<-if(rev.order) lapply(maps,function(x) x<-x[length(x):1]) else maps
	if(!is.null(trans))tip.label<-trans[tip.label]
	# assemble into modified "phylo" object
	tree<-list(edge=edge,Nnode=as.integer(Nnode),tip.label=tip.label,edge.length=edge.length,maps=maps,mapped.edge=mapped.edge)
	class(tree)<-c("simmap","phylo")
	attr(tree,"map.order")<-if(rev.order) "right-to-left" else "left-to-right"
	return(tree)
}

# function reads Nexus taxa block in SIMMAP format trees
# modified from ape:read.nexus
# written by Liam J. Revell 2011-2013

readNexusData<-function(file,version){
	# read modified nexus block
	# this function is adapted from ape:read.nexus (Paradis et al. 2004)
	if(version<=1.0){
		X<-scan(file=file,what="",sep="\n",quiet=TRUE)
		left<-grep("\\[",X)
		right<-grep("\\]",X)
		if(length(left)){
			w<-left==right
			if(any(w)){
				s<-left[w]
				X[s]<-gsub("\\[[^]]*\\]","",X[s])
			}
			w<-!w
			if(any(w)){
			s<-left[w]
			X[s]<-gsub("\\[.*","",X[s])
			sb<-right[w]
			X[sb]<-gsub(".*\\]","",X[sb])
			if(any(s<sb-1))
				X<-X[-unlist(mapply(":",(s+1),(sb-1)))]
			}
		}
		endblock<-grep("END;|ENDBLOCK;",X,ignore.case=TRUE)
		semico<-grep(";",X)
		i1<-grep("begin smptrees;",X,ignore.case=TRUE)
		i2<-grep("translate",X,ignore.case=TRUE)
		translation<-if(length(i2)==1 && i2>i1) TRUE else FALSE
		if(translation){
			end<-semico[semico>i2][1]
			x<-X[(i2+1):end]
			x<-unlist(strsplit(x,"[,;\t]"))
			x<-unlist(strsplit(x," ")) # this is an addition
			x<-x[nzchar(x)]
			trans<-matrix(x,ncol=2,byrow=TRUE)
			trans[,2]<-gsub("['\"]","",trans[,2])
			n<-dim(trans)[1]
		}
		start<-if(translation) semico[semico>i2][1]+1 else semico[semico>i1][1]
		end<-endblock[endblock>i1][1]-1
		tree<-X[start:end]
		tree<-gsub("^.*= *","",tree)
		tree<-tree[tree!=""]
		semico<-grep(";",tree)
		Ntree<-length(semico)
		if(Ntree==1&&length(tree)>1){
			STRING<-paste(tree,collapse="")
		} else {
			if(any(diff(semico)!=1)){
				STRING<-character(Ntree)
				s<-c(1,semico[-Ntree]+1)
				j<-mapply(":",s,semico)
				if(is.list(j)){
					for(i in 1:Ntree) STRING[i]<-paste(tree[j[[i]]],collapse="")
				} else {
					for(i in 1:Ntree) STRING[i]<-paste(tree[j[,i]],collapse="")
				}
			} else STRING<-tree
		}
		text<-STRING
		if(translation==TRUE){
			rownames(trans)<-trans[,1]
			trans<-trans[,2]
			return(list(text=text,trans=trans,Ntree=Ntree))
		} else return(list(text=text,Ntree=Ntree))
	} else if(version>1.0){
		X<-scan(file=file,what="",sep="\n",quiet=TRUE)
		left<-grep("\\[",X)
		right<-grep("\\]",X)
		skip<-if(version<=2.0) grep("\\&map",X) else if(version>2.0) grep("\\&prob",X)
		left<-setdiff(left,skip)
		right<-setdiff(right,skip)
		if(length(left)){
			w<-left==right
			if(any(w)){
				s<-left[w]
				X[s]<-gsub("\\[[^]]*\\]","",X[s])
			}
			w<-!w
			if(any(w)){
			s<-left[w]
			X[s]<-gsub("\\[.*","",X[s])
			sb<-right[w]
			X[sb]<-gsub(".*\\]","",X[sb])
			if(any(s<sb-1))
				X<-X[-unlist(mapply(":",(s+1),(sb-1)))]
			}
		}
		endblock<-grep("END;|ENDBLOCK;",X,ignore.case=TRUE)
		semico<-grep(";",X)
		i1<-grep("begin trees;",X,ignore.case=TRUE)
		i2<-grep("translate",X,ignore.case=TRUE)
		translation<-if(length(i2)==1 && i2>i1) TRUE else FALSE
		if(translation){
			end<-semico[semico>i2][1]
			x<-X[(i2+1):end]
			x<-unlist(strsplit(x,"[,;\t]"))
			x<-unlist(strsplit(x," ")) # this is an addition
			x<-x[nzchar(x)]
			trans<-matrix(x,ncol=2,byrow=TRUE)
			trans[,2]<-gsub("['\"]","",trans[,2])
			n<-dim(trans)[1]
		}
		start<-if(translation) semico[semico>i2][1]+1 else semico[semico>i1][1]
		end<-endblock[endblock>i1][1]-1
		tree<-X[start:end]
		tree<-sub(".* = *","",tree)
		tree<-tree[tree!=""]
		semico<-grep(";",tree)
		Ntree<-length(semico)
		if(Ntree==1&&length(tree)>1){
			STRING<-paste(tree,collapse="")
		} else {
			if(any(diff(semico)!=1)){
				STRING<-character(Ntree)
				s<-c(1,semico[-Ntree]+1)
				j<-mapply(":",s,semico)
				if(is.list(j)){
					for(i in 1:Ntree) STRING[i]<-paste(tree[j[[i]]],collapse="")
				} else {
					for(i in 1:Ntree) STRING[i]<-paste(tree[j[,i]],collapse="")
				}
			} else STRING<-tree
		}
		text<-STRING
		if(translation==TRUE){
			rownames(trans)<-trans[,1]
			trans<-trans[,2]
			return(list(text=text,trans=trans,Ntree=Ntree))
		} else return(list(text=text,Ntree=Ntree))
	}
}
