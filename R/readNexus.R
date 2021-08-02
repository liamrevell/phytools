## function
## written by Liam J. Revell 2013, 2021

readNexus<-function(file="",format=c("standard","raxml")){
	format<-format[1]
	if(tolower(format)=="standard") tree<-read.nexus(file)
	else if(tolower(format)=="raxml"){
		XX<-readNexusData(file,version=3.5)
		text<-XX$text
		trans<-XX$trans
		Ntree<-XX$Ntree
		tree<-lapply(text,modified.text_to_tree,trans=trans)
		if(length(tree)==1) tree<-tree[[1]] else 
			class(tree)<-"multiPhylo"
	}
	else { 
		cat("Do not recognize format\n")
		tree<-NULL
	}
	tree
}

# function gets edge length
# written by Liam J. Revell 2013, 2021

getEdgeLength<-function(text,start){
	i<-start+1
	l<-1
	temp<-vector()
	while(is.na(match(text[i],c(",",")")))){
		temp[l]<-text[i]; l<-l+1; i<-i+1
	}
	list(edge.length=as.numeric(paste(temp,collapse="")),
		end=i)
}

## function gets bootstrap stored as [&label="bootstrap%"]
## written by Liam J. Revell 2021

getBS<-function(text,start){
	i<-start
	if(text[i]=="["){
		j<-1
		xx<-vector()
		while(text[i]!="]"){
			i<-i+1
			xx[j]<-text[i]
			j<-j+1
		}
		label<-paste(xx[1:(length(xx)-1)],collapse="")
		label<-sub("&label=","",label)
	}
	list(label=label,end=i+1)
}

## function translates text-string to tree
## written by Liam J. Revell 2011-2013, 2021

modified.text_to_tree<-function(text,trans){
	text<-unlist(strsplit(text,NULL))
	tip.label<-vector(mode="character") 
	edge<-matrix(c(1,NA),1,2)
	edge.length<-vector()
	node.label<-vector(mode="character")
	currnode<-1
	Nnode<-currnode
	i<-j<-k<-1
	while(text[i]!="(") i<-i+1
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
					temp<-getEdgeLength(text,i)
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
			# is there a node label?
			if(text[i]=="["){
				temp<-getBS(text,i)
				node.label[currnode]<-as.character(temp$label)
				i<-temp$end
			}
			# is there a branch length?
			if(text[i]==":"){
				temp<-getEdgeLength(text,i)
				ii<-match(currnode,edge[,2])
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
					temp<-getEdgeLength(text,i)
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
	if(!is.null(trans)) tip.label<-trans[tip.label]
	# assemble into modified "phylo" object
	ntip<-abs(min(edge))
	edge[which(edge>0)]<-ntip+edge[which(edge>0)]
	edge[which(edge<0)]<-abs(edge[which(edge<0)])
	tree<-list(edge=edge,Nnode=as.integer(Nnode),tip.label=tip.label,
		edge.length=edge.length,node.label=node.label)
	class(tree)<-"phylo"
	tree
}
