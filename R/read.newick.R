# function to read a Newick string with node labels & (possible) singles
# written by Liam J. Revell 2013, 2014, 2015

read.newick<-function(file="",text,...){
	# check to see if reading from file
	if(file!="") text<-scan(file,sep="\n",what="character",...)
	if(length(text)>1){
		tree<-lapply(text,newick)
		class(tree)<-"multiPhylo"
	} else tree<-newick(text)
	return(tree)
}

# main Newick string function
# written by Liam J. Revell 2013, 2014
newick<-function(text){
	text<-unlist(strsplit(text, NULL))
	tip.label<-vector(mode="character")
	node.label<-vector(mode="character") 
	edge<-matrix(NA,sum(text=="(")+sum(text==","),2)
	ei<-vector()
	edge.length<-vector()
	currnode<-1
	Nnode<-currnode
	i<-j<-k<-1 
	while(text[i]!=";"){
		if(text[i]=="("){
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
				ei[currnode]<-j
			}
			j<-j+1
		} else if(text[i]==")"){
			i<-i+1
			# is the next element a label?
			if(is.na(match(text[i],c("(",")",",",":",";")))){
				temp<-getLabel(text,i)
				node.label[currnode]<-temp$label
				i<-temp$end
			}
			ii<-ei[currnode]
			# is there a branch length?
			if(text[i]==":"){
				temp<-getEdgeLength(text,i)
				if(currnode>1) edge.length[ii]<-temp$edge.length
				else root.edge<-temp$edge.length
				i<-temp$end
			}	
			if(currnode>1) currnode<-edge[ii,1] # move down the tree	`
		} else if(text[i]==","){
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
				ei[currnode]<-j
			}
			j<-j+1
		}
	}
	Ntip<-k-1
	edge<-edge[!is.na(edge[,2]),]
	edge[edge>0]<-edge[edge>0]+Ntip
	edge[edge<0]<--edge[edge<0]
	edge.length[is.na(edge.length)]<-0
	if(length(edge.length)==0) edge.length<-NULL
	node.label[is.na(node.label)]<-""
	if(length(node.label)==0) node.label<-NULL
	# assemble into "phylo" object
	tree<-list(edge=edge,Nnode=as.integer(Nnode),tip.label=tip.label,edge.length=edge.length,node.label=node.label)
	class(tree)<-"phylo"
	attr(tree,"order")<-"cladewise"
	return(tree)
}

# function gets label
# written by Liam J. Revell 2011-2013, 2014
getLabel<-function(text,start,stop.char=c(",",":",")",";")){
	i<-0
	while(is.na(match(text[i+start],stop.char))) i<-i+1
	label<-paste(text[0:(i-1)+start],collapse="")	
	return(list(label=paste(label,collapse=""),end=i+start))
}

# function gets branch length
# written by Liam J. Revell 2011-2013, 2014
getEdgeLength<-function(text,start){
	i<-start+1
	stop.char<-c(",",")",";")
	while(is.na(match(text[i],stop.char))) i<-i+1
	edge.length<-as.numeric(paste(text[(start+1):(i-1)],collapse=""))
	return(list(edge.length=edge.length,end=i))
}
