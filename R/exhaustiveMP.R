# function does branch & bound or exhaustive MP tree search
# uses "phangorn (Schliep 2011) & "ape" (Paradis & Strimmer 2004)
# data is a phyDat object; method can be "branch.and.bound" or "exhaustive"
# written by Liam J. Revell 2011, 2013

exhaustiveMP<-function(data,tree=NULL,method="branch.and.bound"){
	if(method=="branch.and.bound"){
		if(length(data)>15) stop("branch and bound only allowed for n<=15")
		if(is.null(tree)){
			if(attr(data,"type")=="DNA"){
				print("no input tree; starting with NJ tree")
				tree<-NJ(dist.dna(as.DNAbin(data)))
			} else {
				print("no input tree; using random starting tree")
				tree<-rtree(n=length(data),tip.label=names(data),br=NULL,rooted=FALSE)
			}
		}
		trees<-branch.and.bound(data,tree)
	} else if(method=="exhaustive"){
		if(length(data)>10) stop("exhaustive search only allowed for n<=10")
		if(!is.null(tree)) print("starting tree not necessary for exhaustive search")
		trees<-exhaustive.search(data)
	}
	if(length(trees)==1) trees<-trees[[1]]
	return(trees)
}

# function performs branch & bound search
# written by Liam J. Revell 2011

branch.and.bound<-function(data,tree){
	# first, compute the parsimony score on the input tree
	bound<-parsimony(tree,data)
	# now pick three species to start
################## this is new and I don't know if it helps
	if(is.null(tree$edge.length)){
		print("starting with 3 species chosen at random")
		new<-list(stree(n=3,tip.label=sample(tree$tip.label,3)))
	} else {
		print("starting with 3 species chosen to maximize distance") 
		mdSp<-names(sort(colSums(cophenetic(tree)),decreasing=TRUE))[1:3]
		mdSp<-c("Tarsier","Chimp","J.Macaque")
		print(mdSp)
		new<-list(stree(n=3,tip.label=mdSp)) 
	}
################## ends here
	class(new)<-"multiPhylo"
	added<-new[[1]]$tip.label; remaining<-setdiff(tree$tip.label,added)
	# branch & bound
	while(length(remaining)>0){
		old<-new; new<-list()
		new.tip<-sample(remaining,1)
		pscores<-vector()
		for(i in 1:length(old)){			
			temp<-add.everywhere(old[[i]],new.tip)
			score<-parsimony(temp,data)
			new<-c(new,temp[score<=bound])
			pscores<-c(pscores,score[score<=bound])
		}
		added<-c(added,new.tip)
		print(paste(length(added),"species added;",length(new),"trees retained",collapse=""))
		remaining<-setdiff(tree$tip.label,added)
	}
	# ok, done, now sort what needs to be returned
	trees<-new[pscores==min(pscores)]
	for(i in 1:length(trees)) attr(trees[[i]],"pscore")<-min(pscores)
	return(trees) # return all mp trees
}

# function does exhaustive tree search
# written by Liam J. Revell 2011

exhaustive.search<-function(data){
	all.trees<-allTrees(n=length(data),tip.label=names(data),rooted=FALSE)
	print(paste("searching",length(all.trees),"trees",collapse=""))
	all.trees = .uncompressTipLabel(all.trees)
	pscores<-parsimony(all.trees,data)
	minscore<-min(pscores); trees<-all.trees[pscores==minscore]
	for(i in 1:length(trees)) attr(trees[[i]],"pscore")<-min(pscores)
	return(trees)
} 
